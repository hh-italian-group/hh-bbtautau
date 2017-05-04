/*! Study of correlation matrix and mutual information of BDT variables
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/ClassInfo.h"
#include "TMVA/VariableTransformBase.h"
#include <TMatrixDEigen.h>
#include <TH1.h>
#include <TH2.h>
#include <fstream>
#include <random>

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(Long64_t, number_events);
};

namespace analysis {

bool SM=false;

struct SampleEntry{
  std::string filename;
  double weight;
  bool issignal;
  SampleEntry() : weight(-1), issignal(false){}
};

std::ostream& operator<<(std::ostream& os, const SampleEntry& entry)
{
    os << entry.filename << " " << entry.weight << " " << std::boolalpha << entry.issignal;
    return os;
}

std::istream& operator>>(std::istream& is, SampleEntry& entry)
{
    is >> entry.filename >> entry.weight >> std::boolalpha >> entry.issignal;
    return is;
}

class MvaVariables{
private:
    std::vector<double> variables;
    std::map<std::string, size_t> names;
    std::shared_ptr<TMVA::DataLoader> loader;
public:
    MvaVariables(std::shared_ptr<TMVA::DataLoader> dloader) : loader(dloader){}

    double& operator[](const std::string& name)
    {
      if (!names.count(name)){
          variables.push_back(0);
          names[name]=variables.size()-1;
          loader->AddVariable(name);
      }
      return variables.at(names.at(name));
    }

    void AddEvent(bool issignal, bool istraining, double weight)
    {
        const std::string samplename = issignal ? "Signal" : "Background";
        const TMVA::Types::ETreeType treetype= istraining ? TMVA::Types::kTraining : TMVA::Types::kTesting;
        loader->AddEvent(samplename, treetype, variables, weight);
    }


};

class MvaData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    TH2D_ENTRY(histo, 100, 0, 0, 100, 0, 0)
    TH2D_ENTRY(mutualinformation, 100, 0, 0, 100, 0, 0)
};

//Create covariance matrix from correlation matrix and standar deviations.
TMatrixD CreateCovarianceMatrix(const TMatrixD corr_matrix, std::vector<double>& std_dev_vector, int N)
{
    TMatrixD cov_matrix(N,N);
    for(int n = 0; n < N; n++) {
        if(corr_matrix(n,n) != 1)
            throw std::runtime_error("Self-correlation not equal to one.");
        for(int k = n; k < N; k++) {
            if(std::abs(corr_matrix(n,k)) > 1)
                throw std::runtime_error("Correlation coefficient is greater than one.");
            if(corr_matrix(n,k) != corr_matrix(k,n))
                throw std::runtime_error("Correlation matrix is not symmetric.");
            cov_matrix(n,k) = corr_matrix(n,k) * std_dev_vector[n] * std_dev_vector[k];
            cov_matrix(k,n) = cov_matrix(n,k);
        }
    }
    return cov_matrix;


    auto sample_vars = mvaVars.GetSampleVariables("signal");

    using VarPair = std::pair<std::string, std::string>;
    std::map<VarPair, double> cov_matrix;

    for(const auto& var_1 : sample_vars) {
        for(const auto& var_2 : sample_vars) {
            const VarPair var_12(var_1.first, var_2.first);
            const VarPair var_21(var_2.first, var_1.first);
            if(/*corr_matrix.count(var_12) ||*/ corr_matrix.count(var_21)) continue;

            const double cov = Covariance(var_1.second, var_2.second);
            corr_matrix[var_12] = corr;
        }
    }
    
    auto corr_matrix = CovToCorr(cov_matrix);
    auto corr_matrix_fix = RemoveDiag(corr_matrix);
    std::vector< std::pair<VarPair, double> > corr_vector(corr_matrix_fix.begin(), corr_matrix_fix.end());

    std::sort(corr_vector.begin(), corr_vector.end(), [](const std::pair<VarPair, double>& a, const std::pair<VarPair, double>& b) {
        return std::abs(a.second) > std::abs(b.second);
    });
}

//Create variance matrix from standard deviation vector
TMatrixD CreateVarianceMatrix(std::vector<double> std_dev_vector, int N)
{
    TMatrixD var_matrix(N, N);
    for(int n = 0; n < N; n++)
        var_matrix(n,n) = std::pow(std_dev_vector[n], 2);
    return var_matrix;
}

TMatrixD DiagonalizeLinearTransformMatrix(const TMatrixD& matrix, TMatrixD* eigen = 0, TMatrixD* eigen_inv = 0)
{
    if(matrix.GetNcols() != matrix.GetNrows())
        throw std::runtime_error("Matrix is not square.");

    TMatrixDEigen eigen_producer(matrix);
    auto _eigen = eigen_producer.GetEigenVectors();
    auto _eigen_inv = _eigen;
    _eigen_inv.Invert();
    if(eigen) *eigen = _eigen;
    if(eigen_inv) *eigen_inv = _eigen_inv;
    return _eigen_inv * matrix * _eigen;
}

template<typename T>
TMatrixD MatrixPower(const TMatrixD& matrix, T power)
{
    TMatrixD eigen(matrix.GetNrows(), matrix.GetNcols()), eigen_inv = eigen;
    auto diag = DiagonalizeLinearTransformMatrix(matrix, &eigen, &eigen_inv);
    for(int n = 0; n < diag.GetNcols(); ++n)
        diag[n][n] = std::pow(diag[n][n], power);
    return eigen * diag * eigen_inv;
}

// Compute matrix W using optimal whitening to decorrelat original variables.
// Variables z=Wx will be decorrelated between each other.
// Used whitening formula defined by A. Kessy, A. Lewin and K. Strimmer in https://arxiv.org/abs/1512.00809.
inline TMatrixD ComputeWhiteningMatrix(const TMatrixD& corr_matrix, std::vector<double> std_dev_vector, int N)
{
    if(corr_matrix.GetNcols() != N)
        throw std::runtime_error("Correlation matrix is not square.");
    if(std_dev_vector.size() != N)
        throw std::runtime_error("Size of the std dev vector do not correspond to the size of the correlation matrix.");

    auto corr_inv_sqrt = MatrixPower(corr_matrix, -0.5);
    TMatrixD var_matrix = CreateVarianceMatrix(std_dev_vector,N);
    auto var_inv_sqrt = MatrixPower(var_matrix, -0.5);
    return corr_inv_sqrt * var_inv_sqrt;
}



//Mutual Information method for non-linear correlations estimates in 2D histogram
double GetMutualInformation(TH2D h){
    double area=h.Integral();
    TH2D hc(h); //copy histogram and rebin to speed up procedure
    hc.RebinX(2);
    hc.RebinY(2);
    int maxbin_x=hc.GetNbinsX();
    int maxbin_y=hc.GetNbinsY();
    double mutualinformation=0;

    for (Int_t x = 1; x <= maxbin_x; x++) {
        for (Int_t y = 1; y <= maxbin_y; y++) {
            double p_xy = hc.GetBinContent(x,y)/area;
            double p_x  = hc.Integral(x,x,1,maxbin_y)/area;
            double p_y  = hc.Integral(1,maxbin_x,y,y)/area;
            if (p_x > 0. && p_y > 0. && p_xy > 0.){
                mutualinformation += p_xy*TMath::Log(p_xy / (p_x * p_y));
            }
        }
    }
    return mutualinformation;
}



//Fill the matrix of mutual information
void GetMatrixMutualInformation(TMatrixD &mutualinformation, MvaData &anaData, std::shared_ptr<TMVA::DataLoader> dataloader, int nvar,std::string className, std::string info1,std::string info2){
    for(Long64_t i = 0; i < nvar; i++){
        for(Long64_t j = i; j < nvar; j++){
            info1=(dataloader->GetDataSetInfo().GetVariableInfo(i).GetTitle());
            info2=(dataloader->GetDataSetInfo().GetVariableInfo(j).GetTitle());
            mutualinformation(i,j)=GetMutualInformation(anaData.histo(className+info1+" "+info2));
            mutualinformation(j,i)=mutualinformation(i,j);
//            std::cout<<className<<i<<" "<<j<<" "<<mutualinformation(i,j)<<std::endl;
        }
    }
}

//Fill a 2D-histo from the  mutual information matrix
void GetHistoMutualInformation(TMatrixD mutualinformation, MvaData &anaData, std::shared_ptr<TMVA::DataLoader> dataloader, int nvar,  std::string className){
    anaData.mutualinformation(className).SetCanExtend(TH1::kAllAxes);
    anaData.mutualinformation(className).SetBins(nvar, 0, nvar, nvar, 0, nvar);
    for (Long64_t ibin=1; ibin<=nvar; ibin++) {
        anaData.mutualinformation(className).GetXaxis()->SetBinLabel(ibin, dataloader->GetDataSetInfo().GetVariableInfo(ibin-1).GetTitle());
        anaData.mutualinformation(className).GetYaxis()->SetBinLabel(ibin, dataloader->GetDataSetInfo().GetVariableInfo(ibin-1).GetTitle());
        for (Long64_t jbin=1; jbin<=nvar; jbin++) {
            anaData.mutualinformation(className).SetBinContent( ibin, jbin, mutualinformation(ibin-1,jbin-1));
        }
    }
}

//Compute the difference of two matrix
TMatrixD DifferenceMatrix(TMatrixD m1,TMatrixD m2, int dim)
{
    TMatrixD difference(dim,dim);
    for (Long64_t i=0; i<dim; i++) {
        for (Long64_t j=0; j<dim; j++) {
            difference(i,j)=m1(i,j)-m2(i,j);
            difference(j,i)=difference(i,j);
        }
    }
    return difference;
}

std::vector<std::pair<double, std::string>> GetElements(TMatrixD matrix, std::shared_ptr<TMVA::DataLoader> dataloader,int nvar){

    std::string info1,info2;
    std::vector<std::pair<double, std::string>> elements;
    for(Long64_t i = 0; i < nvar; i++){
        info1=(dataloader->GetDataSetInfo().GetVariableInfo(i).GetTitle());
        for(Long64_t j = i; j < nvar; j++){
            info2=(dataloader->GetDataSetInfo().GetVariableInfo(j).GetTitle());
            elements.emplace_back(matrix(i,j),info1+" "+info2);
        }
    }
    return elements;
}

//Create plots of Jensen Shannon between signal and background for different values of mass
void KLDPlotSignalBkg(std::map<int, std::map<std::string,std::vector<double>>> sample_signal, std::map<int, std::map<std::string,std::vector<double>>> sample_bkg, std::map<int, std::map<std::string, double>> band_signal, std::map<int, std::map<std::string, double>> band_bkg, std::shared_ptr<TFile> outfile){
    std::map<int,std::vector<std::pair<std::string,double>>> JSDivergence;
    std::map<std::string, TGraph*> plot;
    auto map_bkg = sample_bkg[1];
    auto bandwidth_bkg = band_bkg[1];
    int i = 0;
    for(const auto& var_mass_2 : sample_signal){
        if (var_mass_2.second.size() == 0) continue;
        auto map_signal_2 = var_mass_2.second;
        auto bandwidth_signal = band_signal[var_mass_2.first];
        auto JSDivergence_mass = JSDivergence[var_mass_2.first];
        for (const auto& var : var_mass_2.second){
            double k = stat_estimators::JensenShannonDivergence(map_bkg[var.first], map_signal_2[var.first], bandwidth_bkg[var.first],bandwidth_signal[var.first]);
            JSDivergence_mass.emplace_back(var.first,k);
            if (plot.count(var.first)==0) plot[var.first] = new TGraph();
            plot[var.first]->SetPoint(i,var_mass_2.first, k);
            plot[var.first]->SetLineColor(kGreen+1);
            plot[var.first]->SetLineWidth(1);
            plot[var.first]->SetMarkerColor(1);
            plot[var.first]->SetMarkerSize(1);
            plot[var.first]->SetMarkerStyle(3);
            plot[var.first]->SetTitle(("JSD_"+var.first+"_SignalBkg").c_str());
            plot[var.first]->SetName(("JSD_"+var.first+"_SignalBkg").c_str());
            plot[var.first]->GetHistogram()->GetXaxis()->SetTitle("mass");
            plot[var.first]->GetHistogram()->GetYaxis()->SetTitle("JS Divergence");
        }
        std::sort(JSDivergence_mass.begin(),JSDivergence_mass.end(), [](std::pair<std::string,double> el1, std::pair<std::string,double> el2){
            return el1.second < el2.second;
        });
        JSDivergence[var_mass_2.first] = JSDivergence_mass;
        i++;
    }

    std::ofstream ofs("JSDivergence_signalvsbkg.csv", std::ofstream::out);
    for(const auto& var_mass_2 : sample_signal){
        ofs<<"Massa: "<<","<<var_mass_2.first<<std::endl;
        auto JSDivergence_mass = JSDivergence[var_mass_2.first];
        for (Long64_t j = 0; j < map_bkg.size(); j++){
            ofs<<JSDivergence_mass[j].first<<",";
        }
        ofs<<std::endl;
        for (Long64_t j = 0; j < map_bkg.size(); j++){
            ofs<<JSDivergence_mass[j].second<<",";
        }
        ofs<<std::endl;
    }
    ofs.close();

    for(const auto& var: map_bkg){
        root_ext::WriteObject(*plot[var.first], outfile.get());
    }

}
bool PairCompare (const std::pair<double, std::string> &firstelement, const std::pair<double, std::string> &secondelement){
    return firstelement.first > secondelement.first;
}

void ReadVector(std::vector<std::pair<double, std::string>> vector, Arguments args, std::string classname, std::string type){
    if (SM==true) {
        std::ofstream ofs ("SM"+args.tree_name()+classname+type+".txt", std::ofstream::out);
        for(Long64_t i = 0; i < vector.size(); i++){
            if (vector[i].first!=1)
                ofs <<i<<"    "<<vector[i].second<<"    "<<vector[i].first<<std::endl;
        }
        ofs.close();
    }
    else {
        std::ofstream ofs ("lm"+args.tree_name()+classname+type+".txt", std::ofstream::out);
        for(Long64_t i = 0; i < vector.size(); i++){
            if (vector[i].first!=1)
                ofs <<i<<"    "<<vector[i].second<<"    "<<vector[i].first<<std::endl;
        }
        ofs.close();
    }

}


using SampleEntryCollection = std::vector<SampleEntry>;

class MvaClassification {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    static const std::set<std::string>& GetDisabledBranches()
    {
        static const std::set<std::string> DisabledBranches_read = {
            "dphi_mumet", "dphi_metsv", "dR_taumu", "mT1", "mT2", "dphi_bbmet", "dphi_bbsv", "dR_bb", "m_bb", "n_jets",
            "btag_weight", "ttbar_weight",  "PU_weight", "shape_denominator_weight", "trigger_accepts", "trigger_matches",
            "event.tauId_keys_1","event.tauId_keys_2","event.tauId_values_1","event.tauId_values_2"
        };
        return DisabledBranches_read;
    }


    static SampleEntryCollection ReadConfig(const std::string& cfg_file){
        std::ifstream f(cfg_file);
        SampleEntryCollection collection;
        while(f.good()){
            std::string line;
            std::getline(f, line);
            if (line.size()==0 || line.at(0)=='#') continue;
            std::istringstream s(line);
            SampleEntry entry;
            s>>entry;
            collection.push_back(entry);
        }
        return collection;
    }

    MvaClassification(const Arguments& _args): args(_args), samples(ReadConfig(args.cfg_file())),
        outfile(root_ext::CreateRootFile(args.output_file())),
        factory(new TMVA::Factory ("myFactory", outfile.get(),"!V:!Silent:Color:DrawProgressBar:Transformations=I,D:AnalysisType=Classification")),
        dataloader(new TMVA::DataLoader ("mydataloader")), vars(dataloader),  anaData(),
        gen(rd()), testvstraining(0,1)

    {
    }

    void Run()
    {
        std::string info1,info2;
        for(const SampleEntry& entry:samples)
        {
            auto input_file=root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
            EventTuple tuple(args.tree_name(), input_file.get(), true, GetDisabledBranches());
            std::cout<<entry<<" number of events: "<<std::min(tuple.GetEntries(),args.number_events())<<std::endl;
            std::string samplename = entry.issignal ? "Signal" : "Background";

            for(Long64_t current_entry = 0; current_entry < std::min(tuple.GetEntries(),args.number_events()); ++current_entry) {
                tuple.GetEntry(current_entry);
                const Event& event = tuple.data();
                const bool istraining= testvstraining(gen) == 0;

                if (event.eventEnergyScale!=0 || (event.q_1+event.q_2)!=0 || event.jets_p4.size() < 2
                    || event.extraelec_veto==true || event.extramuon_veto==true) continue;

                LorentzVectorE_Float bb= event.jets_p4[0] + event.jets_p4[1];
                LorentzVectorM_Float leptons= event.p4_1 + event.p4_2;
                LorentzVectorM_Float leptonsMET= event.p4_1 + event.p4_2 + event.pfMET_p4;

                double circular_cut=std::sqrt(pow(event.SVfit_p4.mass()-116.,2)+pow(bb.M()-111,2));
                if (circular_cut>40) continue;
//                if ((args.tree_name=="eTau" && circular_cut>40)||(args.tree_name=="muTau" && circular_cut>30)) continue;

                vars["pt_l1"]=event.p4_1.pt();
                vars["pt_l2"]=event.p4_2.pt();
                vars["pt_jet1"]=event.jets_p4[1].pt();
                vars["pt_jet2"]=event.jets_p4[2].pt();
                vars["pt_leptons"]=leptons.pt();
                vars["pt_sv"]=event.SVfit_p4.pt();
                vars["pt_leptonsmet"]=leptonsMET.pt();
                vars["pt_bb"]=bb.pt();
                vars["pt_met"]=event.pfMET_p4.pt();

//                vars["p_zeta"]=??

                vars["abs(dphi_leptons)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.p4_2));
                vars["dphi_leptons"] = (ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.p4_2));
                vars["abs(dphi_jets)"]=std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.jets_p4[0], event.jets_p4[1]));
                vars["dphi_jets"]=(ROOT::Math::VectorUtil::DeltaPhi(event.jets_p4[0], event.jets_p4[1]));
                vars["abs(dphi_l1met)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.pfMET_p4));
                vars["dphi_l1met"] = (ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.pfMET_p4));
                vars["abs(dphi_l2met)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_2, event.pfMET_p4));
                vars["dphi_l2met"] = (ROOT::Math::VectorUtil::DeltaPhi(event.p4_2, event.pfMET_p4));
                vars["abs(dphi_svmet)"]=std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4));
                vars["dphi_svmet"]=(ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4));
                vars["abs(dphi_bbmet)"]=std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.pfMET_p4));
                vars["dphi_bbmet"]=(ROOT::Math::VectorUtil::DeltaPhi(bb, event.pfMET_p4));
                vars["abs(dphi_bbsv)"]=std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.SVfit_p4));
                vars["dphi_bbsv"]=(ROOT::Math::VectorUtil::DeltaPhi(bb, event.SVfit_p4));

                vars["abs(deta_leptons)"] =std::abs(event.p4_1.eta() - event.p4_2.eta());
                vars["deta_leptons"] =(event.p4_1.eta() - event.p4_2.eta());
                vars["abs(deta_jets)"]=std::abs(event.jets_p4[0].eta() - event.jets_p4[1].eta());
                vars["deta_jets"]=(event.jets_p4[0].eta() - event.jets_p4[1].eta());
                vars["abs(deta_l1met)"] = std::abs(event.p4_1.eta()-event.pfMET_p4.eta());
                vars["deta_l1met"] = (event.p4_1.eta()-event.pfMET_p4.eta());
                vars["abs(deta_l2met)"] = std::abs(event.p4_2.eta()-event.pfMET_p4.eta());
                vars["deta_l2met"] = (event.p4_2.eta()-event.pfMET_p4.eta());
                vars["abs(deta_svmet)"]=std::abs(event.SVfit_p4.eta()-event.pfMET_p4.eta());
                vars["deta_svmet"]=(event.SVfit_p4.eta()-event.pfMET_p4.eta());
                vars["abs(deta_bbmet)"]=std::abs(bb.eta()-event.pfMET_p4.eta());
                vars["deta_bbmet"]=(bb.eta()-event.pfMET_p4.eta());
                vars["abs(deta_bbsv)"]=std::abs(bb.eta()-event.SVfit_p4.eta());
                vars["deta_bbsv"]=bb.eta()-event.SVfit_p4.eta();

                vars["dR_leptons"] = std::abs(ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2));
                vars["dR_jets"]=std::abs(ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]));
                vars["dR_l1met"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.pfMET_p4));
                vars["dR_l2met"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_2, event.pfMET_p4));
                vars["dR_svmet"]=(ROOT::Math::VectorUtil::DeltaR(event.SVfit_p4, event.pfMET_p4));
                vars["dR_bbmet"]=(ROOT::Math::VectorUtil::DeltaR(bb, event.pfMET_p4));
                vars["dR_bbsv"]=(ROOT::Math::VectorUtil::DeltaR(bb, event.SVfit_p4));

                vars["dR_jetsPt_bb"]=std::abs(ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]))*bb.Pt();
                vars["dR_leptonsPt_sv"] = std::abs(ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2))*event.SVfit_p4.Pt();

                vars["mass_leptonsmet"]=leptonsMET.M();//
                vars["mass_leptonsmetinvariant"]=ROOT::Math::VectorUtil::InvariantMass(leptons,event.pfMET_p4); //
                vars["mass_sv"]=event.SVfit_p4.M();
                vars["mass_jets"]=bb.M();
                vars["MT_l1"]=Calculate_MT(event.p4_1,event.pfMET_p4);
                vars["MT_l2"]=Calculate_MT(event.p4_2,event.pfMET_p4);
                vars["MT_tautau"]=std::sqrt(pow(event.p4_1.Et()+event.p4_1.Et(),2)+pow(event.p4_1.px()+event.p4_1.px(),2)+pow(event.p4_1.py()+event.p4_1.py(),2));//
                vars["MT_tautausv"]=Calculate_MT(event.SVfit_p4, event.pfMET_p4); //
                vars["MT_leptons"]=Calculate_MT(leptons, event.pfMET_p4);
                vars["mass_H"]=ROOT::Math::VectorUtil::InvariantMass(bb,event.SVfit_p4);

                vars["theta_1"]=ROOT::Math::VectorUtil::Angle(event.p4_1,event.SVfit_p4);
                vars["theta_2"]=ROOT::Math::VectorUtil::Angle(event.jets_p4[1],bb);

                vars.AddEvent(entry.issignal, istraining, entry.weight);

                for(Long64_t i = 0; i < vars.variables.size(); i++){
                    info1=(dataloader->GetDataSetInfo().GetVariableInfo(i).GetTitle());
                    for(Long64_t j = i; j < vars.variables.size(); j++){
                        info2=(dataloader->GetDataSetInfo().GetVariableInfo(j).GetTitle());
                        anaData.histo(samplename+info1+" "+info2).Fill(vars.variables[i], vars.variables[j]);
                        anaData.histo(samplename+info1+" "+info2).SetCanExtend(TH1::kAllAxes);
                        anaData.histo(samplename+info1+" "+info2).SetYTitle(dataloader->GetDataSetInfo().GetVariableInfo(j).GetTitle());
                        anaData.histo(samplename+info1+" "+info2).SetXTitle(dataloader->GetDataSetInfo().GetVariableInfo(i).GetTitle());
                    }
                }
            }
        }

        dataloader->PrepareTrainingAndTestTree( "","", "SplitMode=Random" );
        dataloader->GetDefaultDataSetInfo().GetDataSet();
        int nvar=vars.variables.size();

        const TMatrixD* matrixcorrelation_signal=dataloader->GetDataSetInfo().CorrelationMatrix("Signal");
        const TH2* histocorrelation_signal(0);
        histocorrelation_signal=dataloader->GetDataSetInfo().CreateCorrelationMatrixHist(matrixcorrelation_signal,"CorrelationMatrixSignal", "Correlation Matrix (Signal)");
        root_ext::WriteObject(*histocorrelation_signal, outfile.get());
        std::vector<std::pair<double, std::string>> elements_correlation_signal;
        elements_correlation_signal=GetElements(*matrixcorrelation_signal,dataloader,nvar);
        std::sort(elements_correlation_signal.begin(),elements_correlation_signal.end(), PairCompare);
        ReadVector(elements_correlation_signal,args,"Signal","correlation");

        const TMatrixD* matrixcorrelation_background=dataloader->GetDataSetInfo().CorrelationMatrix("Background");
        const TH2* histocorrelation_background(0);
        histocorrelation_background=dataloader->GetDataSetInfo().CreateCorrelationMatrixHist(matrixcorrelation_background,"CorrelationMatrixBackground", "Correlation Matrix (Background)");
        root_ext::WriteObject(*histocorrelation_background, outfile.get());
        std::vector<std::pair<double, std::string>> elements_correlation_background;
        elements_correlation_background=GetElements(*matrixcorrelation_background,dataloader,nvar);
        std::sort(elements_correlation_background.begin(),elements_correlation_background.end(),PairCompare);
        ReadVector(elements_correlation_background,args,"Background","correlation");

        TMatrixD matrixcorrelation_difference(nvar,nvar);
        matrixcorrelation_difference=DifferenceMatrix(*matrixcorrelation_signal,*matrixcorrelation_background,nvar);
        const TH2* histocorrelation_difference(0);
        histocorrelation_difference=dataloader->GetDataSetInfo().CreateCorrelationMatrixHist(&matrixcorrelation_difference,"CorrelationMatrixDifference","Correlation Matrix (Difference)");
        root_ext::WriteObject(*histocorrelation_difference, outfile.get());
        std::vector<std::pair<double, std::string>> elements_correlation_difference;
        elements_correlation_difference=GetElements(matrixcorrelation_difference,dataloader,nvar);
        std::sort(elements_correlation_difference.begin(),elements_correlation_difference.end(),PairCompare);
        ReadVector(elements_correlation_difference, args,"Difference","correlation");

//        std::vector<double> std_dev_signal, std_dev_background;
//        for(Long64_t i = 0; i < vars.variables.size(); i++){
//            info1=(dataloader->GetDataSetInfo().GetVariableInfo(i).GetTitle());
//            std_dev_signal.push_back(anaData.histo("Signal"+info1+" "+info1).GetStdDev(1));
//            std_dev_background.push_back(anaData.histo("Background"+info1+" "+info1).GetStdDev(1));
//        }
//        TMatrixD matrixcovariance_signal(nvar,nvar);
//        matrixcovariance_signal=CreateCovarianceMatrix(*matrixcorrelation_signal,std_dev_signal, nvar);
//        whiteningmatrix_signal=ComputeWhiteningMatrix(*matrixcorrelation_signal,std_dev_signal, nvar);
//        std::cout<<"Correlation matrix (Signal)"<<std::endl;
//        matrixcorrelation_signal->Print();
//        std::cout<<"Covariance matrix (Signal)"<<std::endl;
//        matrixcovariance_signal.Print();
//        std::cout<<"Whitening matrix (Signal)"<<std::endl;
//        whiteningmatrix_signal.Print();
//        TMatrixD matrixcovariance_background(nvar,nvar);
//        matrixcovariance_background=CreateCovarianceMatrix(*matrixcorrelation_background,std_dev_background, nvar);
//        TMatrixD whiteningmatrix_signal(nvar,nvar);
//        TMatrixD whiteningmatrix_background(nvar,nvar);
//        whiteningmatrix_background=CreateCovarianceMatrix(*matrixcorrelation_background,std_dev_background, nvar);
//        std::cout<<"Correlation matrix (Background)"<<std::endl;
//        matrixcorrelation_background->Print();
//        std::cout<<"Covariance matrix (Background)"<<std::endl;
//        matrixcovariance_background.Print();
//        std::cout<<"Whitening matrix (Background)"<<std::endl;
//        whiteningmatrix_background.Print();

        TMatrixD mutualinformation_signal(nvar,nvar);
        GetMatrixMutualInformation(mutualinformation_signal,anaData,dataloader,nvar,"Signal", info1, info2);
        GetHistoMutualInformation(mutualinformation_signal,anaData,dataloader,nvar,"Signal");
        root_ext::WriteObject(anaData.mutualinformation("Signal"),outfile.get());
        std::vector<std::pair<double, std::string>> elements_mutualinformation_signal;
        elements_mutualinformation_signal=GetElements(mutualinformation_signal,dataloader,nvar);
        std::sort(elements_mutualinformation_signal.begin(),elements_mutualinformation_signal.end(), PairCompare);
        ReadVector(elements_mutualinformation_signal,args,"Signal","mutualinformation");

        TMatrixD mutualinformation_background(nvar,nvar);
        GetMatrixMutualInformation(mutualinformation_background,anaData,dataloader,nvar,"Background", info1, info2);
        GetHistoMutualInformation(mutualinformation_background,anaData,dataloader,nvar,"Background");
        root_ext::WriteObject(anaData.mutualinformation("Background"),outfile.get());
        std::vector<std::pair<double, std::string>> elements_mutualinformation_background;
        elements_mutualinformation_background=GetElements(mutualinformation_background,dataloader,nvar);
        std::sort(elements_mutualinformation_background.begin(),elements_mutualinformation_background.end(),PairCompare);
        ReadVector(elements_mutualinformation_background,args,"Background","mutualinformation");

    }

    //Create plots of Jensen Shannon between signal and background for different values of mass
    void KLDPlotSignalBkg(std::map<int, std::map<std::string,std::vector<double>>> sample_signal, std::map<int, std::map<std::string,std::vector<double>>> sample_bkg, std::map<int, std::map<std::string, double>> band_signal, std::map<int, std::map<std::string, double>> band_bkg, std::shared_ptr<TFile> outfile){
        std::map<int,std::vector<std::pair<std::string,double>>> JSDivergence;
        std::map<std::string, TGraph*> plot;
        auto map_bkg = sample_bkg[1];
        auto bandwidth_bkg = band_bkg[1];
        int i = 0;
        for(const auto& var_mass_2 : sample_signal){
            if (var_mass_2.second.size() == 0) continue;
            auto map_signal_2 = var_mass_2.second;
            auto bandwidth_signal = band_signal[var_mass_2.first];
            auto JSDivergence_mass = JSDivergence[var_mass_2.first];
            for (const auto& var : var_mass_2.second){
                double k = stat_estimators::JensenShannonDivergence(map_bkg[var.first], map_signal_2[var.first], bandwidth_bkg[var.first],bandwidth_signal[var.first]);
                JSDivergence_mass.emplace_back(var.first,k);
                if (plot.count(var.first)==0) plot[var.first] = new TGraph();
                plot[var.first]->SetPoint(i,var_mass_2.first, k);
                plot[var.first]->SetLineColor(kGreen+1);
                plot[var.first]->SetLineWidth(1);
                plot[var.first]->SetMarkerColor(1);
                plot[var.first]->SetMarkerSize(1);
                plot[var.first]->SetMarkerStyle(3);
                plot[var.first]->SetTitle(("JSD_"+var.first+"_SignalBkg").c_str());
                plot[var.first]->SetName(("JSD_"+var.first+"_SignalBkg").c_str());
                plot[var.first]->GetHistogram()->GetXaxis()->SetTitle("mass");
                plot[var.first]->GetHistogram()->GetYaxis()->SetTitle("JS Divergence");
            }
            std::sort(JSDivergence_mass.begin(),JSDivergence_mass.end(), [](std::pair<std::string,double> el1, std::pair<std::string,double> el2){
                return el1.second < el2.second;
            });
            JSDivergence[var_mass_2.first] = JSDivergence_mass;
            i++;
        }

        std::ofstream ofs("JSDivergence_signalvsbkg.csv", std::ofstream::out);
        for(const auto& var_mass_2 : sample_signal){
            ofs<<"Massa: "<<","<<var_mass_2.first<<std::endl;
            auto JSDivergence_mass = JSDivergence[var_mass_2.first];
            for (Long64_t j = 0; j < map_bkg.size(); j++){
                ofs<<JSDivergence_mass[j].first<<",";
            }
            ofs<<std::endl;
            for (Long64_t j = 0; j < map_bkg.size(); j++){
                ofs<<JSDivergence_mass[j].second<<",";
            }
            ofs<<std::endl;
        }
        ofs.close();

        for(const auto& var: map_bkg){
            root_ext::WriteObject(*plot[var.first], outfile.get());
        }

    }



private:
    Arguments args;
    SampleEntryCollection samples;
    std::shared_ptr<TFile> outfile;
    std::shared_ptr<TMVA::Factory> factory;
    std::shared_ptr<TMVA::DataLoader> dataloader;
    MvaVariables vars;
    MvaData anaData;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<> testvstraining;

};


}


PROGRAM_MAIN(analysis::MvaClassification, Arguments) // definition of the main program function

//./run.sh MvaVariableCorrelation --input_path ~/Desktop/tuples --output_file corr_lm_eTau.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name eTau --number_events 1000000
/*! Study of correlation matrix and mutual information of BDT variables
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/ClassInfo.h"
#include "TMVA/VariableTransformBase.h"
#include <TMatrixDEigen.h>
#include <TH1.h>
#include <TH2.h>
#include <fstream>
#include <random>

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(Long64_t, number_events);
};

namespace analysis{

class MvaData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    TH2D_ENTRY(histo, 100, 0, 0, 100, 0, 0)
    TH2D_ENTRY(mutualinformation, 100, 0, 0, 100, 0, 0)
};

struct SampleEntry{
  std::string filename;
  double weight;
  bool issignal;
  SampleEntry() : weight(-1), issignal(false){}
};

struct Variable{
  std::string varname;
  double weight;
  bool issignal;
  SampleEntry() : weight(-1), issignal(false){}
};

std::ostream& operator<<(std::ostream& os, const SampleEntry& entry)
{
    os << entry.filename << " " << entry.weight << " " << std::boolalpha << entry.issignal;
    return os;
}

std::istream& operator>>(std::istream& is, SampleEntry& entry)
{
    is >> entry.filename >> entry.weight >> std::boolalpha >> entry.issignal;
    return is;
}

class MvaVariables{
private:
    std::map<std::string, size_t> names;
public:
    std::vector<double> variables;
    MvaVariables(std::vector<SampleEntry> vector){}

    double& operator[](const std::string& name)
    {
      if (!names.count(name)){
          variables.push_back(0);
          names[name]=variables.size()-1;
          loader->AddVariable(name);
      }
      return variables.at(names.at(name));
    }

    void AddEvent(bool issignal, bool istraining, double weight)
    {
        const std::string samplename = issignal ? "Signal" : "Background";
        const TMVA::Types::ETreeType treetype= istraining ? TMVA::Types::kTraining : TMVA::Types::kTesting;
        loader->AddEvent(samplename, treetype, variables, weight);
    }

};


}

auto it = JSDivergence_vector.begin();
int count = 0;
while  (it != JSDivergence_vector.end()){
    std::cout<<"count "<<count<<std::endl;
    auto entry = JSDivergence_vector[count];
    std::cout<<"entry size "<<entry.first.names.size()<<std::endl;
    if (!entry.first.names.size()) break;
    bool check = true;
    if(entry.first.names.count(name)){
        if (*entry.first.names.begin() == name && *entry.first.names.rbegin() == name)  {
            ++it;
            std::cout<<"continue"<<std::endl;
            count++;
            continue;
        }
        std::string other;
        if (*entry.first.names.begin() == name )   other = *entry.first.names.rbegin();
        if (*entry.first.names.rbegin() == name )  other = *entry.first.names.begin();
        std::cout<<"other: "<<other<<std::endl;
        if (not_corrected.count(other)) {
            count++;
            continue;
        }
        const VarPair var_pair = name < other ?
                    VarPair(name, other) : VarPair(other, name);
        const double MI_sgn = mutual_matrix_signal.at(var_pair);
        const double MI_bkg = mutual_matrix_bkg.at(var_pair);
        if(MI_sgn < threashold || MI_bkg < threashold) {
            not_corrected.insert(other);
            std::cout<<"not corrected size "<<not_corrected.size()<<" "<<other<<std::endl;
            count--;
            JSDivergence_vector.erase(it);
            check = false;
        }
        else ++it;
    }
    else ++it;
    count++;
}
