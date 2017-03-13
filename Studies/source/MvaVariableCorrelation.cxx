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
    std::map<std::string, std::map<std::string, std::vector<double>>> all_variables;
    std::map<std::string, double> variables;

public:
    MvaVariables(){}

    double& operator[](const std::string& name)
    {
        return variables[name];
    }

    void AddEvent(const std::string& sample)
    {
        std::map<std::string, std::vector<double>>& sample_vars = all_variables[sample];
        for(const auto& name_value : variables) {
            const std::string& name = name_value.first;
            const double value = name_value.second;
            sample_vars[name].push_back(value);
        }
    }

    const std::map<std::string, std::vector<double>>& GetSampleVariables(const std::string& sample) const
    {
        return all_variables.at(sample);
    }

};

class MvaData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    TH1D_ENTRY(histo, 100, 0, 0)
    TH2D_ENTRY(matrix, 100, 0, 0, 100, 0, 0)
};

using VarPair = std::pair<std::string, std::string>;

//Create an histogram for each variable. Check if there are at least 10 events for each bin and if not rebin
// the histos.Compute chi2test for signal-background histos for each variable.
std::map<std::string, double> Chi2Test(std::map<std::string, std::vector<double>> sample_signal, std::map<std::string, std::vector<double>> sample_bkg, std::shared_ptr<TFile> outfile){

    std::map<std::string, double> pi_value;
    for(const auto& var_1 : sample_signal) {
        std::vector<double> vector_s = var_1.second;
        std::vector<double> vector_b = sample_bkg[var_1.first];
        int nbin = 50;
        TH1D* h_s = new TH1D((var_1.first+"Signal").c_str(),(var_1.first).c_str(),nbin,0,0);
        TH1D* h_b = new TH1D((var_1.first+"Bkg").c_str(),(var_1.first).c_str(),nbin,0,0);
        h_s->SetCanExtend(TH1::kAllAxes);
        h_b->SetCanExtend(TH1::kAllAxes);
        for(Long64_t i = 0; i < vector_s.size(); i++){
            h_s->Fill(vector_s[i]);
        }
        for(Long64_t i = 0; i < vector_b.size(); i++){
            h_b->Fill(vector_b[i]);
        }

        root_ext::WriteObject(*h_s, outfile.get());
        root_ext::WriteObject(*h_b, outfile.get());

        int new_nbin_s = 0;
        double rem_bin_s = 0;
        double bin_s[51] = {};
        bin_s[0]= h_s->GetBinLowEdge(1);
        for (Long64_t i=1; i<nbin; i++)
        {
            double entry = h_s->GetBinContent(i);
            while (entry < 10. && i<=50)
            {
                i++;
                entry = entry + h_s->GetBinContent(i);
                if (entry > 10) {
                    rem_bin_s = i;
                }
            }
            if (entry > 10){
                new_nbin_s++;
                bin_s[new_nbin_s]=h_s->GetBinLowEdge(i+1);
            }
        }

        int new_nbin_b = 0;
        double rem_bin_b = 0;
        double bin_b[51] = {};
        bin_b[0]= h_b->GetBinLowEdge(1);
        for (Long64_t i=1; i<nbin; i++)
        {
            double entry = h_b->GetBinContent(i);
            while (entry < 10. && i<50)
            {
                i++;
                entry = entry + h_b->GetBinContent(i);
                if (entry > 10) {
                    rem_bin_b = i;
                }
            }
            if (entry > 10){
                new_nbin_b++;
                bin_b[new_nbin_b]=h_b->GetBinLowEdge(i+1);
            }
        }

        std::cout<<"ciao"<<std::endl;
        TH1D* h_s_n = new TH1D((var_1.first+"_Signal").c_str(),(var_1.first).c_str(), new_nbin_s, bin_s);
        std::cout<<"!"<<std::endl;
        TH1D* h_b_n = new TH1D((var_1.first+"_Bkg").c_str(),(var_1.first).c_str(),new_nbin_b, bin_b);
        h_s_n->SetCanExtend(TH1::kAllAxes);
        h_b_n->SetCanExtend(TH1::kAllAxes);
        for(Long64_t i = 0; i < vector_s.size(); i++){
            h_s_n->Fill(vector_s[i]);
        }
        for(Long64_t i = 0; i < vector_b.size(); i++){
            h_b_n->Fill(vector_b[i]);
        }
        root_ext::WriteObject(*h_s_n, outfile.get());
        root_ext::WriteObject(*h_b_n, outfile.get());
//        h_s_n->Scale(1/h_s_n->Integral());
//        h_b_n->Scale(1/h_b_n->Integral());

//        double a = h_s_n->Chi2Test(h_b_n, "WW");
//        std::cout<<a<<std::endl;
//        pi_value[var_1.first] = a;
        h_s->Delete();
        h_b->Delete();
        h_s_n->Delete();
        h_b_n->Delete();
     }
     return pi_value;
}


//Create covariance of two variables
double Covariance (const std::vector<double> vec_1, const std::vector<double> vec_2){
    double mean_1 = 0, mean_2 = 0;
    for (Long64_t i=0; i<vec_1.size(); i++ ){
        mean_1 = mean_1 + vec_1[i];
        mean_2 = mean_2 + vec_2[i];
    }
    mean_1 = mean_1/vec_1.size();
    mean_2 = mean_2/vec_1.size();
    double cov = 0;
    for (Long64_t i=0; i<vec_1.size(); i++ ){
        cov = cov + (vec_1[i] - mean_1) * (vec_2[i] - mean_2);
    }
    return cov/(vec_1.size()-1);
}

//Create correlation of two variables from covariance
std::map<VarPair, double> CovToCorr(std::map<VarPair, double> cov_matrix){
    std::map<VarPair, double> corr;
    for(const auto& elements : cov_matrix) {
        auto name_1 = elements.first;
        const VarPair var_11(name_1.first, name_1.first);
        double sigma_1 = std::sqrt(cov_matrix[var_11]);
        const VarPair var_22(name_1.second, name_1.second);
        double sigma_2 = std::sqrt(cov_matrix[var_22]);
        const VarPair var_12(name_1.first, name_1.second);
        corr[var_12] = cov_matrix[var_12]/(sigma_1*sigma_2);
    }
    return corr;
}


//Create covariance(correlation) histo matrix
void CreateMatrixHistos(std::map<std::string, std::vector<double>> sample_vars, std::map<VarPair, double> element, MvaData &anaData, std::string type,  std::string class_sample){
    int bin =  sample_vars.size();
    anaData.matrix(type+"_"+class_sample).SetCanExtend(TH1::kAllAxes);
    anaData.matrix(type+"_"+class_sample).SetBins(bin, 0, bin, bin, 0, bin);
    int i = 1;
    for(const auto& var_1 : sample_vars) {
        int j = 1;
        anaData.matrix(type+"_"+class_sample).GetXaxis()->SetBinLabel(i, (var_1.first).c_str());
        anaData.matrix(type+"_"+class_sample).GetYaxis()->SetBinLabel(i, (var_1.first).c_str());
        for(const auto& var_2 : sample_vars) {
            if (j>=i){
            const VarPair var_12(var_1.first, var_2.first);
            anaData.matrix(type+"_"+class_sample).SetBinContent(i, j, element[var_12]);
            anaData.matrix(type+"_"+class_sample).SetBinContent(j, i, element[var_12]);
            }
            j++;
        }
        i++;
    }
}


//Remove diagonal and under diagonal elements from a matrix
std::map<VarPair, double> RemoveDiagonal(std::map<VarPair, double> matrix){
    std::map<VarPair, double> matrix_diagonal;
    matrix_diagonal = matrix;
    for(const auto& elements : matrix) {
        auto name = elements.first;
        const VarPair var_11(name.first, name.first);
        matrix_diagonal.erase(var_11);
    }
    return matrix_diagonal;
}

//Difference between two vectors of pairs
std::vector< std::pair<VarPair, double>> Difference(std::vector< std::pair<VarPair, double>> vec_1, std::vector< std::pair<VarPair, double>> vec_2){
    std::vector< std::pair<VarPair, double>> diff;
    for(Long64_t i = 0; i < vec_1.size(); i++){
        std::pair<VarPair, double> element_1 = vec_1[i];
        std::pair<VarPair, double> element_2 = vec_2[i];
        if (element_1.first == element_2.first){
            VarPair nome=element_1.first;
            VarPair nome2=element_2.first;
            double difference = std::abs(element_1.second - element_2.second);
            diff.push_back(std::make_pair(element_1.first,difference));
        }
    }
    return diff;
}

//Compute Kolmogorov test
std::map< std::string, double> KolmogorovTest(std::map<std::string, std::vector<double>> sample_signal, std::map<std::string, std::vector<double>> sample_bkg){
    std::map< std::string, double> kolmogorov;
    double k;
    for(const auto& var_1 : sample_signal) {
        std::vector<double> vector_signal = sample_signal[var_1.first];
        std::sort(vector_signal.begin(), vector_signal.end());
        std::vector<double> vector_bkg = sample_bkg[var_1.first];
        std::sort(vector_bkg.begin(), vector_bkg.end());
        Double_t v_s[vector_signal.size()], v_b[vector_bkg.size()];
        for(Long64_t i = 0; i < vector_signal.size(); i++){
            v_s[i]=vector_signal[i];
            v_b[i]=vector_bkg[i];
        }
        k = TMath::KolmogorovTest(vector_signal.size(), v_s, vector_bkg.size(), v_b, "M");
        kolmogorov[var_1.first] = k;
    }
    return kolmogorov;
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
        outfile(root_ext::CreateRootFile(args.output_file())), vars(),  anaData(outfile)
    {
    }

    void Run()
    {
        const double mass_top = 173.21;
        for(const SampleEntry& entry:samples)
        {
            auto input_file=root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
            EventTuple tuple(args.tree_name(), input_file.get(), true, GetDisabledBranches());
            std::cout<<entry<<" number of events: "<<std::min(tuple.GetEntries(),args.number_events())<<std::endl;
            std::string samplename = entry.issignal ? "Signal" : "Background";

            for(Long64_t current_entry = 0; current_entry < std::min(tuple.GetEntries(),args.number_events()); ++current_entry) {
                tuple.GetEntry(current_entry);
                const Event& event = tuple.data();

                if (event.eventEnergyScale!=0 || (event.q_1+event.q_2)!=0 || event.jets_p4.size() < 2
                    || event.extraelec_veto==true || event.extramuon_veto==true) continue;

                LorentzVectorE_Float bb= event.jets_p4[0] + event.jets_p4[1];
                LorentzVectorM_Float leptons= event.p4_1 + event.p4_2;
                LorentzVectorM_Float leptonsMET= event.p4_1 + event.p4_2 + event.pfMET_p4;

                double circular_cut=std::sqrt(pow(event.SVfit_p4.mass()-116.,2)+pow(bb.M()-111,2));
                if (circular_cut>40) continue;
//                if ((args.tree_name=="eTau" && circular_cut>40)||(args.tree_name=="muTau" && circular_cut>30)) continue;

//                vars["pt_l1"]=event.p4_1.pt();
                vars["pt_l2"]=event.p4_2.pt();
//                vars["pt_jet1"]=event.jets_p4[1].pt();
                vars["pt_jet2"]=event.jets_p4[2].pt();
//                vars["pt_leptons"]=leptons.pt();
//                vars["pt_sv"]=event.SVfit_p4.pt();
//                vars["pt_leptonsmet"]=leptonsMET.pt();
//                vars["pt_bb"]=bb.pt();
//                vars["pt_met"]=event.pfMET_p4.pt();

                vars["p_zeta"]=

//                vars["abs(dphi_leptons)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.p4_2));
                vars["dphi_leptons"] = (ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.p4_2));
//                vars["abs(dphi_jets)"]=std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.jets_p4[0], event.jets_p4[1]));
                vars["dphi_jets"]=(ROOT::Math::VectorUtil::DeltaPhi(event.jets_p4[0], event.jets_p4[1]));
//                vars["abs(dphi_l1met)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.pfMET_p4));
                vars["dphi_l1met"] = (ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.pfMET_p4));
//                vars["abs(dphi_l2met)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_2, event.pfMET_p4));
//                vars["dphi_l2met"] = (ROOT::Math::VectorUtil::DeltaPhi(event.p4_2, event.pfMET_p4));
//                vars["abs(dphi_svmet)"]=std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4));
//                vars["dphi_svmet"]=(ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4));
//                vars["abs(dphi_bbmet)"]=std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.pfMET_p4));
                vars["dphi_bbmet"]=(ROOT::Math::VectorUtil::DeltaPhi(bb, event.pfMET_p4));
//                vars["abs(dphi_bbsv)"]=std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.SVfit_p4));
                vars["dphi_bbsv"]=(ROOT::Math::VectorUtil::DeltaPhi(bb, event.SVfit_p4));

                vars["abs(deta_leptons)"] =std::abs(event.p4_1.eta() - event.p4_2.eta());
//                vars["deta_leptons"] =(event.p4_1.eta() - event.p4_2.eta());
                vars["abs(deta_jets)"]=std::abs(event.jets_p4[0].eta() - event.jets_p4[1].eta());
                vars["deta_jets"]=(event.jets_p4[0].eta() - event.jets_p4[1].eta());
//                vars["abs(deta_l1met)"] = std::abs(event.p4_1.eta()-event.pfMET_p4.eta());
//                vars["deta_l1met"] = (event.p4_1.eta()-event.pfMET_p4.eta());
//                vars["abs(deta_l2met)"] = std::abs(event.p4_2.eta()-event.pfMET_p4.eta());
//                vars["deta_l2met"] = (event.p4_2.eta()-event.pfMET_p4.eta());
//                vars["abs(deta_svmet)"]=std::abs(event.SVfit_p4.eta()-event.pfMET_p4.eta());
//                vars["deta_svmet"]=(event.SVfit_p4.eta()-event.pfMET_p4.eta());
//                vars["abs(deta_bbmet)"]=std::abs(bb.eta()-event.pfMET_p4.eta());
//                vars["deta_bbmet"]=(bb.eta()-event.pfMET_p4.eta());
                vars["abs(deta_bbsv)"]=std::abs(bb.eta()-event.SVfit_p4.eta());
//                vars["deta_bbsv"]=bb.eta()-event.SVfit_p4.eta();

//                vars["dR_leptons"] = std::abs(ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2));
//                vars["dR_jets"]=std::abs(ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]));
//                vars["dR_l1met"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.pfMET_p4));
//                vars["dR_l2met"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_2, event.pfMET_p4));
//                vars["dR_svmet"]=(ROOT::Math::VectorUtil::DeltaR(event.SVfit_p4, event.pfMET_p4));
//                vars["dR_bbmet"]=(ROOT::Math::VectorUtil::DeltaR(bb, event.pfMET_p4));
//                vars["dR_bbsv"]=(ROOT::Math::VectorUtil::DeltaR(bb, event.SVfit_p4));

//                vars["dR_jetsPt_bb"]=std::abs(ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]))*bb.Pt();
//                vars["dR_leptonsPt_sv"] = std::abs(ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2))*event.SVfit_p4.Pt();

//                vars["mass_leptonsmetinvariant"]=ROOT::Math::VectorUtil::InvariantMass(leptons,event.pfMET_p4); //
                vars["mass_sv"]=event.SVfit_p4.M();
                vars["mass_jets"]=bb.M();
//                vars["MT_l1"]=Calculate_MT(event.p4_1,event.pfMET_p4);
//                vars["MT_l2"]=Calculate_MT(event.p4_2,event.pfMET_p4);
//                vars["MT_tautau"]=std::sqrt(pow(event.p4_1.Et()+event.p4_1.Et(),2)+pow(event.p4_1.px()+event.p4_1.px(),2)+pow(event.p4_1.py()+event.p4_1.py(),2));//
//                vars["MT_tautausv"]=Calculate_MT(event.SVfit_p4, event.pfMET_p4); //
//                vars["MT_leptons"]=Calculate_MT(leptons, event.pfMET_p4);
//                vars["mass_H"]=ROOT::Math::VectorUtil::InvariantMass(bb,event.SVfit_p4);

                LorentzVectorM_Float t1 = event.p4_1 + event.jets_p4[1] + event.pfMET_p4;
                LorentzVectorM_Float t2 = event.p4_2 + event.jets_p4[1] + event.pfMET_p4;
                LorentzVectorM_Float t3 = event.p4_1 + event.jets_p4[2] + event.pfMET_p4;
                LorentzVectorM_Float t4 = event.p4_2 + event.jets_p4[4] + event.pfMET_p4;
                double d1= std::abs(t1.mass() - mass_top);
                double d2= std::abs(t2.mass() - mass_top);
                double d3= std::abs(t3.mass() - mass_top);
                double d4= std::abs(t4.mass() - mass_top);
                if (d1<d2 && d1<d3 && d1<d4) vars["mass_top"] = d1;
                if (d2<d1 && d2<d3 && d2<d4) vars["mass_top"] = d2;
                if (d3<d1 && d3<d2 && d3<d4) vars["mass_top"] = d3;
                if (d4<d1 && d4<d3 && d4<d2) vars["mass_top"] = d4;

//                vars["theta_1"]=ROOT::Math::VectorUtil::Angle(event.p4_1,event.SVfit_p4);
//                vars["theta_2"]=ROOT::Math::VectorUtil::Angle(event.jets_p4[1],bb);

                vars.AddEvent(samplename);
            }
        }

        std::map<std::string, std::vector<double>> sample_vars_signal = vars.GetSampleVariables("Signal");
        std::map<std::string, std::vector<double>> sample_vars_bkg = vars.GetSampleVariables("Background");

        std::map<VarPair, double> cov_matrix_signal, cov_matrix_bkg;
        for(const auto& var_1 : sample_vars_signal) {
            for(const auto& var_2 : sample_vars_signal) {
                const VarPair var_12(var_1.first, var_2.first);
                const VarPair var_21(var_2.first, var_1.first);
                if(cov_matrix_signal.count(var_21)) continue;
                const double cov = Covariance(var_1.second, var_2.second);
                cov_matrix_signal[var_12] = cov;
            }
        }

        for(const auto& var_1 : sample_vars_bkg) {
            for(const auto& var_2 : sample_vars_bkg) {
                const VarPair var_12(var_1.first, var_2.first);
                const VarPair var_21(var_2.first, var_1.first);
                if(cov_matrix_bkg.count(var_21)) continue;
                const double cov = Covariance(var_1.second, var_2.second);
                cov_matrix_bkg[var_12] = cov;
            }
        }

//        CreateMatrixHistos(sample_vars_signal,cov_matrix_signal,anaData,"covariance","Signal");
//        CreateMatrixHistos(sample_vars_bkg,cov_matrix_bkg,anaData,"covariance","Background");

        std::map<VarPair, double> corr_matrix_signal, corr_matrix_bkg;
        corr_matrix_signal = CovToCorr(cov_matrix_signal);
        corr_matrix_bkg = CovToCorr(cov_matrix_bkg);

        CreateMatrixHistos(sample_vars_signal,corr_matrix_signal,anaData,"correlation","Signal");
        CreateMatrixHistos(sample_vars_bkg,corr_matrix_bkg,anaData,"correlation","Background");

        std::map<VarPair, double> corr_matrix_signal_fix, corr_matrix_bkg_fix;
        corr_matrix_signal_fix = RemoveDiagonal(corr_matrix_signal);
        corr_matrix_bkg_fix = RemoveDiagonal(corr_matrix_bkg);

        std::vector< std::pair<VarPair, double>> corr_vector_signal(corr_matrix_signal_fix.begin(), corr_matrix_signal_fix.end());
        std::vector< std::pair<VarPair, double>> corr_vector_bkg(corr_matrix_bkg_fix.begin(), corr_matrix_bkg_fix.end());
        std::vector< std::pair<VarPair, double>> corr_vector_difference = Difference(corr_vector_signal,corr_vector_bkg);
        std::map< std::string, double> kolmogorov = KolmogorovTest(sample_vars_signal, sample_vars_bkg);
        //std::map< std::string, double> p_value = Chi2Test(sample_vars_signal, sample_vars_bkg, outfile);

        std::sort(corr_vector_difference.begin(), corr_vector_difference.end(), [](const std::pair<VarPair, double>& a, const std::pair<VarPair, double>& b) {
            return std::abs(a.second) < std::abs(b.second);
        });

        std::map<std::string, int> eliminate;
        std::ofstream ofs(args.tree_name()+".txt", std::ofstream::out);
        ofs << "Variable    " << "  Variable2   " << "  corr_s  " << "  corr_b    " << "  corr_d    "<<"  pvalue_1    "<<"  pvalue_2    "<<"  KS_1    "<<"  KS_2"<<std::endl;
        for(Long64_t i = 0; i < corr_vector_difference.size(); i++){
            std::pair<VarPair, double> element_difference = corr_vector_difference[i];
            VarPair pair = element_difference.first;

            if (/*(std::abs(corr_matrix_signal_fix[pair])>0.2 && std::abs(corr_matrix_bkg_fix[pair])>0.2) &&*/ element_difference.second<0.2){
                ofs <<pair.first<<"    "<<pair.second<<"    "<<std::abs(corr_matrix_signal_fix[pair])<<"    "<<std::abs(corr_matrix_bkg_fix[pair])<<"    "<<element_difference.second<</*"   "<<p_value[pair.first]<<"    "<<p_value[pair.second]<<*/"  "<<kolmogorov[pair.first]<<"    "<<kolmogorov[pair.second]<<std::endl;
                if (kolmogorov[pair.first]<=kolmogorov[pair.second] && eliminate.count(pair.first)!=0) eliminate[pair.first]++;
                if (kolmogorov[pair.first]<=kolmogorov[pair.second] && eliminate.count(pair.first)==0) eliminate[pair.first]=1;
                if (kolmogorov[pair.first]>kolmogorov[pair.second] && eliminate.count(pair.second)!=0) eliminate[pair.second]++;
                if (kolmogorov[pair.first]>kolmogorov[pair.second] && eliminate.count(pair.second)==0) eliminate[pair.second]=1;

            }
        }
        ofs.close();
        std::ofstream off("elimination.txt", std::ofstream::out);
        for(const auto& i : eliminate){
            off << i.first << "     " << i.second<<std::endl;
        }
        off.close();
    }
private:
    Arguments args;
    SampleEntryCollection samples;
    std::shared_ptr<TFile> outfile;
    MvaVariables vars;
    MvaData anaData;

};


}


PROGRAM_MAIN(analysis::MvaClassification, Arguments) // definition of the main program function

//dies/config/mva_config.cfg --tree_name eTau --number_events 1000000








