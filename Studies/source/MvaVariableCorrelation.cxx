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
#include <algorithm>

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
    std::map<std::string, std::map<int, std::map<std::string, std::vector<double>>>> all_variables;
    std::map<std::string, double> variables;

public:
    MvaVariables(){}

    double& operator[](std::string name)
    {
        return variables[name];
    }

    void AddEvent(const std::string& sample, int mass)
           {
               std::map<int, std::map<std::string, std::vector<double>>>& sample_map = all_variables[sample];
               std::map<std::string, std::vector<double>>& sample_vars = sample_map[mass];
               for(const auto& name_value : variables) {
                   const std::string& name = name_value.first;
                   const double value = name_value.second;
                   sample_vars[name].push_back(value);
               }
           }

    const std::map<int, std::map<std::string,std::vector<double>>>& GetSampleVariables(const std::string& sample) const
    {
        return all_variables.at(sample);
    }
};

using VarPair = std::pair<std::string,std::string>;

/*Create a pair of histogram for each variable. The first one is for signal, the second for background.
They are binned in the same way, to have at least 10 entries for bin.*/
std::map<int, std::map<std::string, std::pair<TH1D*,TH1D*>>> GetHistos(std::map<int, std::map<std::string,std::vector<double>>> sample_signal, std::map<int, std::map<std::string,std::vector<double>>> sample_bkg, std::shared_ptr<TFile> outfile){

    std::map<int, std::map<std::string, std::pair<TH1D*,TH1D*>>> histogram;

    for(const auto& var_m : sample_signal){
        std::map<std::string, std::pair<TH1D*,TH1D*>> histos;
        std::map<std::string,std::vector<double>> v_b = sample_bkg[1];
        std::string mass;
        if (var_m.first!=1) mass = "_mass"+std::to_string(var_m.first);
        else mass = "";
        std::cout<<"MASSA: "<<var_m.first<<std::endl;
        for(const auto& var_1 : var_m.second) {
            std::vector<double> vector_s = var_1.second;
            std::vector<double> vector_b = v_b[var_1.first];
            std::cout<<"variabile: "<<var_1.first<<std::endl;
            std::cout<<"eventi segnale: "<<vector_s.size()<<std::endl;
            std::cout<<"eventi fondo: "<<vector_b.size()<<std::endl;

            auto max_s = std::max_element(vector_s.begin(),vector_s.end());
            auto min_s = std::min_element(vector_s.begin(),vector_s.end());
            auto max_b = std::max_element(vector_b.begin(),vector_b.end());
            auto min_b = std::min_element(vector_b.begin(),vector_b.end());
            double max, min;
            if (*min_s<=*min_b) min = *min_s;
            else min = *min_b;
            if (*max_s>=*max_b) max = *max_s;
            else max = *max_b;

            const int nbin = 50;
            TH1D* h_s = new TH1D((var_1.first+"Signal"+mass).c_str(),(var_1.first+"Signal"+mass).c_str(),nbin,min,max);
            h_s->SetCanExtend(TH1::kAllAxes);
            h_s->SetXTitle((var_1.first).c_str());
            for(Long64_t i = 0; i < vector_s.size(); i++){
                h_s->Fill(vector_s[i]);
            }
//            root_ext::WriteObject(*h_s, outfile.get());
            int new_nbin_s = 0;
            double bin_s[nbin] = {};
            for (Long64_t i=0; i<=nbin; i++)
            {
                double entry = h_s->GetBinContent(i);
                Long64_t rebin = i;
                while (entry <= 10 && i<=nbin)
                {
                    i++;
                    entry = entry + h_s->GetBinContent(i);
                }
                if (entry > 10){
                    new_nbin_s++;
                    bin_s[new_nbin_s]=h_s->GetBinLowEdge(rebin);
                }
            }
            new_nbin_s++;
            bin_s[new_nbin_s]=h_s->GetBinLowEdge(nbin+1);
            bin_s[0]=bin_s[1]-1;
            h_s->Delete();

            histos[var_1.first].first = new TH1D((var_1.first+"Signal1"+mass).c_str(),(var_1.first+"Signal1"+mass).c_str(), new_nbin_s, bin_s);
            histos[var_1.first].first->SetXTitle((var_1.first).c_str());
            for(Long64_t i = 0; i < vector_s.size(); i++){
                histos[var_1.first].first->Fill(vector_s[i]);
            }
            histos[var_1.first].second = new TH1D((var_1.first+"Bkg1"+mass).c_str(),(var_1.first+"Bkg1"+mass).c_str(),new_nbin_s, bin_s);
            histos[var_1.first].second->SetXTitle((var_1.first).c_str());
            for(Long64_t i = 0; i < vector_b.size(); i++){
                histos[var_1.first].second->Fill(vector_b[i]);
            }
//            root_ext::WriteObject(*histos[var_1.first].first, outfile.get());
//            root_ext::WriteObject(*histos[var_1.first].second, outfile.get());

            int new_nbin_b = 0;
            double bin_b[nbin] = {};
            for (Long64_t i=0; i<=new_nbin_s; i++)
            {
                double entry = histos[var_1.first].second->GetBinContent(i);
                Long64_t rebin = i;
                while (entry <= 10. && i<=new_nbin_s)
                {
                    i++;
                    entry = entry + histos[var_1.first].second->GetBinContent(i);
                }
                if (entry > 10){
                    new_nbin_b++;
                    bin_b[new_nbin_b]=histos[var_1.first].second->GetBinLowEdge(rebin);
                }
            }
            new_nbin_b++;
            bin_b[new_nbin_b]=histos[var_1.first].second->GetBinLowEdge(new_nbin_s+1);
            bin_b[0]=bin_b[1]-1;

            bool check = false;
            int count = 0;
            while(check==false && count<nbin){
                if (bin_s[count]!=bin_b[count]) check=true;
                count++;
            }
            if (check==true){
                histos[var_1.first].first->Delete();
                histos[var_1.first].second->Delete();
                histos[var_1.first].first = new TH1D((var_1.first+"Signal2"+mass).c_str(),(var_1.first+"Signal2"+mass).c_str(), new_nbin_b, bin_b);
                histos[var_1.first].first->SetXTitle((var_1.first).c_str());
                for(Long64_t i = 0; i < vector_s.size(); i++){
                    histos[var_1.first].first->Fill(vector_s[i]);
                }
                histos[var_1.first].second = new TH1D((var_1.first+"Bkg2"+mass).c_str(),(var_1.first+"Bkg2"+mass).c_str(),new_nbin_b, bin_b);
                histos[var_1.first].second->SetXTitle((var_1.first).c_str());
                for(Long64_t i = 0; i < vector_b.size(); i++){
                    histos[var_1.first].second->Fill(vector_b[i]);
                }
//                root_ext::WriteObject(*histos[var_1.first].first, outfile.get());
//                root_ext::WriteObject(*histos[var_1.first].second, outfile.get());
                histos[var_1.first].first->Scale(1/histos[var_1.first].first->Integral());
                histos[var_1.first].second->Scale(1/histos[var_1.first].second->Integral());
            }
            else{
//                root_ext::WriteObject(*histos[var_1.first].first, outfile.get());
//                root_ext::WriteObject(*histos[var_1.first].second, outfile.get());
                histos[var_1.first].first->Scale(1/histos[var_1.first].first->Integral());
                histos[var_1.first].second->Scale(1/histos[var_1.first].second->Integral());
            }
         }

        histogram[var_m.first] = histos;
    }
    return histogram;
}

//Chi2Test for signal vs background histo for each variable
std::map< int, std::map<std::string,double>> Chi2Test(std::map<int, std::map<std::string, std::pair<TH1D*,TH1D*>>> h){

    std::map< int, std::map<std::string,double>> chi;
    for (const auto& var_m : h){
        std::cout<<var_m.first<<std::endl;
        std::map<std::string,double> chi2;
        for(const auto& var_1 : var_m.second) {
            std::cout<<var_1.first<<std::endl;
            std::map<std::string, std::pair<TH1D*,TH1D*>> histos = h[var_m.first];
            std::pair<TH1D*,TH1D*> pair_histo = histos[var_1.first];
            double a = pair_histo.first->Chi2Test(pair_histo.second, "WWCHI2/NDF");
            chi2[var_1.first] = a;
        }
        chi[var_m.first] = chi2;
    }
    return chi;
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

std::map<int, std::map<VarPair,double>> Cov(std::map<int, std::map<std::string,std::vector<double>>> sample_vars){
    std::map<int, std::map<VarPair,double>> cov_matrix;
    for (const auto& var: sample_vars){
        if (var.second.size() == 0) continue;
        const auto map = var.second;
        std::map<VarPair,double> cov;
        for(const auto& var_1 : map) {
            for(const auto& var_2 : map) {
                const VarPair var_12(var_1.first, var_2.first);
                const VarPair var_21(var_2.first, var_1.first);
                if(cov.count(var_21)) continue;
                cov[var_12] = Covariance(var_1.second, var_2.second);
            }
        }
        cov_matrix[var.first] = cov;
    }
    return cov_matrix;
}

//Create correlation of two variables from covariance
std::map<int, std::map<VarPair,double>> CovToCorr(std::map<int, std::map<VarPair,double>> cov_matrix){
    std::map<int, std::map<VarPair,double>> corr;

    for(const auto& var : cov_matrix){
        if (var.second.size() == 0) continue;
        std::map<VarPair,double> correlation;
        auto covariance = var.second;
        for(const auto& elements : var.second) {
            auto name_1 = elements.first;
            const VarPair var_11(name_1.first, name_1.first);
            double sigma_1 = std::sqrt(covariance[var_11]);
            const VarPair var_22(name_1.second, name_1.second);
            double sigma_2 = std::sqrt(covariance[var_22]);
            const VarPair var_12(name_1.first, name_1.second);
            correlation[var_12] = covariance[var_12]/(sigma_1*sigma_2);
        }
        corr[var.first] = correlation;
    }
    return corr;
}

//Create covariance(correlation) histo matrix
void CreateMatrixHistos(std::map<int, std::map<std::string,std::vector<double>>> sample_vars, std::map<int, std::map<VarPair,double>> element, std::shared_ptr<TFile> outfile, std::string type,  std::string class_sample){
    for(const auto var_m: sample_vars){
        if (var_m.second.size() == 0) continue;
        int bin =  var_m.second.size();
        auto el = element[var_m.first];
        std::string mass;
        if (var_m.first!=1) mass = "_mass"+std::to_string(var_m.first);
        else mass = "";
        TH2D* matrix= new TH2D((type+"_"+class_sample+mass).c_str(),(type+"_"+class_sample+mass).c_str(),bin,0,0,bin,0,0);
        matrix->SetCanExtend(TH1::kAllAxes);
        matrix->SetBins(bin, 0, bin, bin, 0, bin);
        int i = 1;
        for(const auto& var_1 : var_m.second) {
            int j = 1;
            matrix->GetXaxis()->SetBinLabel(i, (var_1.first).c_str());
            matrix->GetYaxis()->SetBinLabel(i, (var_1.first).c_str());
            for(const auto& var_2 : var_m.second) {
                if (j>=i){
                const VarPair var_12(var_1.first, var_2.first);
                matrix->SetBinContent(i, j, el[var_12]);
                matrix->SetBinContent(j, i, el[var_12]);
                }
                j++;
            }
            i++;
        }
        root_ext::WriteObject(*matrix, outfile.get());
    }
}

//Remove diagonal elements from a symmetric matrix
std::map<int, std::map<VarPair,double>> RemoveDiagonal(std::map<int, std::map<VarPair,double>> matrix){
    std::map<int, std::map<VarPair,double>> matrix_diagonal;
    matrix_diagonal = matrix;
    for(const auto& var_m: matrix){
        if (var_m.second.size() == 0) continue;
        std::map<VarPair,double> mat = matrix_diagonal[var_m.first];
        for(const auto& elements : var_m.second) {
            const auto name = elements.first;
            const VarPair var_11(name.first, name.first);
            mat.erase(var_11);
        }
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

//Kolmogorov test return the maximum distance
std::map<int, std::map<std::string,double>> KolmogorovTest(std::map<int, std::map<std::string,std::vector<double>>> sample_signal, std::map<int, std::map<std::string,std::vector<double>>> sample_bkg){
    std::map<int, std::map<std::string,double>> kolmogorov;
    for(const auto& var_mass : sample_signal) {
        if (var_mass.second.size() == 0) continue;
        auto map_signal = var_mass.second;
        auto map_bkg = sample_bkg[1];
        std::map<std::string,double> k;
        for (const auto& var : var_mass.second){
            std::vector<double> vector_signal = map_signal[var.first];
            std::sort(vector_signal.begin(), vector_signal.end());
            std::vector<double> vector_bkg = map_bkg[var.first];
            std::sort(vector_bkg.begin(), vector_bkg.end());
            Double_t v_s[vector_signal.size()], v_b[vector_bkg.size()];
            for(Long64_t i = 0; i < vector_signal.size(); i++){
                v_s[i] = vector_signal[i];
            }
            for(Long64_t i = 0; i < vector_bkg.size(); i++){
                v_b[i] = vector_bkg[i];
            }
            k[var.first] = TMath::KolmogorovTest(vector_signal.size(), v_s, vector_bkg.size(), v_b, "");
        }
        kolmogorov[var_mass.first] = k;
    }
    return kolmogorov;
}

//Create plots of Kolmogorv probability in function of mass for each variable
void KolmogorvPlot(std::map<int, std::map<std::string,double>> kolmogorov, std::shared_ptr<TFile> outfile){
    std::map<std::string, TGraph*> plot;
    std::map<std::string,double> map_mass;
    int i = 0;
    for(const auto& var_mass : kolmogorov) {
        map_mass = var_mass.second;
        if (var_mass.second.size() == 0) continue;
        auto map = var_mass.second;
        for (const auto& var : map){
            if (plot.count(var.first)==0) plot[var.first] = new TGraph(kolmogorov.size());
            plot[var.first]->SetPoint(i,var_mass.first,var.second);
            plot[var.first]->SetLineColor(1);
            plot[var.first]->SetLineWidth(1);
            plot[var.first]->SetMarkerColor(1);
            plot[var.first]->SetMarkerSize(0.5);
            plot[var.first]->SetMarkerStyle(3);
            plot[var.first]->SetTitle((var.first).c_str());
            plot[var.first]->SetName((var.first).c_str());
        }
        i++;
    }
    for(const auto& var: map_mass){
        root_ext::WriteObject(*plot[var.first], outfile.get());
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
        outfile(root_ext::CreateRootFile(args.output_file())), vars()
    {
    }

    void Run()
    {
        const double mass_top = 173.21;
//        int m[19] = {250,260,270,280,300,320,340,350,400,450,500,550,600,650,700,750,800,900,1};
//        int i=0;
        for(const SampleEntry& entry:samples)
        {
//            int mass = m[i];
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

                vars["pt_l1"] = event.p4_1.pt();
                vars["pt_l2"] = event.p4_2.pt();
                vars["pt_b1"] = event.jets_p4[0].pt();
                vars["pt_b2"] = event.jets_p4[1].pt();
                vars["pt_l1+l2"] = leptons.pt();
                vars["pt_htautau"] = event.SVfit_p4.pt();
                vars["pt_l1+l2+met"] = leptonsMET.pt();
                vars["pt_hbb"] = bb.pt();
                vars["pt_met"] = event.pfMET_p4.pt();

                vars["p_zeta"] = Calculate_Pzeta(event.p4_1, event.p4_2, event.pfMET_p4);
                vars["p_zeta_visible"] = Calculate_visiblePzeta(event.p4_1,event.p4_2);

                vars["abs(dphi_l1-l2)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.p4_2));
                vars["dphi_l1-l2"] = (ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.p4_2));
                vars["abs(dphi_b1-b2)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.jets_p4[0], event.jets_p4[1]));
                vars["dphi_b1-b2"] = (ROOT::Math::VectorUtil::DeltaPhi(event.jets_p4[0], event.jets_p4[1]));
                vars["abs(dphi_l1-met)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.pfMET_p4));
                vars["dphi_l1-met"] = (ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.pfMET_p4));
                vars["abs(dphi_l2-met)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_2, event.pfMET_p4));
                vars["dphi_l2-met"] = (ROOT::Math::VectorUtil::DeltaPhi(event.p4_2, event.pfMET_p4));
                vars["abs(dphi_l1+l2-met)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(leptons, event.pfMET_p4));
                vars["dphi_l1+l2-met"] = (ROOT::Math::VectorUtil::DeltaPhi(leptons, event.pfMET_p4));
                vars["abs(dphi_htautau-met)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4));
                vars["dphi_htautau-met"] = (ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4));
                vars["abs(dphi_hbb-met)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.pfMET_p4));
                vars["dphi_hbb-met"] = (ROOT::Math::VectorUtil::DeltaPhi(bb, event.pfMET_p4));
                vars["abs(dphi_hbb-hatutau)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.SVfit_p4));
                vars["dphi_hbb-htautau"] = (ROOT::Math::VectorUtil::DeltaPhi(bb, event.SVfit_p4));

                vars["abs(deta_l1+l2)"] = std::abs(event.p4_1.eta() - event.p4_2.eta());
                vars["deta_l1+l2"] = (event.p4_1.eta() - event.p4_2.eta());
                vars["abs(deta_b1+b2)"] = std::abs(event.jets_p4[0].eta() - event.jets_p4[1].eta());
                vars["deta_b1+b2"] = (event.jets_p4[0].eta() - event.jets_p4[1].eta());
                vars["abs(deta_l1+met)"] = std::abs(event.p4_1.eta()-event.pfMET_p4.eta());
                vars["deta_l1+met"] = (event.p4_1.eta()-event.pfMET_p4.eta());
                vars["abs(deta_l2+met)"] = std::abs(event.p4_2.eta()-event.pfMET_p4.eta());
                vars["deta_l2+met"] = (event.p4_2.eta()-event.pfMET_p4.eta());
                vars["abs(deta_l1+l2-met)"] = std::abs(leptons.eta()-event.pfMET_p4.eta());
                vars["deta_l1+l2-met"] = (leptons.eta()-event.pfMET_p4.eta());
                vars["abs(deta_htautau-met)"] = std::abs(event.SVfit_p4.eta()-event.pfMET_p4.eta());
                vars["deta_hatutau-met"] = (event.SVfit_p4.eta()-event.pfMET_p4.eta());
                vars["abs(deta_hbb-met)"] = std::abs(bb.eta()-event.pfMET_p4.eta());
                vars["deta_hbb-met"] = (bb.eta()-event.pfMET_p4.eta());
                vars["abs(deta_hbb-htautau)"] = std::abs(bb.eta()-event.SVfit_p4.eta());
                vars["deta_hbb-hatutau"] = bb.eta()-event.SVfit_p4.eta();

                vars["dR_l1+l2"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2));
                vars["dR_b1+b2"] = (ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]));
                vars["dR_l1-met"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.pfMET_p4));
                vars["dR_l2-met"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_2, event.pfMET_p4));
                vars["dR_l1+l2-met"] = (ROOT::Math::VectorUtil::DeltaR(leptons, event.pfMET_p4));
                vars["dR_htautau-met"] = (ROOT::Math::VectorUtil::DeltaR(event.SVfit_p4, event.pfMET_p4));
                vars["dR_hbb-met"] = (ROOT::Math::VectorUtil::DeltaR(bb, event.pfMET_p4));
                vars["dR_hbb-htautau"] = (ROOT::Math::VectorUtil::DeltaR(bb, event.SVfit_p4));

                vars["dR_b1+b2Pt_hbb"] = (ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]))*bb.Pt();
                vars["dR_l1+l2Pt_htautau"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2))*event.SVfit_p4.Pt();

                vars["mass_l1+l2-met"] = ROOT::Math::VectorUtil::InvariantMass(leptons,event.pfMET_p4); //
                vars["mass_htautau"] = event.SVfit_p4.M();
                vars["mass_l1+l2"] = std::sqrt(pow(event.p4_1.Et()+event.p4_2.Et(),2)-pow(event.p4_1.px()+event.p4_2.px(),2)+pow(event.p4_1.py()+event.p4_2.py(),2));//
                vars["mass_hbb"] = bb.M();
                vars["MT_l1"] = Calculate_MT(event.p4_1,event.pfMET_p4);
                vars["MT_l2"] = Calculate_MT(event.p4_2,event.pfMET_p4);
                vars["MT_hatutau"] = Calculate_MT(event.SVfit_p4, event.pfMET_p4);
                vars["MT_l1+l2"] = Calculate_MT(leptons, event.pfMET_p4);
                vars["MT_tot"] = Calculate_TotalMT(event.p4_1,event.p4_2,event.pfMET_p4); //Total transverse mass
                vars["mass_H"] = ROOT::Math::VectorUtil::InvariantMass(bb,event.SVfit_p4);

                LorentzVectorM_Float t1 = event.p4_1 + event.jets_p4[0] + event.pfMET_p4;
                LorentzVectorM_Float t2 = event.p4_2 + event.jets_p4[0] + event.pfMET_p4;
                LorentzVectorM_Float t3 = event.p4_1 + event.jets_p4[1] + event.pfMET_p4;
                LorentzVectorM_Float t4 = event.p4_2 + event.jets_p4[1] + event.pfMET_p4;
                double d1 = std::abs(t1.mass() - mass_top);
                double d2 = std::abs(t2.mass() - mass_top);
                double d3 = std::abs(t3.mass() - mass_top);
                double d4 = std::abs(t4.mass() - mass_top);
                if (d1<d2 && d1<d3 && d1<d4) vars["mass_top(met)"] = d1;
                if (d2<d1 && d2<d3 && d2<d4) vars["mass_top(met)"] = d2;
                if (d3<d1 && d3<d2 && d3<d4) vars["mass_top(met)"] = d3;
                if (d4<d1 && d4<d3 && d4<d2) vars["mass_top(met)"] = d4;

                LorentzVectorM_Float t11 = event.p4_1 + event.jets_p4[0];
                LorentzVectorM_Float t22 = event.p4_2 + event.jets_p4[0];
                LorentzVectorM_Float t33 = event.p4_1 + event.jets_p4[1];
                LorentzVectorM_Float t44 = event.p4_2 + event.jets_p4[1];
                double d11 = std::abs(t11.mass() - mass_top);
                double d22 = std::abs(t22.mass() - mass_top);
                double d33 = std::abs(t33.mass() - mass_top);
                double d44 = std::abs(t44.mass() - mass_top);
                if (d11<d22 && d11<d33 && d11<d44) vars["mass_top"] = d11;
                if (d22<d11 && d22<d33 && d22<d44) vars["mass_top"] = d22;
                if (d33<d11 && d33<d22 && d33<d44) vars["mass_top"] = d33;
                if (d44<d11 && d44<d33 && d44<d22) vars["mass_top"] = d44;

                const analysis::LorentzVectorXYZ sv(event.SVfit_p4.px(),event.SVfit_p4.py(),event.SVfit_p4.pz(),event.SVfit_p4.e());
                const auto boosted_tau1 = ROOT::Math::VectorUtil::boost(event.p4_1, sv.BoostToCM());
                vars["theta_l1(h)"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_tau1, sv)); //theta angle between the first final state lepton and the direction of flight of h_tautau in the h_tautau rest frame
                vars["phi_l1(h)"] = std::atan(boosted_tau1.py()/boosted_tau1.px()); //phi angle between the first final state lepton and the direction of flight of h_tautau in the h_tautau rest frame
                const auto boosted_tau2 = ROOT::Math::VectorUtil::boost(event.p4_2, sv.BoostToCM());
                vars["theta_l2(h)"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_tau2, sv)); //angle between the second final state lepton and the direction of flight of h_tautau in the h_tautau rest frame
                vars["phi_l2(h)"] = std::atan(boosted_tau2.py()/boosted_tau2.px()); //phi angle between the second final state lepton and the direction of flight of h_tautau in the h_tautau rest frame
                vars["R_l1-l2(h)"] = ROOT::Math::VectorUtil::DeltaR(boosted_tau1, boosted_tau2); // R between the two final state leptons in the h_tautau rest frame

                const analysis::LorentzVectorXYZ hbb(bb.px(),bb.py(),bb.pz(),bb.e());
                const auto boosted_b1 = ROOT::Math::VectorUtil::boost(event.jets_p4[0], hbb.BoostToCM());
                vars["theta_b1(h)"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_b1, hbb)); //angle between the first final state bjet and the direction of flight of h_bb in the h_bb rest frame
                vars["phi_b1(h)"] = std::atan(boosted_b1.py()/boosted_b1.px()); //phi angle between the first final state bjet and the direction of flight of h_bb in the h_bb rest frame
                const auto boosted_b2 = ROOT::Math::VectorUtil::boost(event.jets_p4[2], hbb.BoostToCM());
                vars["theta_b2(h)"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_b2, hbb)); //angle between the second final state bjet and the direction of flight of h_bb in the h_bb rest frame
                if (boosted_b2.px()!=0) vars["phi_b2(h)"] = std::atan(boosted_b2.py()/boosted_b2.px()); //phi angle between the second final state bjet and the direction of flight of h_bb in the h_bb rest frame
                vars["R_b1-b2(h)"] = ROOT::Math::VectorUtil::DeltaR(boosted_b1, boosted_b2); // R between the two final state b-jetsin the h_bb rest frame

                LorentzVectorE_Float H = bb + event.SVfit_p4;
                const analysis::LorentzVectorXYZ vec_H(H.px(),H.py(),H.pz(),H.e());
                const auto boosted_l1 = ROOT::Math::VectorUtil::boost(event.p4_1, vec_H.BoostToCM());
                const auto boosted_l2 = ROOT::Math::VectorUtil::boost(event.p4_2, vec_H.BoostToCM());
                const auto boosted_j1 = ROOT::Math::VectorUtil::boost(event.jets_p4[0], vec_H.BoostToCM());
                const auto boosted_j2 = ROOT::Math::VectorUtil::boost(event.jets_p4[1], vec_H.BoostToCM());
                const TVector3 vec_l1(boosted_l1.px(),boosted_l1.py(),boosted_l1.pz());
                const TVector3 vec_l2(boosted_l2.px(),boosted_l2.py(),boosted_l2.pz());
                const TVector3 vec_j1(boosted_j1.px(),boosted_j1.py(),boosted_j1.pz());
                const TVector3 vec_j2(boosted_j2.px(),boosted_j2.py(),boosted_j2.pz());
                const auto n1 = vec_l1.Cross(vec_l2);
                const auto n2 = vec_j1.Cross(vec_j2);
                vars["phi(H)"] = ROOT::Math::VectorUtil::Angle(n1, n2); //angle between the decay planes of the four final state elements expressed in the H rest frame

                const auto boosted_htautau = ROOT::Math::VectorUtil::boost(event.SVfit_p4, vec_H.BoostToCM());
                vars["theta_star1(H)"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_htautau, ROOT::Math::Cartesian3D<>(0, 0, 1))); // Is the production angle of the h_tautau defined in the H rest frame

                const auto boosted_hbb = ROOT::Math::VectorUtil::boost(bb, vec_H.BoostToCM());
                vars["theta_star2(H)"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_hbb, ROOT::Math::Cartesian3D<>(0, 0, 1)));// Is the production angle of the h_bb defined in the H rest frame

                vars["phi_hbb-htautau(H)"] = ROOT::Math::VectorUtil::DeltaPhi(boosted_htautau,boosted_hbb); //Phi angle between hbb and htautau in the H rest frame
                vars["theta_hbb-htautau(H)"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_htautau, boosted_hbb)); //Theta angle between hbb and htautau in the H rest frame

                const TVector3 vec_htautau(boosted_htautau.px(),boosted_htautau.py(),boosted_htautau.pz());
                TVector3 z_axis(0,0,1);
                const auto n3 = vec_htautau.Cross(z_axis);
                vars["phi_1(H)"] = ROOT::Math::VectorUtil::Angle(n1,n3); //Angle between the decay plane of the lepton pair and a plane defined by the vector of the h_tautau in the H rest frame and the positive direction of z axis

                const TVector3 vec_hbb(boosted_hbb.px(),boosted_hbb.py(),boosted_hbb.pz());
                const auto n4 = vec_hbb.Cross(z_axis);
                vars["phi_2(H)"] = ROOT::Math::VectorUtil::Angle(n2,n4); //Angle between the decay plane of the b-jets pair and a plane defined by the vector of the h_bb in the H rest frame and the positive direction of z axis

                vars["phi_htautau(H)"] = std::atan(boosted_htautau.py()/boosted_htautau.px()); //Phi angle of h_tautau in the H rest frame
                vars["phi_hbb(H)"] = std::atan(vec_hbb.Y()/vec_hbb.X()); //Phi angle of h_bb in the H rest frame

                vars.AddEvent(samplename,1); //vars.AddEvent(samplename,1) to loop over all the masses, vars.AddEvent(samplename,mass) to look at every mass individually
            }
//            i++;
        }

        std::map<int, std::map<std::string,std::vector<double>>> sample_vars_signal = vars.GetSampleVariables("Signal");
        std::map<int, std::map<std::string,std::vector<double>>> sample_vars_bkg = vars.GetSampleVariables("Background");
        std::cout<<"n.variabili: "<<sample_vars_bkg[1].size()<<std::endl;
        std::cout<<"n.masse: "<<sample_vars_signal.size()<<" + "<<sample_vars_bkg.size()<<" fondo."<<std::endl;
        for (const auto& var: sample_vars_signal){
            if (var.second.size() == 0) continue;
            int mass = var.first;
            auto map = var.second;
            auto vector = map["pt_l1"];
            std::cout<<"massa: "<<mass<<", eventi segnale: "<<vector.size()<<std::endl;
        }
        for (const auto& var: sample_vars_bkg){
            if (var.second.size() == 0) continue;
            auto map = var.second;
            auto vector = map["pt_l1"];
            std::cout<<"eventi fondo: "<<vector.size()<<std::endl;
        }

        std::map<int, std::map<std::string, std::pair<TH1D*,TH1D*>>> histos = GetHistos(sample_vars_signal, sample_vars_bkg, outfile);

        std::map<int, std::map<VarPair,double>> cov_matrix_signal, cov_matrix_bkg;
        cov_matrix_signal = Cov(sample_vars_signal);
        cov_matrix_bkg = Cov(sample_vars_bkg);
        std::map<int, std::map<VarPair,double>> corr_matrix_signal, corr_matrix_bkg;
        corr_matrix_signal = CovToCorr(cov_matrix_signal);
        corr_matrix_bkg = CovToCorr(cov_matrix_bkg);
        CreateMatrixHistos(sample_vars_signal,corr_matrix_signal,outfile,"correlation","Signal");
        CreateMatrixHistos(sample_vars_bkg,corr_matrix_bkg,outfile,"correlation","Background");
        std::map<int, std::map<VarPair,double>> corr_matrix_signal_fix, corr_matrix_bkg_fix;
        corr_matrix_signal_fix = RemoveDiagonal(corr_matrix_signal);
        corr_matrix_bkg_fix = RemoveDiagonal(corr_matrix_bkg);

        std::map<int,std::vector< std::pair<VarPair, double>>> corr_vector_signal, corr_vector_bkg, corr_vector_difference;
        for (const auto& var: sample_vars_signal){
            int mass = var.first;
            std::vector<std::pair<VarPair, double>> cvs(corr_matrix_signal_fix[mass].begin(), corr_matrix_signal_fix[mass].end());
            std::vector<std::pair<VarPair, double>> cvb(corr_matrix_bkg_fix[1].begin(), corr_matrix_bkg_fix[1].end());
            std::vector<std::pair<VarPair, double>> cvd = Difference(cvs,cvb);
            std::sort(cvd.begin(), cvd.end(), [](const std::pair<VarPair, double>& a, const std::pair<VarPair, double>& b) {
                return std::abs(a.second) < std::abs(b.second);
            });
            corr_vector_signal[mass] = cvs;
            corr_vector_bkg[mass] = cvb;
            corr_vector_difference[mass] = cvd;
        }

        std::map<int, std::map<std::string,double>> kolmogorov = KolmogorovTest(sample_vars_signal, sample_vars_bkg);
        std::map< int, std::map<std::string,double>> chi2 = Chi2Test(histos);
       KolmogorvPlot(kolmogorov,outfile);



//        std::map<std::string, int> eliminate;
//        std::ofstream ofs(args.tree_name()+".txt", std::ofstream::out);
//        ofs << "Variable    " << "  Variable2   " << "  corr_s  " << "  corr_b    " << "  corr_d    "<<"  chi2_1    "<<"  chi2_2    "<<"  KS_1    "<<"  KS_2"<<std::endl;
//        auto corr_vec_diff = corr_vector_difference[1];
//        auto corr_mat_signal = corr_matrix_signal_fix[1];
//        auto corr_mat_bkg = corr_matrix_bkg_fix[1];
//        auto kolmo = kolmogorov[1];
//        auto chi = chi2[1];
//        for(Long64_t i = 0; i < corr_vec_diff.size(); i++){
//            std::pair<VarPair, double> element_difference = corr_vec_diff[i];
//            VarPair pair = element_difference.first;
//            if ((std::abs(corr_mat_signal[pair])>0.5 && std::abs(corr_mat_bkg[pair])>0.5) && element_difference.second<0.3){
//                ofs <<pair.first<<"    "<<pair.second<<"    "<<std::abs(corr_mat_signal[pair])<<"    "<<std::abs(corr_mat_bkg[pair])<<"    "<<element_difference.second<<"   "<<chi[pair.first]<<"    "<<chi[pair.second]<<"  "<<kolmo[pair.first]<<"    "<<kolmo[pair.second]<<std::endl;
//                if (kolmo[pair.second]<=kolmo[pair.first] && chi[pair.first]<=chi[pair.second]) eliminate[pair.first]++;
//                if (kolmo[pair.second]>kolmo[pair.first] && chi[pair.first]>chi[pair.second]) eliminate[pair.second]++;
//            }
//        }
//        ofs.close();
//        std::ofstream off("elimination.txt", std::ofstream::out);
//        for(const auto& i : eliminate){
//            off << i.first << "     " << i.second<<std::endl;
//        }
//        off.close();
    }
private:
    Arguments args;
    SampleEntryCollection samples;
    std::shared_ptr<TFile> outfile;
    MvaVariables vars;
};
}

PROGRAM_MAIN(analysis::MvaClassification, Arguments) // definition of the main program function









