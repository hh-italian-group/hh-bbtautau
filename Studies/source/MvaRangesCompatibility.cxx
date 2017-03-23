/*! Study of compatibility of range masses for each variable
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/StatEstimators.h"
#include "hh-bbtautau/Analysis/include/Lester_mt2_bisect.h"
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

template<typename LVector1, typename LVector2, typename LVector3, typename LVector4, typename LVector5 >
double Calculate_MT2(const LVector1& lepton1_p4, const LVector2& lepton2_p4, const LVector3& bjet_1, const LVector4& bjet_2, const LVector5& met_p4)
{
    asymm_mt2_lester_bisect::disableCopyrightMessage();
    const double mVisA = (lepton1_p4+bjet_1).mass();
    const double pxA = (lepton1_p4+bjet_1).px();
    const double pyA = (lepton1_p4+bjet_1).py();
    const double mVisB = (lepton2_p4+bjet_2).mass();
    const double pxB = (lepton2_p4+bjet_2).px();
    const double pyB = (lepton2_p4+bjet_2).py();
    const double pxMet = met_p4.px();
    const double pyMet = met_p4.py();
    double chiA = 0.; // hypothesised mass of invisible on side A.  Must be >=0.
    double chiB = 0.; // hypothesised mass of invisible on side B.  Must be >=0.
    double MT2 =  asymm_mt2_lester_bisect::get_mT2(mVisA, pxA, pyA,mVisB, pxB, pyB,pxMet, pyMet,chiA, chiB,0);
    return MT2;
}


//Create a pair of histogram for each variable. The first one is for signal, the second for background.
std::map<int, std::map<std::string, std::pair<TH1D*,TH1D*>>> GetHistos(std::map<int, std::map<std::string,std::vector<double>>> sample_signal, std::map<int, std::map<std::string,std::vector<double>>> sample_bkg, std::shared_ptr<TFile> outfile){

    std::map<int, std::map<std::string, std::pair<TH1D*,TH1D*>>> histogram;
    for(const auto& var_m : sample_signal){
        std::map<std::string, std::pair<TH1D*,TH1D*>> histos;
        std::map<std::string,std::vector<double>> v_b = sample_bkg[1];
        std::string mass;
        if (var_m.first != 1) mass = "_mass"+std::to_string(var_m.first);
        else mass = "";
        if (var_m.first == 2000) mass = "_SM";
        for(const auto& var_1 : var_m.second) {
            std::vector<double> vector_s = var_1.second;
            std::vector<double> vector_b = v_b[var_1.first];

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
            TH1D* h_s = new TH1D((var_1.first+"Signal"+mass).c_str(),(var_1.first+"Signal"+mass).c_str(),nbin,min-1,max+1);
            h_s->SetXTitle((var_1.first).c_str());
            for(Long64_t i = 0; i < vector_s.size(); i++){
                h_s->Fill(vector_s[i]);
            }
            h_s->Scale(1/h_s->Integral());
            root_ext::WriteObject(*h_s, outfile.get());

            TH1D* h_b = new TH1D((var_1.first+"Bkg"+mass).c_str(),(var_1.first+"Bkg"+mass).c_str(),nbin,min-1,max+1);
            h_b->SetXTitle((var_1.first).c_str());
            for(Long64_t i = 0; i < vector_b.size(); i++){
                h_b->Fill(vector_b[i]);
            }
            h_b->Scale(1/h_b->Integral());
            root_ext::WriteObject(*h_b, outfile.get());
            h_b->Delete();
         }
        histogram[var_m.first] = histos;
    }
    return histogram;
}

//Optimal Bandwidth
std::map<int, std::map<std::string, double>> OptimalBandwidth(std::map<int, std::map<std::string,std::vector<double>>> sample_signal){
    std::map<int, std::map<std::string, double>> bandwidth;
    for(const auto& var_mass : sample_signal) {
        if (var_mass.second.size() == 0) continue;
        std::map<std::string, double> bandwidth_int;
        for (const auto& var : var_mass.second){
            bandwidth_int[var.first] = stat_estimators::OptimalBandwith(var.second);
        }
        bandwidth[var_mass.first] = bandwidth_int;
    }
    return bandwidth;
}

//Kolmogorov test
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
            k[var.first] = TMath::KolmogorovTest(vector_signal.size(), v_s, vector_bkg.size(), v_b, "M");
        }
        kolmogorov[var_mass.first] = k;
    }
    return kolmogorov;

}

//Create plots of Kolmogorv probability between signal and background for different values of mass
void KolmogorovPlotSignalBkg(std::map<int, std::map<std::string,double>> kolmogorov, std::shared_ptr<TFile> outfile){
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
            plot[var.first]->SetLineColor(kGreen+1);
            plot[var.first]->SetLineWidth(1);
            plot[var.first]->SetMarkerColor(1);
            plot[var.first]->SetMarkerSize(1);
            plot[var.first]->SetMarkerStyle(3);
            plot[var.first]->SetTitle(("Kolmogorov_"+ var.first+"_SignalBkg").c_str());
            plot[var.first]->SetName(("Kolmogorov_"+ var.first+"_SignalBkg").c_str());
            plot[var.first]->GetHistogram()->GetXaxis()->SetTitle("mass");
            plot[var.first]->GetHistogram()->GetYaxis()->SetTitle("KS Distance");
        }
        i++;
    }
    for(const auto& var: map_mass){
        root_ext::WriteObject(*plot[var.first], outfile.get());
    }
}

//Create plots of Jensen Shannon Divergence Divergence between signal and background for different values of mass
void KLDPlotSignalBkg(std::map<int, std::map<std::string,std::vector<double>>> sample_signal, std::map<int, std::map<std::string,std::vector<double>>> sample_bkg, std::map<int, std::map<std::string, double>> band_signal, std::map<int, std::map<std::string, double>> band_bkg, std::shared_ptr<TFile> outfile){
    std::map<std::string, TGraph*> plot;
    auto map_bkg = sample_bkg[1];
    int i = 0;
    auto bandwidth_bkg = band_bkg[1];
    for(const auto& var_mass_2 : sample_signal){
        if (var_mass_2.second.size() == 0) continue;
        auto map_signal_2 = var_mass_2.second;
        auto bandwidth_signal = band_signal[var_mass_2.first];
        for (const auto& var : var_mass_2.second){
            double k = stat_estimators::JensenShannonDivergence(map_bkg[var.first], map_signal_2[var.first], bandwidth_bkg[var.first],bandwidth_signal[var.first]);
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
        i++;
    }
    for(const auto& var: map_bkg){
        root_ext::WriteObject(*plot[var.first], outfile.get());
    }
}

//Create plots of Kolmogorv probability for compatibility of signal at different mass
void CompatibilitySignalPlotKolmogorov(std::map<int, std::map<std::string,std::vector<double>>> sample_signal, std::shared_ptr<TFile> outfile){
    std::map<std::string, TGraph*> plot;

    auto var_mass = *sample_signal.begin();
    auto map_signal = var_mass.second;
    int i = 0;
    for(const auto& var_mass_2 : sample_signal){
        if (var_mass_2.second.size() == 0) continue;
        auto map_signal_2 = var_mass_2.second;
        for (const auto& var : var_mass.second){
            std::vector<double> vector_signal = map_signal[var.first];
            std::sort(vector_signal.begin(), vector_signal.end());
            std::vector<double> vector_signal_2 = map_signal_2[var.first];
            std::sort(vector_signal_2.begin(), vector_signal_2.end());
            Double_t v_s[vector_signal.size()], v_s_2[vector_signal_2.size()];
            for(Long64_t i = 0; i < vector_signal.size(); i++){
                v_s[i] = vector_signal[i];
            }
            for(Long64_t i = 0; i < vector_signal_2.size(); i++){
                v_s_2[i] = vector_signal_2[i];
            }
            double k = TMath::KolmogorovTest(vector_signal.size(), v_s, vector_signal_2.size(), v_s_2, "");
            if (plot.count(var.first)==0) plot[var.first] = new TGraph();
            plot[var.first]->SetPoint(i,var_mass_2.first,k);
            plot[var.first]->SetLineColor(kGreen+1);
            plot[var.first]->SetLineWidth(1);
            plot[var.first]->SetMarkerColor(1);
            plot[var.first]->SetMarkerSize(1);
            plot[var.first]->SetMarkerStyle(3);
            plot[var.first]->SetTitle(("Kolmogorov_"+var.first+"_Signal").c_str());
            plot[var.first]->SetName(("Kolmogorov_"+var.first+"_Signal").c_str());
            plot[var.first]->GetHistogram()->GetXaxis()->SetTitle("mass");
            plot[var.first]->GetHistogram()->GetYaxis()->SetTitle("KS Probability");
        }
        i++;
    }
    for(const auto& var: map_signal){
        root_ext::WriteObject(*plot[var.first], outfile.get());
    }
}

//Create plots of Jensen Shannon Divergence for compatibility of signal at different mass
void CompatibilitySignalPlotKLD(std::map<int, std::map<std::string,std::vector<double>>> sample_signal, std::map<int, std::map<std::string, double>> bandwidth, std::shared_ptr<TFile> outfile){
    std::map<std::string, TGraph*> plot;

    auto var_mass = *sample_signal.begin();
    auto bandwidth_signal_1 = bandwidth[var_mass.first];
    auto map_signal = var_mass.second;
    int i = 0;
    for(const auto& var_mass_2 : sample_signal){
        if (var_mass_2.second.size() == 0) continue;
        auto map_signal_2 = var_mass_2.second;
        auto bandwidth_signal_2 = bandwidth[var_mass_2.first];
        for (const auto& var : var_mass.second){
            double k = stat_estimators::JensenShannonDivergence(map_signal_2[var.first], map_signal[var.first], bandwidth_signal_2[var.first], bandwidth_signal_1[var.first]);
            if (plot.count(var.first)==0) plot[var.first] = new TGraph();
            plot[var.first]->SetPoint(i,var_mass_2.first,k);
            plot[var.first]->SetLineColor(kGreen+1);
            plot[var.first]->SetLineWidth(1);
            plot[var.first]->SetMarkerColor(1);
            plot[var.first]->SetMarkerSize(1);
            plot[var.first]->SetMarkerStyle(3);
            plot[var.first]->SetTitle(("JSD_"+var.first+"_Signal").c_str());
            plot[var.first]->SetName(("JSD_"+var.first+"_Signal").c_str());
            plot[var.first]->GetHistogram()->GetXaxis()->SetTitle("mass");
            plot[var.first]->GetHistogram()->GetYaxis()->SetTitle("JS Divergence");
        }
        i++;
    }
    for(const auto& var: map_signal){
        root_ext::WriteObject(*plot[var.first], outfile.get());
    }
}

//Create 2Dhisto of Kolmogorov Distance for compatibility of signal at different mass
std::map<std::string,TH2D*> CompatibilitySignalHistoKolmogorov(std::map<int, std::map<std::string,std::vector<double>>> sample_signal, std::shared_ptr<TFile> outfile){

    int bin = sample_signal.size();
    std::map<std::string,TH2D*> map_histo;
    int k = 1;
    auto var_mass = *sample_signal.begin();
    for (const auto& var : var_mass.second){
         TH2D* histo= new TH2D(("Kolmogorov"+var.first+"CompatibilitySignal").c_str(),("Kolmogorov"+var.first+"CompatibilitySignal").c_str(),bin,0,0,bin,0,0);
         histo->SetCanExtend(TH1::kAllAxes);
         histo->SetBins(bin, 0, bin, bin, 0, bin);
         histo->GetXaxis()->SetTitle("mass");
         histo->GetYaxis()->SetTitle("mass");
         map_histo[var.first] = histo;
    }
    for(const auto& var_mass : sample_signal){
        for (const auto& var : var_mass.second){
            if (var_mass.first == 2000 ) {
                map_histo[var.first]->GetXaxis()->SetBinLabel(k, "SM");
                map_histo[var.first]->GetYaxis()->SetBinLabel(k, "SM");
            }
            else{
                map_histo[var.first]->GetXaxis()->SetBinLabel(k, (std::to_string(var_mass.first)).c_str());
                map_histo[var.first]->GetYaxis()->SetBinLabel(k, (std::to_string(var_mass.first)).c_str());
            }
        }
        k++;
    }
    int i = 1;
    for(const auto& var_mass : sample_signal){
        auto map_signal = var_mass.second;
        if (map_signal.size() == 0) continue;
        int j = 1;
        for(const auto& var_mass_2 : sample_signal){

            auto map_signal_2 = var_mass_2.second;
            if (map_signal_2.size() == 0) continue;
            for (const auto& var : var_mass.second){
                std::vector<double> vector_signal = map_signal[var.first];
                std::sort(vector_signal.begin(), vector_signal.end());
                std::vector<double> vector_signal_2 = map_signal_2[var.first];
                std::sort(vector_signal_2.begin(), vector_signal_2.end());
                Double_t v_s[vector_signal.size()], v_s_2[vector_signal_2.size()];
                for(Long64_t i = 0; i < vector_signal.size(); i++){
                    v_s[i] = vector_signal[i];
                }
                for(Long64_t i = 0; i < vector_signal_2.size(); i++){
                    v_s_2[i] = vector_signal_2[i];
                }
                double k = TMath::KolmogorovTest(vector_signal.size(), v_s, vector_signal_2.size(), v_s_2, "");
                map_histo[var.first]->SetBinContent(i, j, k);
            }
            j++;
        }
        i++;
    }
    for(const auto& var: var_mass.second){
        root_ext::WriteObject(*map_histo[var.first], outfile.get());
    }
    return map_histo;
}

//Create 2Dhisto of Jensen Shannon Divergence Divergence for compatibility of signal at different mass
std::map<std::string,TH2D*> CompatibilitySignalHistoKLD(std::map<int, std::map<std::string,std::vector<double>>> sample_signal, std::map<int, std::map<std::string, double>> bandwidth, std::shared_ptr<TFile> outfile){
    int bin = sample_signal.size();
    std::map<std::string,TH2D*> map_histo;
    int k = 1;
    auto var_mass = *sample_signal.begin();
    for (const auto& var : var_mass.second){
         TH2D* histo= new TH2D(("JSD_"+var.first+"CompatibilitySignal").c_str(),("JSD_"+var.first+"CompatibilitySignal").c_str(),bin,0,0,bin,0,0);
         histo->SetCanExtend(TH1::kAllAxes);
         histo->SetBins(bin, 0, bin, bin, 0, bin);
         histo->GetXaxis()->SetTitle("mass");
         histo->GetYaxis()->SetTitle("mass");
         map_histo[var.first] = histo;
    }
    for(const auto& var_mass : sample_signal){
        for (const auto& var : var_mass.second){
            if (var_mass.first == 2000 ) {
                map_histo[var.first]->GetXaxis()->SetBinLabel(k, "SM");
                map_histo[var.first]->GetYaxis()->SetBinLabel(k, "SM");
            }
            else{
                map_histo[var.first]->GetXaxis()->SetBinLabel(k, (std::to_string(var_mass.first)).c_str());
                map_histo[var.first]->GetYaxis()->SetBinLabel(k, (std::to_string(var_mass.first)).c_str());
            }
        }
        k++;
    }
    int i = 1;
    for(const auto& var_mass : sample_signal){
        auto map_signal = var_mass.second;
        if (map_signal.size() == 0) continue;
        int j = 1;
        auto bandwidth_signal_1 = bandwidth[var_mass.first];
        for(const auto& var_mass_2 : sample_signal){
            auto map_signal_2 = var_mass_2.second;
            if (map_signal_2.size() == 0) continue;
            auto bandwidth_signal_2 = bandwidth[var_mass_2.first];
            for (const auto& var : var_mass.second){
            double k = stat_estimators::JensenShannonDivergence(map_signal_2[var.first],map_signal[var.first], bandwidth_signal_1[var.first],bandwidth_signal_2[var.first]);
            map_histo[var.first]->SetBinContent(i, j, k);
            }
            j++;
        }
        i++;
    }
    for(const auto& var: var_mass.second){
        root_ext::WriteObject(*map_histo[var.first], outfile.get());
    }
    return map_histo;
}


//For each variable study ranges compatibility and for each of them compute the kolmogorov test between signal and background
void CompatibilityRangeKolmogorov(std::map<std::string,TH2D*> map_histo, std::map<int, std::map<std::string,std::vector<double>>> sample_signal, std::map<int, std::map<std::string,std::vector<double>>> sample_bkg, std::map<int, std::map<std::string, double>> bandwidth_bkg, std::shared_ptr<TFile> outfile){

    sample_signal.erase(--sample_signal.rbegin().base());
    double masse[sample_signal.size()];
    int i = 0;
    for (const auto& var : sample_signal){
        if (i<(sample_signal.size())) masse[i] = var.first;
        i++;
    }

    std::map<std::string, double> band_bkg = bandwidth_bkg[1];
    auto map_bkg = sample_bkg[1];
    auto var_mass = *sample_signal.begin();

    struct element{
        std::map<int, std::vector<int>> ranges;
        std::map<int, double> kolmogorov;
        std::map<int, double> kullback;
        std::map<int, double> entropy;
        std::map<int, double> kullback_ss;
    };
    std::vector<std::pair<std::string, element>> vector_pair;
    int count = 0;

    TH1D* histo_nrange =  new TH1D("Number of ranges", "Number of ranges",20,0,20);
    TH1D* histo_end_1 =  new TH1D("End of first range", "End of first range",sample_signal.size(),0,1000);
    TH1D* histo_end_2 =  new TH1D("End of second range", "End of second range",sample_signal.size(),0,1000);
    histo_end_1->SetBins(sample_signal.size()-1,masse);
    histo_end_2->SetBins(sample_signal.size()-1,masse);
    histo_nrange->SetXTitle("# range");
    histo_end_1->SetXTitle("mass");
    histo_end_2->SetXTitle("mass");


    for (const auto& var : var_mass.second){
        element el;
        std::vector<int> range;
        int j;
        for (Long64_t i = 1; i <=(sample_signal.size()); i++){
            j = i;
            double bin_content = 1;
            int k = 0;
            while(bin_content>0.05 && j<=(sample_signal.size()) && bin_content!=0){
                j++;
                bin_content = map_histo[var.first]->GetBinContent(i,j);
                k++;
            }
            range.push_back(k);
            i = j-1;
        }
        auto vec_bkg = map_bkg[var.first];
        double vector_bkg[vec_bkg.size()];
        for (Long64_t i = 1; i <vec_bkg.size(); i++) vector_bkg[i] = vec_bkg[i];

        int z = 0;
        std::map<int, std::vector<int>> map_ranges;
        std::map<int, double> map_kolmogorov;
        std::map<int, double> map_kullback;
        std::map<int, double> map_kullback_ss;
        std::map<int, double> map_entropy;
        std::map<int, double> map_bandwidth;
        std::map<int,std::vector<double>> total_mass;
        for (Long64_t i = 0; i <range.size(); i++){
            std::vector<int> vec_masse;
            std::vector<double> p_sum;
            int k = range[i];
            int j = 0;
            while(j<k){
                int mass = masse[z];
                vec_masse.push_back(mass);
                auto map_signal = sample_signal[mass];
                auto vec_signal = map_signal[var.first];
                std::copy(std::begin(vec_signal),  std::end(vec_signal), std::back_inserter(p_sum));
                z++;
                j++;
            }
            map_ranges[i] = vec_masse;
            total_mass[i] = p_sum;

            double vector_signal[p_sum.size()];
            for (Long64_t i = 1; i <p_sum.size(); i++) vector_signal[i] = p_sum[i];
            map_kolmogorov[i] = TMath::KolmogorovTest(p_sum.size(), vector_signal, vec_bkg.size(), vector_bkg, "M");
            map_bandwidth[i] = stat_estimators::OptimalBandwith(p_sum,0.01);
            double kull = stat_estimators::JensenShannonDivergence(vec_bkg, total_mass[i], map_bandwidth[i],band_bkg[var.first]);
            map_entropy[i] = stat_estimators::Entropy(total_mass[i], map_bandwidth[i]);
            map_kullback[i] = kull;
            if (i>0) {
                double kull_ss = stat_estimators::JensenShannonDivergence(total_mass[i-1], p_sum, map_bandwidth[i],map_bandwidth[i-1]);
                map_kullback_ss[i] = kull_ss;
            }
        }
        el.entropy = map_entropy;
        el.kullback = map_kullback;
        el.kullback_ss = map_kullback_ss;
        el.kolmogorov = map_kolmogorov;
        el.ranges = map_ranges;
        vector_pair.emplace_back(var.first, el);
        count++;
    }

    std::sort(vector_pair.begin(), vector_pair.end(), [](std::pair<std::string, element>& el1,  std::pair<std::string, element>& el2 ){
        auto first = el1.second;
        auto second = el2.second;
        return first.ranges.size()<second.ranges.size();
    });

    std::ofstream ofs("compatibility_no250.csv", std::ofstream::out);
    for (Long64_t i = 0; i<vector_pair.size(); i++){

        auto element = vector_pair[i].second;
        std::map<int, std::vector<int>> element_range = element.ranges;
        std::map<int, double> element_kolmogorov = element.kolmogorov;
        std::map<int, double> element_kullback = element.kullback;
        std::map<int, double> element_kullback_ss = element.kullback_ss;
        std::map<int, double> element_entropy = element.entropy;
        ofs<<vector_pair[i].first<<","<<" n.range:"<<element_range.size()<<std::endl;
        ofs<<",";
        histo_nrange->Fill(element_range.size());
        for (Long64_t j=1; j<=element_range.size();j++){
            ofs<<"Range"<<j<<",";
        }
        ofs<<std::endl;
        ofs<<"Mass: "<<",";
        int j = 0;
        for (const auto& var : element_range){
            auto vector_element_range = element_range[var.first];
            if (vector_element_range.size()>1) ofs<<vector_element_range.front()<<"-"<<vector_element_range.back()<<",";
            else if (vector_element_range.size()==1) ofs<<vector_element_range.front()<<",";
            if (j == 0) histo_end_1->Fill(vector_element_range.back());
            if (j == 1) histo_end_2->Fill(vector_element_range.back());
            j++;
        }

        ofs<<std::endl;
        ofs<<"KS-distance: ";
        for (const auto& var : element_range){
            double k = element_kolmogorov[var.first];
            ofs<<","<<k;
        }
        ofs<<std::endl;
        ofs<<"JSD: ";
        for (const auto& var : element_range){
            double k = element_kullback[var.first];
            ofs<<","<<k;
        }
        ofs<<std::endl;
        ofs<<"Entropy: ";
        for (const auto& var : element_range){
            double k = element_entropy[var.first];
            ofs<<","<<k;
        }
        std::string name = vector_pair[i].first;
        double s = stat_estimators::Entropy(map_bkg[name], band_bkg[name]);
        ofs<<","<<s;
        ofs<<std::endl;
        ofs<<"JSD-prev: ";
        for (const auto& var : element_range){
            double k = element_kullback_ss[var.first];
            ofs<<","<<k;
        }
        ofs<<std::endl;
    }
    root_ext::WriteObject(*histo_nrange, outfile.get());
    root_ext::WriteObject(*histo_end_1, outfile.get());
    root_ext::WriteObject(*histo_end_2, outfile.get());
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
        int m[20] = {250,260,270,280,300,320,340,350,400,450,500,550,600,650,700,750,800,900,2000,1};
        int i=0;
        for(const SampleEntry& entry:samples)
        {
            int mass = m[i];
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

                vars["pt_l1"] = event.p4_1.pt();
                vars["pt_l2"] = event.p4_2.pt();
                vars["pt_b1"] = event.jets_p4[0].pt();
                vars["pt_b2"] = event.jets_p4[1].pt();
                vars["pt_l1l2"] = leptons.pt();
                vars["pt_htautau"] = event.SVfit_p4.pt();
                vars["pt_l1l2met"] = leptonsMET.pt();
                vars["pt_hbb"] = bb.pt();
                vars["pt_met"] = event.pfMET_p4.pt();

                vars["p_zeta"] = Calculate_Pzeta(event.p4_1, event.p4_2, event.pfMET_p4);
                vars["p_zetavisible"] = Calculate_visiblePzeta(event.p4_1,event.p4_2);

                vars["abs(dphi_l1l2)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.p4_2));
                vars["dphi_l1l2"] = (ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.p4_2));
                vars["abs(dphi_b1b2)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.jets_p4[0], event.jets_p4[1]));
                vars["dphi_b1b2"] = (ROOT::Math::VectorUtil::DeltaPhi(event.jets_p4[0], event.jets_p4[1]));
                vars["abs(dphi_l1met)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.pfMET_p4));
                vars["dphi_l1met"] = (ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.pfMET_p4));
                vars["abs(dphi_l2met)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_2, event.pfMET_p4));
                vars["dphi_l2met"] = (ROOT::Math::VectorUtil::DeltaPhi(event.p4_2, event.pfMET_p4));
                vars["abs(dphi_l1l2met)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(leptons, event.pfMET_p4));
                vars["dphi_l1l2met"] = (ROOT::Math::VectorUtil::DeltaPhi(leptons, event.pfMET_p4));
                vars["abs(dphi_htautaumet)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4));
                vars["dphi_htautaumet"] = (ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4));
                vars["abs(dphi_hbbmet)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.pfMET_p4));
                vars["dphi_hbbmet"] = (ROOT::Math::VectorUtil::DeltaPhi(bb, event.pfMET_p4));
                vars["abs(dphi_hbbhatutau)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.SVfit_p4));
                vars["dphi_hbbhtautau"] = (ROOT::Math::VectorUtil::DeltaPhi(bb, event.SVfit_p4));

                vars["abs(deta_l1l2)"] = std::abs(event.p4_1.eta() - event.p4_2.eta());
                vars["deta_l1l2"] = (event.p4_1.eta() - event.p4_2.eta());
                vars["abs(deta_b1b2)"] = std::abs(event.jets_p4[0].eta() - event.jets_p4[1].eta());
                vars["deta_b1b2"] = (event.jets_p4[0].eta() - event.jets_p4[1].eta());
                vars["abs(deta_l1met)"] = std::abs(event.p4_1.eta()-event.pfMET_p4.eta());
                vars["deta_l1met"] = (event.p4_1.eta()-event.pfMET_p4.eta());
                vars["abs(deta_l2met)"] = std::abs(event.p4_2.eta()-event.pfMET_p4.eta());
                vars["deta_l2met"] = (event.p4_2.eta()-event.pfMET_p4.eta());
                vars["abs(deta_l1l2met)"] = std::abs(leptons.eta()-event.pfMET_p4.eta());
                vars["deta_l1l2met"] = (leptons.eta()-event.pfMET_p4.eta());
                vars["abs(deta_htautaumet)"] = std::abs(event.SVfit_p4.eta()-event.pfMET_p4.eta());
                vars["deta_hatutaumet"] = (event.SVfit_p4.eta()-event.pfMET_p4.eta());
                vars["abs(deta_hbbmet)"] = std::abs(bb.eta()-event.pfMET_p4.eta());
                vars["deta_hbbmet"] = (bb.eta()-event.pfMET_p4.eta());
                vars["abs(deta_hbbhtautau)"] = std::abs(bb.eta()-event.SVfit_p4.eta());
                vars["deta_hbbhatutau"] = bb.eta()-event.SVfit_p4.eta();

                vars["dR_l1l2"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2));
                vars["dR_b1b2"] = (ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]));
                vars["dR_l1met"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.pfMET_p4));
                vars["dR_l2met"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_2, event.pfMET_p4));
                vars["dR_l1l2met"] = (ROOT::Math::VectorUtil::DeltaR(leptons, event.pfMET_p4));
                vars["dR_htautaumet"] = (ROOT::Math::VectorUtil::DeltaR(event.SVfit_p4, event.pfMET_p4));
                vars["dR_hbbmet"] = (ROOT::Math::VectorUtil::DeltaR(bb, event.pfMET_p4));
                vars["dR_hbbhtautau"] = (ROOT::Math::VectorUtil::DeltaR(bb, event.SVfit_p4));

                vars["dR_b1b2Pt_hbb"] = (ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]))*bb.Pt();
                vars["dR_l1l2Pt_htautau"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2))*event.SVfit_p4.Pt();

                vars["mass_l1l2met"] = ROOT::Math::VectorUtil::InvariantMass(leptons,event.pfMET_p4); //
                vars["mass_htautau"] = event.SVfit_p4.M();
                vars["mass_l1l2"] = std::sqrt(pow(event.p4_1.Et()+event.p4_2.Et(),2)-pow(event.p4_1.px()+event.p4_2.px(),2)+pow(event.p4_1.py()+event.p4_2.py(),2));//
                vars["mass_hbb"] = bb.M();
                vars["MT_l1"] = Calculate_MT(event.p4_1,event.pfMET_p4);
                vars["MT_l2"] = Calculate_MT(event.p4_2,event.pfMET_p4);
                vars["MT_hatutau"] = Calculate_MT(event.SVfit_p4, event.pfMET_p4);
                vars["MT_l1l2"] = Calculate_MT(leptons, event.pfMET_p4);
                vars["MT_tot"] = Calculate_TotalMT(event.p4_1,event.p4_2,event.pfMET_p4); //Total transverse mass
                vars["MT2"] = std::min(Calculate_MT2(event.p4_1,event.p4_2,event.jets_p4[0], event.jets_p4[1], event.pfMET_p4),Calculate_MT2(event.p4_1, event.jets_p4[1], event.p4_2,event.jets_p4[0], event.pfMET_p4)); //Stransverse mass
                vars["mass_H"] = ROOT::Math::VectorUtil::InvariantMass(bb,event.SVfit_p4);

                LorentzVectorM_Float a1 = event.p4_1 + event.jets_p4[0] + event.pfMET_p4;
                LorentzVectorM_Float b1 = event.p4_2 + event.jets_p4[1];
                LorentzVectorM_Float a2 = event.p4_1 + event.jets_p4[0];
                LorentzVectorM_Float b2 = event.p4_2 + event.jets_p4[1] + event.pfMET_p4;
                LorentzVectorM_Float a3 = event.p4_1 + event.jets_p4[1] + event.pfMET_p4;
                LorentzVectorM_Float b3 = event.p4_2 + event.jets_p4[0];
                LorentzVectorM_Float a4 = event.p4_1 + event.jets_p4[1];
                LorentzVectorM_Float b4 = event.p4_2 + event.jets_p4[0] + event.pfMET_p4;

                double d1 = pow(std::abs(a1.mass() - mass_top),2) + pow (std::abs(b1.mass() - mass_top),2);
                double d2 = pow(std::abs(a2.mass() - mass_top),2) + pow (std::abs(b2.mass() - mass_top),2);
                double d3 = pow(std::abs(a3.mass() - mass_top),2) + pow (std::abs(b3.mass() - mass_top),2);
                double d4 = pow(std::abs(a4.mass() - mass_top),2) + pow (std::abs(b4.mass() - mass_top),2);

                if (d1<d2 && d1<d3 && d1<d4) {
                    vars["Mass_top1"] = a1.mass();
                    vars["Mass_top2"] = b1.mass();
                }
                if (d2<d1 && d2<d3 && d2<d4) {
                    vars["Mass_top1"] = a2.mass();
                    vars["Mass_top2"] = b2.mass();
                }
                if (d3<d1 && d3<d2 && d3<d4) {
                    vars["Mass_top1"] = a3.mass();
                    vars["Mass_top2"] = b3.mass();
                }
                if (d4<d1 && d4<d3 && d4<d2) {
                    vars["Mass_top1"] = a4.mass();
                    vars["Mass_top2"] = b4.mass();
                }

                const analysis::LorentzVectorXYZ sv(event.SVfit_p4.px(),event.SVfit_p4.py(),event.SVfit_p4.pz(),event.SVfit_p4.e());
                const auto boosted_tau1 = ROOT::Math::VectorUtil::boost(event.p4_1, sv.BoostToCM());
//                vars["theta_l1"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_tau1, sv)); //theta angle between the first final state lepton and the direction of flight of h_tautau in the h_tautau rest frame
//                vars["phi_l1"] = std::atan(boosted_tau1.py()/boosted_tau1.px()); //phi angle between the first final state lepton and the direction of flight of h_tautau in the h_tautau rest frame
                const auto boosted_tau2 = ROOT::Math::VectorUtil::boost(event.p4_2, sv.BoostToCM());
//                vars["theta_l2"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_tau2, sv)); //angle between the second final state lepton and the direction of flight of h_tautau in the h_tautau rest frame
//                vars["phi_l2"] = std::atan(boosted_tau2.py()/boosted_tau2.px()); //phi angle between the second final state lepton and the direction of flight of h_tautau in the h_tautau rest frame
                vars["R_l1l2"] = ROOT::Math::VectorUtil::DeltaR(boosted_tau1, boosted_tau2); // R between the two final state leptons in the h_tautau rest frame

                const analysis::LorentzVectorXYZ hbb(bb.px(),bb.py(),bb.pz(),bb.e());
                const auto boosted_b1 = ROOT::Math::VectorUtil::boost(event.jets_p4[0], hbb.BoostToCM());
//                vars["theta_b1"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_b1, hbb)); //angle between the first final state bjet and the direction of flight of h_bb in the h_bb rest frame
//                vars["phi_b1"] = std::atan(boosted_b1.py()/boosted_b1.px()); //phi angle between the first final state bjet and the direction of flight of h_bb in the h_bb rest frame
                const auto boosted_b2 = ROOT::Math::VectorUtil::boost(event.jets_p4[2], hbb.BoostToCM());
//                vars["theta_b2"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_b2, hbb)); //angle between the second final state bjet and the direction of flight of h_bb in the h_bb rest frame
//                if (boosted_b2.px()!=0) vars["phi_b2(h)"] = std::atan(boosted_b2.py()/boosted_b2.px()); //phi angle between the second final state bjet and the direction of flight of h_bb in the h_bb rest frame
                vars["R_b1b2"] = ROOT::Math::VectorUtil::DeltaR(boosted_b1, boosted_b2); // R between the two final state b-jetsin the h_bb rest frame

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
                vars["phi"] = ROOT::Math::VectorUtil::Angle(n1, n2); //angle between the decay planes of the four final state elements expressed in the H rest frame

                const auto boosted_htautau = ROOT::Math::VectorUtil::boost(event.SVfit_p4, vec_H.BoostToCM());
                vars["theta_star1"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_htautau, ROOT::Math::Cartesian3D<>(0, 0, 1))); // Is the production angle of the h_tautau defined in the H rest frame

                const auto boosted_hbb = ROOT::Math::VectorUtil::boost(bb, vec_H.BoostToCM());
                vars["theta_star2"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_hbb, ROOT::Math::Cartesian3D<>(0, 0, 1)));// Is the production angle of the h_bb defined in the H rest frame

                const TVector3 vec_htautau(boosted_htautau.px(),boosted_htautau.py(),boosted_htautau.pz());
                TVector3 z_axis(0,0,1);
                const auto n3 = vec_htautau.Cross(z_axis);
                vars["phi_1"] = ROOT::Math::VectorUtil::Angle(n1,n3); //Angle between the decay plane of the lepton pair and a plane defined by the vector of the h_tautau in the H rest frame and the positive direction of z axis

                const TVector3 vec_hbb(boosted_hbb.px(),boosted_hbb.py(),boosted_hbb.pz());
                const auto n4 = vec_hbb.Cross(z_axis);
                vars["phi_2"] = ROOT::Math::VectorUtil::Angle(n2,n4); //Angle between the decay plane of the b-jets pair and a plane defined by the vector of the h_bb in the H rest frame and the positive direction of z axis

//                vars["phi_htautau(H)"] = std::atan(boosted_htautau.py()/boosted_htautau.px()); //Phi angle of h_tautau in the H rest frame
//                vars["phi_hbb(H)"] = std::atan(boosted_hbb.py()/boosted_hbb.px()); //Phi angle of h_bb in the H rest frame

                vars.AddEvent(samplename,mass); //vars.AddEvent(samplename,1) to loop over all the masses, vars.AddEvent(samplename,mass) to look at every mass individually
            }
            i++;
        }

        std::map<int, std::map<std::string,std::vector<double>>> sample_vars_signal = vars.GetSampleVariables("Signal");
        sample_vars_signal.erase(250);
        std::map<int, std::map<std::string,std::vector<double>>> sample_vars_bkg = vars.GetSampleVariables("Background");
        std::cout<<"n.variabili: "<<sample_vars_bkg[1].size()<<std::endl;
        std::cout<<"n.masse: "<<sample_vars_signal.size()<<" + "<<sample_vars_bkg.size()<<" fondo."<<std::endl;
        for (const auto& var: sample_vars_signal){
            if (var.second.size() == 0) continue;
            int mass = var.first;
            auto map = var.second;
            auto vector = map["pt_l1"];
            if (mass != 1 && mass != 2000) {
                std::cout<<"massa: "<<mass<<", eventi segnale: "<<vector.size()<<std::endl;
            }
            else {
                if (mass == 2000) std::cout<<"eventi SM: "<<vector.size()<<std::endl;
                else std::cout<<"eventi segnale: "<<vector.size()<<std::endl;
            }
        }

        std::map<int, std::map<std::string, double>> bandwith_signal = OptimalBandwidth(sample_vars_signal);
        std::map<int, std::map<std::string, double>> bandwidth_bkg = OptimalBandwidth(sample_vars_bkg);

        for (const auto& var: sample_vars_bkg){
            if (var.second.size() == 0) continue;
            auto map = var.second;
            auto vector = map["pt_l1"];
            std::cout<<"eventi fondo: "<<vector.size()<<std::endl;
        }

        std::map<int, std::map<std::string, std::pair<TH1D*,TH1D*>>> histos = GetHistos(sample_vars_signal, sample_vars_bkg, outfile);

        std::map<int, std::map<std::string,double>> kolmogorov = KolmogorovTest(sample_vars_signal, sample_vars_bkg);
        KolmogorovPlotSignalBkg(kolmogorov, outfile);
        KLDPlotSignalBkg(sample_vars_signal,sample_vars_bkg,bandwith_signal,bandwidth_bkg,outfile);
        std::map<std::string,TH2D*> map_histo_kolmogorov = CompatibilitySignalHistoKolmogorov(sample_vars_signal,outfile);
        std::map<std::string,TH2D*> map_histo_KLD = CompatibilitySignalHistoKLD(sample_vars_signal,bandwith_signal,outfile);
        CompatibilityRangeKolmogorov(map_histo_kolmogorov,sample_vars_signal,sample_vars_bkg, bandwidth_bkg, outfile);
        CompatibilitySignalPlotKolmogorov(sample_vars_signal,outfile);
        CompatibilitySignalPlotKLD(sample_vars_signal, bandwith_signal,outfile);

    }
private:
    Arguments args;
    SampleEntryCollection samples;
    std::shared_ptr<TFile> outfile;
    MvaVariables vars;
};
}

PROGRAM_MAIN(analysis::MvaClassification, Arguments) // definition of the main program function
//./run.sh MvaVariableCorrelation --input_path ~/Desktop/tuples --output_file ranges.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name muTau --number_events 50000
