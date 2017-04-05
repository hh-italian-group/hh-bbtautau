/*! Study of correlation matrix and mutual information of BDT variables
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <TMatrixDEigen.h>
#include <TH1.h>
#include <TH2.h>
#include <fstream>
#include <random>
#include <algorithm>
#include <functional>
#include <future>
#include <initializer_list>
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
#include "hh-bbtautau/Analysis/include/MvaVariables.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(Long64_t, number_events);
    REQ_ARG(Long64_t, number_variables);
};

#define MY_TREE_DATA() \
    VAR(Float_t, var) \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(ntuple, MyEvent, MyTuple, MY_TREE_DATA, "my")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(ntuple, MyTuple, MY_TREE_DATA)
#undef VAR
#undef MY_TREE_DATA


namespace analysis {

//using namespace std::placeholders;

struct SampleEntry{
  std::string filename;
  double weight;
  bool issignal;
  int mass;
  SampleEntry() : weight(-1), issignal(false){}
};

std::ostream& operator<<(std::ostream& os, const SampleEntry& entry)
{
    os << entry.filename << " " << entry.mass << " " << entry.weight << " " << std::boolalpha << entry.issignal;
    return os;
}

std::istream& operator>>(std::istream& is, SampleEntry& entry)
{
    is >> entry.filename >> entry.mass >> entry.weight >> std::boolalpha >> entry.issignal;
    return is;
}

using DataVector = std::vector<double>;
using VarData = std::map<std::string, DataVector>;
using MassVar = std::map<int, VarData>;

class MvaVariablesStudy : public MvaVariables {
private:
    std::map<std::string, MassVar> all_variables;
    std::map<std::string, double> variables;

public:
    using MvaVariables::MvaVariables;

    virtual void SetValue(const std::string& name, double value) override
    {
        variables[name] = value;
    }

    virtual void AddEventVariables(bool issignal, bool istraining, int mass, double weight) override
    {
        const std::string sample = issignal ? "Signal" : "Background";
        MassVar& sample_map = all_variables[sample];
        VarData& sample_vars = sample_map[mass];
        for(const auto& name_value : variables) {
            const std::string& name = name_value.first;
            const double value = name_value.second;
            sample_vars[name].push_back(value);
        }
    }

    const MassVar& GetSampleVariables(const std::string& sample) const
    {
        return all_variables.at(sample);
    }

};

struct Name_ND{
    std::set<std::string> names;
    bool operator<(const Name_ND& x) const
    {
        if (names.size() != x.names.size()) return names.size()<x.names.size();
        if (names.size() == 0) return false;
        auto x_iter = x.names.begin();
        for(auto iter = names.begin(); iter != names.end(); ++iter, ++x_iter){
            if(*iter != *x_iter) return *iter < *x_iter;
        }
        return false;
    }
    Name_ND(std::initializer_list<std::string> _name):names(_name.begin(), _name.end()){}

    bool IsSubset(const std::set<std::string>& other) const
    {
        std::set<std::string> difference;
        std::set_difference(names.begin(),names.end(), other.begin(), other.end(), std::inserter(difference, difference.begin()));
        return !difference.size();
    }
};

using VarPair = std::pair<std::string,std::string>;

//Create a pair of histogram for each variable and for each mass. The first one is for signal, the second for background.
std::map<int, std::map<std::string, std::pair<TH1D*,TH1D*>>> GetHistos(std::map<int, std::map<std::string,std::vector<double>>> sample_signal,  std::map<int, std::map<std::string,std::vector<double>>> sample_bkg, std::shared_ptr<TFile> outfile){
    std::map<int, std::map<std::string, std::pair<TH1D*,TH1D*>>> histogram;
    for(const auto& mass_entry : sample_signal){
        std::map<std::string, std::pair<TH1D*,TH1D*>>& histos = histogram[mass_entry.first];
        std::map<std::string,std::vector<double>> v_b = sample_bkg[1];
        std::string mass;
        if (mass_entry.first != 1) mass = "_mass"+std::to_string(mass_entry.first);
        else mass = "";
        if (mass_entry.first == 2000) mass = "_SM";
        for(const auto& map_entry : mass_entry.second) {
            std::vector<double> vector_s = map_entry.second;
            std::vector<double> vector_b = v_b[map_entry.first];
            const int nbin = 50;
            std::pair<TH1D*,TH1D*>& h = histos[map_entry.first];
            h.first = new TH1D((map_entry.first+"Signal"+mass).c_str(),(map_entry.first+"Signal"+mass).c_str(),nbin,0,0);
            h.first->SetCanExtend(TH1::kAllAxes);
            h.first->SetXTitle((map_entry.first).c_str());
            for(Long64_t i = 0; i < vector_s.size(); i++){
                h.first->Fill(vector_s[i]);
            }
            root_ext::WriteObject(*h.first, outfile.get());
            h.first->Delete();

            h.second = new TH1D((map_entry.first+"Bkg"+mass).c_str(),(map_entry.first+"Bkg"+mass).c_str(),nbin,0,0);
            h.second->SetCanExtend(TH1::kAllAxes);
            h.second->SetXTitle((map_entry.first).c_str());
            for(Long64_t i = 0; i < vector_b.size(); i++){
                h.second->Fill(vector_b[i]);
            }
            root_ext::WriteObject(*h.second, outfile.get());
            h.second->Delete();
         }
    }
    return histogram;
}

using VarBand = std::map<std::string, double>;
using VarPairMI =  std::map<VarPair,double>;
using MassVarBand = std::map<int, VarBand>;
using MassVarPairMI = std::map<int, VarPairMI>;

//Creat optimal bandwidth for a single value of mass
VarBand OptimalBandwidth(const VarData& sample){
    VarBand bandwidth;
    std::map<std::string, std::future<double>> bandwidth_future;
        for (const auto& var : sample){
            bandwidth_future[var.first] = std::async(std::launch::async,stat_estimators::OptimalBandwith<double>,var.second, 0.01);
        }
        for(auto& var : bandwidth_future) {
            if(!var.second.valid())
                throw exception("future not valid");
            bandwidth[var.first] = var.second.get();
        }
    return bandwidth;
}

//Create elements of mutual information matrix for a single value of mass
VarPairMI Mutual(const VarData& sample, const VarBand& bandwidth){
    VarPairMI mutual_matrix;
    std::map<VarPair, std::future<double>> mutual_matrix_future;
    for (const auto& var_1: sample){
        for(const auto& var_2 : sample) {
            const VarPair var_12(var_1.first, var_2.first);
            const VarPair var_21(var_2.first, var_1.first);
            if(mutual_matrix_future.count(var_21)) continue;
            mutual_matrix_future[var_12] = std::async(std::launch::async,stat_estimators::ScaledMutualInformation<double>, var_1.second, var_2.second, bandwidth.at(var_1.first), bandwidth.at(var_2.first));
        }
    }
    for(auto& pair_entry : mutual_matrix_future) {
        if(!pair_entry.second.valid())
            throw exception("future not valid");
        mutual_matrix[pair_entry.first] = pair_entry.second.get();
    }
    return mutual_matrix;
}

void MutualPlot(int mass, const VarPairMI& mutual_matrix_signal,
                const VarPairMI& mutual_matrix_bkg, std::shared_ptr<TFile> outfile){
    auto plot = std::make_shared<TH2D>(("MutualInformation2D_mass"+std::to_string(mass)).c_str(),("MutualInformation2D"+std::to_string(mass)).c_str(),20,0,1,20,0,1);
    for (const auto& entry : mutual_matrix_signal){
        plot->Fill(entry.second, mutual_matrix_bkg.at(entry.first));
        plot->SetXTitle("MI Signal");
        plot->SetYTitle("MI Bkg");
    }
    root_ext::WriteObject(*plot, outfile.get());
}


std::map<int, std::map<VarPair,stat_estimators::EstimatedQuantity>> EstimateCovariance(const std::map<int, std::set<std::string>>& selected, const std::map<int, std::map<std::string,std::vector<double>>>& sample_vars, uint_fast32_t seed){
    std::map<int, std::map<VarPair,stat_estimators::EstimatedQuantity>> covariance_matrix;
    for (const auto& var: sample_vars){
        if (var.second.size() == 0) continue;
        for(const auto& var_1 : selected.at(250)) {
            for(const auto& var_2 : selected.at(250)) {
                const VarPair var_12(var_1, var_2);
                const VarPair var_21(var_2, var_1);
                if(covariance_matrix[var.first].count(var_21)) continue;
                covariance_matrix[var.first][var_12] = stat_estimators::EstimateWithErrorsByResampling(stat_estimators::Covariance<double>,sample_vars.at(var.first).at(var_1),sample_vars.at(var.first).at(var_2),true,true,1000,0.31731, seed);
            }
        }
    }
    return covariance_matrix;
}

std::map<int, std::map<VarPair,double>> GetValue(const std::map<int, std::map<VarPair,stat_estimators::EstimatedQuantity>>& cov_matrix_signal){
    std::map<int, std::map<VarPair,double>> map;
    for (const auto& var: cov_matrix_signal){
        for(const auto& var_1 : var.second) {
            auto map_cov_varpair = cov_matrix_signal.at(var.first).at(var_1.first);
            double value = map_cov_varpair.value;
            map[var.first][var_1.first] = value;
        }
    }
    return map;
}

std::map<int, std::map<VarPair,double>> CovToCorr(const std::map<int, std::map<VarPair,double>>& cov_matrix){
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

//Create covariance(correlation/mutual information) histo matrix
void CreateMatrixHistos(const std::map<int, std::map<std::string,std::vector<double>>>& sample_vars, const std::map<int, std::map<VarPair,double>>& element, std::shared_ptr<TFile> outfile, std::string type,  std::string class_sample){
        for(const auto var_m: sample_vars){
            if (var_m.second.size() == 0) continue;
            int bin =  var_m.second.size();
            auto el = element.at(var_m.first);
            std::string mass= "";
            if (var_m.first != 1 && var_m.first != 2000) mass = "_mass"+std::to_string(var_m.first);
            else if (var_m.first == 2000) mass = "_SM";
            auto matrix = std::make_shared<TH2D>((type+"_"+class_sample+mass).c_str(),(type+"_"+class_sample+mass).c_str(),bin,0,0,bin,0,0);
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
                        matrix->SetBinContent(i, j, el.at(var_12));
                        matrix->SetBinContent(j, i, el.at(var_12));
                    }
                    j++;
                }
                i++;
            }
            root_ext::WriteObject(*matrix, outfile.get());
        }
    }

//Create covariance(correlation/mutual information) histo matrix
void CreateMatrixHistosSelected(const std::map<int, std::map<std::string,std::vector<double>>>& sample_vars, const std::map<int, std::set<std::string>>& selected, const MassVarPairMI& element, std::shared_ptr<TFile> outfile,
                        const std::string& type, const std::string& class_sample){
    for(const auto var_m: sample_vars){
        if (var_m.second.size() == 0) continue;
        int bin =  var_m.second.size();
        std::string mass;
        if (var_m.first != 1 && var_m.first != 2000) mass = "_mass"+std::to_string(var_m.first);
        else if (var_m.first == 2000) mass = "_SM";
        auto matrix = std::make_shared<TH2D>((type+"_"+class_sample+mass).c_str(),(type+"_"+class_sample+mass).c_str(),bin,0,bin,bin,0,bin);
        int i = 1;
        for(const auto& var_1 : selected.at(250)) {
            int j = 1;
            matrix->GetXaxis()->SetBinLabel(i, (var_1).c_str());
            matrix->GetYaxis()->SetBinLabel(i, (var_1).c_str());
            for(const auto& var_2 : selected.at(250)) {
                if (j >= i){
                const VarPair var_12(var_1, var_2);
                matrix->SetBinContent(i, j, (int)(element.at(var_m.first).at(var_12)*100));
                matrix->SetBinContent(j, i, (int)(element.at(var_m.first).at(var_12)*100));
                }
                j++;
            }
            i++;
        }
        root_ext::WriteObject(*matrix, outfile.get());
    }
}


using Name_NDDistance = std::map<Name_ND, double>;
using MassName_ND = std::map<int, Name_NDDistance>;
using VectorName_ND = std::deque<std::pair<Name_ND, double>>;
using MassVectorName_ND = std::map<int, VectorName_ND>;

Name_NDDistance JensenDivergenceND(const VarData& sample_signal, const VarData& sample_bkg,
                                   const VarBand& bandwidth_signal, const VarBand& bandwidth_bkg){
    Name_NDDistance  JSDivergenceND;
    std::map<Name_ND,std::future<double>> JSDivergenceND_future;
    for (const auto& entry_1 : sample_signal){
        std::vector<const DataVector*> x;
        std::vector<const DataVector*> y;
        DataVector band_x;
        DataVector band_y;
        x.push_back(&sample_signal.at(entry_1.first));
        y.push_back(&sample_bkg.at(entry_1.first));
        band_x.push_back(bandwidth_signal.at(entry_1.first));
        band_y.push_back(bandwidth_bkg.at(entry_1.first));
        JSDivergenceND_future[Name_ND{entry_1.first}] = std::async(std::launch::async,stat_estimators::JensenShannonDivergence_ND<double>, x, y, band_x, band_y);
        for (const auto& entry_2 : sample_signal){
            if ( entry_1.first == entry_2.first ) continue;
            if ( JSDivergenceND_future.count(Name_ND({entry_2.first,entry_1.first})) ) continue;
            x.push_back(&sample_signal.at(entry_2.first));
            y.push_back(&sample_bkg.at(entry_2.first));
            band_x.push_back(bandwidth_signal.at(entry_2.first));
            band_y.push_back(bandwidth_bkg.at(entry_2.first));
            JSDivergenceND_future[Name_ND{entry_1.first,entry_2.first}] = std::async(std::launch::async,stat_estimators::JensenShannonDivergence_ND<double>, x, y, band_x, band_y);
            x.erase(x.end() - 1);
            y.erase(y.end() - 1);
            band_x.erase(band_x.end() - 1);
            band_y.erase(band_y.end() - 1);
        }
    }

    for(auto& entry : JSDivergenceND_future) {
        JSDivergenceND[entry.first] = entry.second.get();
    }
    return JSDivergenceND;
}

std::set<std::string> CheckList(const Name_NDDistance& JSDivergenceND,const VarPairMI& mutual_matrix_signal,
               const VarPairMI& mutual_matrix_bkg, const Long64_t& num_variables){

    std::set<std::string> selected, checked;
    VectorName_ND  JSDivergence_vector(JSDivergenceND.begin(), JSDivergenceND.end());
    while(selected.size() < num_variables) {
        std::sort(JSDivergence_vector.begin(), JSDivergence_vector.end(), [](const std::pair<Name_ND, double>& el1, const std::pair<Name_ND, double>& el2){
            return el1.second > el2.second;
        });
        const auto& best_entry = JSDivergence_vector.front();
        while(JSDivergence_vector.size() && best_entry.first.IsSubset(checked))
            JSDivergence_vector.pop_front();
        if(!JSDivergence_vector.size()) break;
        for(const auto& name : best_entry.first.names) {
            checked.insert(name);
            if(selected.count(name)) continue;

            const double JS_1d = JSDivergenceND.at(Name_ND({name}));
            for(auto& entry : JSDivergence_vector) {
                if(entry.first.names.count(name))
                    entry.second -= JS_1d;
            }

            bool passed = true;
            for(const auto& selected_name : selected) {
                static const double threashold = 0.5;
                const VarPair var_pair = name < selected_name ?
                            VarPair(name, selected_name) : VarPair(selected_name, name);
                const double MI_sgn = mutual_matrix_signal.at(var_pair);
                const double MI_bkg = mutual_matrix_bkg.at(var_pair);
                if(MI_sgn < threashold && MI_bkg < threashold) {
                    passed = false;
                    break;
                }
            }
            if(passed)
                selected.insert(name);
        }
    }
    return selected;
}

void Plot2D(const std::map<int, std::set<std::string>>& selected, const MassVar& sample_signal,
            const  MassVar& sample_bkg, std::shared_ptr<TFile> outfile){
    std::set<VarPair> checked;
    for (const auto& entry : selected){
        for (const auto & set_entry_1 :  entry.second){
            for (const auto & set_entry_2 :  entry.second){
                checked.insert(VarPair(set_entry_1, set_entry_2));
                if (checked.count(VarPair(set_entry_2,set_entry_1))) continue;
                auto histo_signal = std::make_shared<TH2D>((set_entry_1+set_entry_2+"Signal").c_str(), (set_entry_1+set_entry_2+"Signal").c_str(), 20,0,0, 20,0,0);
                histo_signal->SetCanExtend(TH1::kAllAxes);
                histo_signal->SetXTitle((set_entry_1).c_str());
                histo_signal->SetYTitle((set_entry_2).c_str());
                auto histo_bkg = std::make_shared<TH2D>((set_entry_1+set_entry_2+"Bkg").c_str(), (set_entry_1+set_entry_2+"Bkg").c_str(), 20,0,0, 20, 0, 0);
                histo_bkg->SetCanExtend(TH1::kAllAxes);
                histo_bkg->SetXTitle((set_entry_1).c_str());
                histo_bkg->SetYTitle((set_entry_2).c_str());
                for(Long64_t i = 0; i< sample_signal.at(entry.first).at(set_entry_1).size(); i++){
                    histo_signal->Fill(sample_signal.at(entry.first).at(set_entry_1)[i],sample_signal.at(entry.first).at(set_entry_2)[i]);
                }
                for(Long64_t i = 0; i< sample_bkg.at(1).at(set_entry_1).size(); i++){
                    histo_bkg->Fill(sample_bkg.at(1).at(set_entry_1)[i],sample_bkg.at(1).at(set_entry_2)[i]);
                }
                root_ext::WriteObject(*histo_signal, outfile.get());
                root_ext::WriteObject(*histo_bkg, outfile.get());
            }
        }
    }
}

void PlotJensenShannon(const MassName_ND& JSDivergenceND, std::shared_ptr<TFile> outfile){
    std::map<std::string, std::shared_ptr<TGraph>> plot;
    int i = 0;
    for (const auto& entry: JSDivergenceND){
        for (const auto& var : entry.second){
            if (var.first.names.size()!=1) continue;
            if (plot.count(*var.first.names.begin()) == 0) {
                plot[*var.first.names.begin()] = std::make_shared<TGraph>();
            }
            plot[*var.first.names.begin()]->SetPoint(i,entry.first,JSDivergenceND.at(entry.first).at(Name_ND({*var.first.names.begin()})));
            plot[*var.first.names.begin()]->SetLineColor(kGreen+1);
            plot[*var.first.names.begin()]->SetLineWidth(1);
            plot[*var.first.names.begin()]->SetMarkerColor(1);
            plot[*var.first.names.begin()]->SetMarkerSize(1);
            plot[*var.first.names.begin()]->SetMarkerStyle(3);
            plot[*var.first.names.begin()]->SetTitle(("JensenShannon_"+*var.first.names.begin()+"_SignalBkg").c_str());
            plot[*var.first.names.begin()]->SetName(("JensenShannon_"+*var.first.names.begin()+"_SignalBkg").c_str());
            plot[*var.first.names.begin()]->GetHistogram()->GetXaxis()->SetTitle("mass");
            plot[*var.first.names.begin()]->GetHistogram()->GetYaxis()->SetTitle("JSDivergence");
        }
        i++;
    }
    for(const auto& var: plot){
        root_ext::WriteObject(*var.second, outfile.get());
    }
}


using VecVariables =  std::vector<std::string>;
std::vector<VecVariables> VariableSelection( const MassVar& sample_signal, const  MassVar& sample_bkg,
                        const MassVarBand& band_signal, const MassVarBand& band_bkg,
                        const MassVarPairMI& mutual_matrix_signal,
                        const MassVarPairMI& mutual_matrix_bkg, std::map<int, std::set<std::string>>& selected,
                        std::shared_ptr<TFile> outfile, const Long64_t& num_variables){

    MassName_ND JSDivergenceND;
    for (const auto& sample: sample_signal){
        JSDivergenceND[sample.first] = JensenDivergenceND(sample_signal.at(sample.first), sample_bkg.at(1), band_signal.at(sample.first), band_bkg.at(1));        
    }
    PlotJensenShannon(JSDivergenceND, outfile);

    std::ofstream ofs("VariablesND.csv", std::ofstream::out);
    for(const auto& mass_entry : sample_signal){
        selected[mass_entry.first] =  CheckList(JSDivergenceND.at(mass_entry.first), mutual_matrix_signal.at(mass_entry.first), mutual_matrix_bkg.at(1), num_variables);
        ofs<<mass_entry.first<<std::endl;
        for(auto& entry_selected: selected[mass_entry.first]){
                ofs<<entry_selected<<",";
        }
        ofs<<std::endl;
        auto histo = std::make_shared<TH1D>(("JSDivergenceND"+std::to_string(mass_entry.first)).c_str(), ("JSDivergenceND"+std::to_string(mass_entry.first)).c_str(), 30,0,1);
        histo->SetXTitle("JSD");
        std::ofstream of("JSDivergenceND"+(std::to_string(mass_entry.first))+".csv", std::ofstream::out);
        of<<"var_1"<<","<<"var_2"<<","<<"JSDivergence"<<","<<"D_12-(D_1+D_2)"<<","<<"ScaledMIsignal_12" <<","<<"ScaledMIbkg_12"<<","<<"selected"<<std::endl;
        for(auto& entry : JSDivergenceND.at(mass_entry.first)){
             histo->Fill(entry.second);
             if ((entry.first).names.size() == 1){
                 auto name_1 = ((entry.first).names).begin();
                 bool selected_1 = selected[mass_entry.first].count(*name_1);
                 of<<*name_1<<","<<"   "<<","<<JSDivergenceND.at(mass_entry.first).at(Name_ND({*name_1}))<<","<<","<<","<<","<<selected_1<<std::endl;
             }
             else  if ((entry.first).names.size() == 2){
                 auto name_1 = ((entry.first).names).begin();
                 auto name_2 = ((entry.first).names).rbegin();
                 bool selected_1 = selected[mass_entry.first].count(*name_1);
                 bool selected_2 = selected[mass_entry.first].count(*name_2);
                 of<<*name_1<<","<<*name_2<<","<<JSDivergenceND.at(mass_entry.first).at(Name_ND({*name_1,*name_2}))<<","<<JSDivergenceND.at(mass_entry.first).at(Name_ND({*name_1,*name_2}))-(JSDivergenceND.at(mass_entry.first).at(Name_ND({*name_1}))+JSDivergenceND.at(mass_entry.first).at(Name_ND({*name_2})))<<","<<mutual_matrix_signal.at(mass_entry.first).at(VarPair(*name_1,*name_2))<<","<<mutual_matrix_bkg.at(1).at(VarPair(*name_1,*name_2))<<","<<selected_1<<","<<selected_2<<std::endl;
             }
        }
        root_ext::WriteObject(*histo, outfile.get());
        of.close();
    }
    ofs.close();

    int bin = selected.size();
    auto matrix_intersection = std::make_shared<TH2D>("IntersectionND","Intersection of variables",bin, 0, bin, bin, 0, bin);
    matrix_intersection->SetXTitle("mass");
    matrix_intersection->SetYTitle("mass");
    int k = 1;
    VecVariables union_1, union_2, union_3;
    for(const auto& mass_1: selected){
        int j = 1;
        std::string mass;
        if (mass_1.first != 1 && mass_1.first != 2000) mass = std::to_string(mass_1.first);
        else if (mass_1.first == 2000) mass = "SM";
        matrix_intersection->GetXaxis()->SetBinLabel(k, (mass).c_str());
        matrix_intersection->GetYaxis()->SetBinLabel(k, (mass).c_str());
        for(const auto& mass_2: selected){
            VecVariables intersection;
            std::set_intersection(mass_1.second.begin(), mass_1.second.end(), mass_2.second.begin(), mass_2.second.end(),
                                    std::back_inserter(intersection));
            matrix_intersection->SetBinContent(k, j, intersection.size());
            j++;

            if (mass_1.first<300 && mass_2.first<300){
                std::set_union(mass_1.second.begin(), mass_1.second.end(), mass_2.second.begin(), mass_2.second.end(),
                                                    std::back_inserter(union_1));
            }
            else if (mass_1.first<=450 && mass_2.first<=450){
                std::set_union(mass_1.second.begin(), mass_1.second.end(), mass_2.second.begin(), mass_2.second.end(),
                                                    std::back_inserter(union_2));
            }
            else if (mass_1.first<=900 && mass_2.first<=900){
                std::set_union(mass_1.second.begin(), mass_1.second.end(), mass_2.second.begin(), mass_2.second.end(),
                                                    std::back_inserter(union_3));
            }
        }
        k++;
    }
    root_ext::WriteObject(*matrix_intersection, outfile.get());

    Plot2D(selected, sample_signal, sample_bkg, outfile);
    sort(union_1.begin(),union_1.end() );
    union_1.erase( unique( union_1.begin(), union_1.end() ), union_1.end() );
    sort( union_2.begin(),union_2.end() );
    union_2.erase( unique( union_2.begin(), union_2.end() ), union_2.end() );
    sort( union_3.begin(),union_3.end() );
    union_3.erase( unique( union_3.begin(), union_3.end() ), union_3.end() );
    std::cout<<union_1.size()<<" "<<union_2.size()<<" "<<union_3.size()<<std::endl;

    std::vector<VecVariables> vector_union;
    vector_union.emplace_back(union_1);
    vector_union.emplace_back(union_2);
    vector_union.emplace_back(union_3);
    std::ofstream of("UnionVariables.csv", std::ofstream::out);
    of<<"Massa 250-280"<<std::endl;
    for(const auto& entry: union_1){
        of<<entry<<",";
    }
    of<<std::endl;
    of<<"Massa 300-450"<<std::endl;
    for(const auto& entry: union_2){
        of<<entry<<",";
    }
    of<<std::endl;
    of<<"Massa 500-900"<<std::endl;
    for(const auto& entry: union_3){
        of<<entry<<",";
    }
    of.close();
    return vector_union;
}


void KolmogorovSignalPlot(const MassVar& sample_signal, const VecVariables& selected,
                          int min, int max, int medium, std::shared_ptr<TFile> outfile){
    std::map<std::string, std::shared_ptr<TGraph>> plot;
     int i = 0;
    for (const auto& entry: sample_signal){
        if (entry.first < min || entry.first > max) continue;
        std::cout<<entry.first<<std::endl;
        for (const auto& var : selected){
            std::vector<double> vector_signal = sample_signal.at(entry.first).at(var);
            std::sort(vector_signal.begin(), vector_signal.end());
            std::vector<double> vector_signal_2 = sample_signal.at(medium).at(var);
            std::sort(vector_signal_2.begin(), vector_signal_2.end());
            Double_t* v_s = vector_signal.data(), *v_s_2 = vector_signal_2.data();
            double k = TMath::KolmogorovTest(vector_signal.size(), v_s, vector_signal_2.size(), v_s_2, "");
            if (plot.count(var) == 0) {
                plot[var] = std::make_shared<TGraph>();
                plot[var]->SetLineColor(kGreen+1);
                plot[var]->SetLineWidth(1);
                plot[var]->SetMarkerColor(1);
                plot[var]->SetMarkerSize(1);
                plot[var]->SetMarkerStyle(3);
                plot[var]->SetTitle(("Kolmogorov_"+var+"_Signal"+std::to_string(min)+"_"+std::to_string(max)).c_str());
                plot[var]->SetName(("Kolmogorov_"+var+"_Signal"+std::to_string(min)+"_"+std::to_string(max)).c_str());
                plot[var]->GetHistogram()->GetXaxis()->SetTitle("mass");
                plot[var]->GetHistogram()->GetYaxis()->SetTitle("KS Probability");
            }
            plot[var]->SetPoint(i,entry.first,k);
        }
        i++;
    }
    for(const auto& var: selected){
        root_ext::WriteObject(*plot[var], outfile.get());
    }
}


using SampleEntryCollection = std::vector<SampleEntry>;

class MvaClassification {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    static const std::set<std::string>& GetEnabledBranches()
    {

        static const std::set<std::string> EnabledBranches_read = {
            "eventEnergyScale", "q_1", "q_2", "jets_p4", "extraelec_veto", "extramuon_veto ", "SVfit_p4",
            "pfMET_p4", "p4_1", "p4_2"
        };
        return EnabledBranches_read;
    }

//    static const std::set<std::string>& GetDisabledBranches()
//    {
//        static const std::set<std::string> DisabledBranches_read = {
//            "dphi_mumet", "dphi_metsv", "dR_taumu", "mT1", "mT2", "dphi_bbmet", "dphi_bbsv", "dR_bb", "m_bb", "n_jets",
//            "btag_weight", "ttbar_weight",  "PU_weight", "shape_denominator_weight", "trigger_accepts", "trigger_matches",
//            "event.tauId_keys_1","event.tauId_keys_2","event.tauId_values_1","event.tauId_values_2"
//        };
//        return DisabledBranches_read;
//    }

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
        outfile(root_ext::CreateRootFile(args.output_file())), split_training_testing(false),
        vars(split_training_testing,UINT_FAST32_MAX, enabled_vars, disabled_vars)
    {
    }

    void Run()
    {

        for(const SampleEntry& entry:samples)
        {
            auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
            EventTuple tuple(args.tree_name(), input_file.get(), true, {} , GetEnabledBranches());
            std::cout << entry << " number of events: " << std::min(tuple.GetEntries(),args.number_events()) << std::endl;

            for(Long64_t current_entry = 0; current_entry < std::min(tuple.GetEntries(),args.number_events()); ++current_entry) {

                tuple.GetEntry(current_entry);
                const Event& event = tuple.data();

                if (event.eventEnergyScale != 0 || (event.q_1+event.q_2) != 0 || event.jets_p4.size() < 2
                    || event.extraelec_veto == true || event.extramuon_veto == true || event.jets_p4[0].eta() > 2.4
                    || event.jets_p4[1].eta() > 2.4) continue;

                LorentzVectorE_Float bb = event.jets_p4[0] + event.jets_p4[1];
                double ellipse_cut = pow(event.SVfit_p4.mass()-116,2)/pow(35.,2) + pow(bb.mass()-111,2)/pow(45.,2);
                if (ellipse_cut>1) continue;

                vars.AddEvent(event, entry.issignal, entry.mass, entry.weight);
            }
        }

        MassVar sample_vars_signal = vars.GetSampleVariables("Signal");
        MassVar sample_vars_bkg = vars.GetSampleVariables("Background");
        std::cout << "n.variabili: " << sample_vars_bkg[1].size() << std::endl;
        std::cout << "n.masse segnale: " << sample_vars_signal.size() << " + " << sample_vars_bkg.size() << " fondo."<<std::endl;

        for (const auto& var: sample_vars_signal){
            auto map = var.second;
            if (var.first != 1 && var.first != 2000) {
                std::cout << "massa: " << var.first << ", eventi segnale: " << map["pt_l1"].size() << std::endl;
            }
            else {
                if (var.first == 2000) std::cout<<"eventi SM: "<<map["pt_l1"].size()<<std::endl;
                else std::cout << "eventi segnale: "<< map["pt_l1"].size() << std::endl;
            }
        }
        for (const auto& var: sample_vars_bkg){
            if (var.second.size() == 0) continue;
            auto map = var.second;
            std::cout<<"eventi fondo: "<<map["pt_l1"].size()<<std::endl;
        }

        MassVarBand bandwidth_signal, bandwidth_bkg;
        MassVarPairMI mutual_matrix_signal, mutual_matrix_bkg;
        std::cout << "bandwidth and mutual information" << std::endl;
        for (const auto& sample: sample_vars_bkg){
            bandwidth_bkg[sample.first] = OptimalBandwidth(sample.second);
            mutual_matrix_bkg[sample.first] = Mutual(sample.second, bandwidth_bkg[sample.first]);
        }
        for (const auto& sample: sample_vars_signal){
            bandwidth_signal[sample.first] = OptimalBandwidth(sample.second);
            mutual_matrix_signal[sample.first] = Mutual(sample.second, bandwidth_signal[sample.first]);
            MutualPlot(sample.first, mutual_matrix_signal[sample.first], mutual_matrix_bkg[1], outfile);
        }
        CreateMatrixHistos(sample_vars_signal,mutual_matrix_signal,outfile,"mutual","Signal");
        CreateMatrixHistos(sample_vars_bkg,mutual_matrix_bkg,outfile,"mutual","Background");

        std::cout<<"selection variable"<<std::endl;
        std::map<int, std::set<std::string>> selected;
        std::vector<VecVariables> union_selected = VariableSelection(sample_vars_signal,sample_vars_bkg,bandwidth_signal,bandwidth_bkg,mutual_matrix_signal,mutual_matrix_bkg, selected, outfile, args.number_variables());
        std::cout<<"Kolmogorov"<<std::endl;
        KolmogorovSignalPlot(sample_vars_signal,union_selected[0],250,280,260,outfile);
        KolmogorovSignalPlot(sample_vars_signal,union_selected[1],300,450,350,outfile);
        KolmogorovSignalPlot(sample_vars_signal,union_selected[2],500,900,700,outfile);

        auto matrix_covariance_signal = EstimateCovariance(selected, sample_vars_signal, UINT_FAST32_MAX);
        auto matrix_covariance_signal_value = GetValue(matrix_covariance_signal);
        auto matrix_corr_signal = CovToCorr(matrix_covariance_signal_value);
        CreateMatrixHistosSelected(sample_vars_signal,selected,matrix_corr_signal,outfile,"correlation","Signal");
        auto matrix_covariance_bkg = EstimateCovariance(selected, sample_vars_bkg, UINT_FAST32_MAX);
        auto matrix_covariance_bkg_value = GetValue(matrix_covariance_bkg);
        auto matrix_corr_bkg = CovToCorr(matrix_covariance_bkg_value);
        CreateMatrixHistosSelected(sample_vars_bkg,selected,matrix_corr_bkg,outfile,"correlation","Background");

    }
private:
    Arguments args;
    SampleEntryCollection samples;
    std::shared_ptr<TFile> outfile;
    bool split_training_testing;
    std::mt19937 gen;
    std::uniform_int_distribution<> test_vs_training;
    MvaVariables::VarNameSet enabled_vars, disabled_vars;
    MvaVariablesStudy vars;
};
}

PROGRAM_MAIN(analysis::MvaClassification, Arguments) // definition of the main program function
//./run.sh MvaVariableCorrelation --input_path ~/Desktop/tuples --output_file ranges.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name muTau --number_events 50000
