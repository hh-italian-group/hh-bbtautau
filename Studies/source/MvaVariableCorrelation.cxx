/*! Study of correlation matrix and mutual information of BDT variables
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

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
#include "hh-bbtautau/Analysis/include/MvaVariables.h"
#include "AnalysisTools/Core/include/ProgressReporter.h"
#include "AnalysisTools/Run/include/MultiThread.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(Long64_t, number_events);
    REQ_ARG(size_t, number_variables);
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

static constexpr int Bkg_mu = -1;
static constexpr int Bke_e = -2;
static constexpr int Signal_SM = 0;

using clock = std::chrono::system_clock;

struct SampleEntry{
  std::string filename;
  double weight;
  int mass;
  SampleEntry() : weight(-1){}
};

std::ostream& operator<<(std::ostream& os, const SampleEntry& entry)
{
    os << entry.filename << " " << entry.mass << " " << entry.weight ;
    return os;
}

std::istream& operator>>(std::istream& is, SampleEntry& entry)
{
    is >> entry.filename >> entry.mass >> entry.weight ;
    return is;
}

using DataVector = std::vector<double>;
using VarData = std::map<std::string, DataVector>;
using MassVar = std::map<int, VarData>;
using VarPair = std::pair<std::string,std::string>;
using VarBand = std::map<std::string, double>;
using VarPairMI =  std::map<VarPair,double>;
using MassVarBand = std::map<int, VarBand>;
using MassVarPairMI = std::map<int, VarPairMI>;
using MassSetSelected =  std::map<int, std::set<std::string>>;
using VarPairEstCorr = std::map<VarPair,stat_estimators::EstimatedQuantity>;
using VarPairCorr =  std::map<VarPair,double>;

class MvaVariablesStudy : public MvaVariables {
private:
    std::map<std::string, double> variables;
    std::map<std::string, MassVar> all_variables;

public:
    using MvaVariables::MvaVariables;

    virtual void SetValue(const std::string& name, double value) override
    {
        variables[name] = value;
    }

    virtual void AddEventVariables(bool istraining, int mass, std::string tree, double weight) override
    {
        VarData& sample_vars = all_variables[tree][mass];
        for(const auto& name_value : variables) {
            const std::string& name = name_value.first;
            const double value = name_value.second;
            sample_vars[name].push_back(value);
        }
    }
    const MassVar& GetSampleVariables(std::string tree) const
    {
        return all_variables.at(tree);
    }
};

struct Name_ND{
    std::set<std::string> names;
    using const_iterator = std::set<std::string>::const_iterator;

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
        std::set_symmetric_difference(names.begin(),names.end(), other.begin(), other.end(), std::inserter(difference, difference.begin()));
        return !difference.size();
    }

    bool Count(const std::pair<Name_ND, double> other) const
    {
        std::set<std::string> difference;
        std::set_symmetric_difference(names.begin(),names.end(), other.first.names.begin(), other.first.names.end(), std::inserter(difference, difference.begin()));
        return !difference.size();
    }

    const_iterator begin() const { return names.begin(); }
    const_iterator end() const { return names.end(); }
    size_t size() const { return names.size(); }

    std::pair<std::string, std::string> ToPair() const
    {
        if(names.size() != 2)
            throw exception("Can't convert to pair.");
        return std::make_pair(*names.begin(), *names.rbegin());
    }
};

struct Range {
    int min, med, max;
    Range() : min(0), med(0), max(0){}
    bool Contains(const int& v) const { return v >= min && v <= max; }
    Range(Double_t _min, Double_t _med, Double_t _max ) : min(_min),med(_med), max(_max) {}
};

Range low_mass(250,270,280);
Range medium_mass(300,350,400);
Range high_energy(450,700,900);
std::vector<Range> ranges{low_mass, medium_mass, high_energy};

//Create optimal bandwidth for each variable for a single value of mass
VarBand OptimalBandwidth(const VarData& sample){
    VarBand bandwidth;
    std::map<std::string, std::future<double>> bandwidth_future;
        for (const auto& var : sample){
//            bandwidth_future[var.first] = std::async(std::launch::async,stat_estimators::OptimalBandwith<double>,var.second, 0.01);
            bandwidth_future[var.first] = run::async(stat_estimators::OptimalBandwith<double>,std::cref(var.second), 0.01);
        }
        for(auto& var : bandwidth_future) {
            var.second.wait();
            if(!var.second.valid())
                throw exception("future not valid");
            bandwidth[var.first] = var.second.get();
        }
    return bandwidth;
}
//VarBand OptimalBandwidth(const VarData& sample){
//    VarBand bandwidth;
//        for (const auto& var : sample){
//            bandwidth[var.first] = stat_estimators::OptimalBandwith(var.second, 0.01);
//        }
//    return bandwidth;
//}

//Create elements of mutual information matrix for a single value of mass


VarPairMI Mutual(const VarData& sample, const VarBand& bandwidth){
    VarPairMI mutual_matrix;
    std::map<VarPair,std::future<double>> matrix_future;
    for (const auto& var_1: sample){
        for(const auto& var_2 : sample) {
            const VarPair var_12(var_1.first, var_2.first);
            const VarPair var_21(var_2.first, var_1.first);
            if (var_2.first < var_1.first) continue;
            matrix_future[var_12] = run::async(stat_estimators::ScaledMutualInformation<double>, std::cref(var_1.second), std::cref(var_2.second), bandwidth.at(var_1.first), bandwidth.at(var_2.first));
        }
    }

    for(auto& var : matrix_future) {
        var.second.wait();
        if(!var.second.valid())
            throw exception("future not valid");
        mutual_matrix[var.first] = var.second.get();
    }

    return mutual_matrix;
}

//VarPairMI Mutual(const VarData& sample, const VarBand& bandwidth){
//    VarPairMI mutual_matrix;

//    for (const auto& var_1: sample){
//        for(const auto& var_2 : sample) {
//            const VarPair var_12(var_1.first, var_2.first);
//            const VarPair var_21(var_2.first, var_1.first);
//            if(mutual_matrix.count(var_21))
//                continue;
//            mutual_matrix[var_12] = stat_estimators::ScaledMutualInformation(var_1.second, var_2.second, bandwidth.at(var_1.first), bandwidth.at(var_2.first));
//        }
//    }

//    return mutual_matrix;
//}

//Create 2D plot of mutual information for signal and background
void MutualPlot(int mass, const VarPairMI& mutual_matrix_signal, const VarPairMI& mutual_matrix_bkg,
               TDirectory* directory){
    auto directory_s = directory->GetDirectory("1D");
    if (directory_s == nullptr) {
        directory->mkdir("1D");
        directory_s = directory->GetDirectory("1D");
        auto histo = std::make_shared<TH1D>("Background","MutualInformation_Background",50,0,1);
        histo->SetXTitle("MI");
        for (const auto& entry : mutual_matrix_bkg){
            histo->Fill(entry.second);
        }
        root_ext::WriteObject(*histo, directory_s);
    }
    auto directory_2d = directory->GetDirectory("2D");
    if (directory_2d == nullptr) {
        directory->mkdir("2D");
        directory_2d = directory->GetDirectory("2D");
    }
    std::string m;
    if ( mass ==  Signal_SM) m = "SM";
    else m = std::to_string(mass);
    auto histo2d = std::make_shared<TH2D>(("Signal_mass"+m+"_Background").c_str(),("MutualInformation_Signal_mass"+m+"_Background").c_str(),50,0,1,50,0,1);
    histo2d->SetXTitle("MI Signal");
    histo2d->SetYTitle("MI Background");
    auto histo = std::make_shared<TH1D>(("Signal_mass"+m).c_str(),("MutualInformation_Signal_mass"+m).c_str(),50,0,1);
    histo->SetXTitle("MI");
    for (const auto& entry : mutual_matrix_signal){
        histo2d->Fill(entry.second, mutual_matrix_bkg.at(entry.first));
        histo->Fill(entry.second);
    }
    root_ext::WriteObject(*histo, directory_s);
    root_ext::WriteObject(*histo2d, directory_2d);
}

//Estimate elements of covariance matrix for selected variables
VarPairEstCorr EstimateCovariance(const std::set<std::string>& selected, const VarData& sample_vars, uint_fast32_t seed){
    VarPairEstCorr covariance_matrix;

    for(const auto& var_1 : selected) {
        for(const auto& var_2 : selected) {
            const VarPair var_12(var_1, var_2);
            const VarPair var_21(var_2, var_1);
            if(covariance_matrix.count(var_21))
                continue;
            covariance_matrix[var_12] = stat_estimators::EstimateWithErrorsByResampling(stat_estimators::Covariance<double>,sample_vars.at(var_1),sample_vars.at(var_2),true,true,1000,0.31731, seed);
        }
    }

    return covariance_matrix;
}

VarPairCorr GetValue(const VarPairEstCorr & cov_matrix_signal){
    VarPairCorr map;
    for(const auto& var_1 : cov_matrix_signal) {
        auto map_cov_varpair = cov_matrix_signal.at(var_1.first);
        double value = map_cov_varpair.value;
        map[var_1.first] = value;
    }
    return map;
}

VarPairCorr CovToCorr(const VarPairCorr& cov_matrix){
    VarPairCorr correlation;
    for(const auto& elements : cov_matrix) {
        auto name_1 = elements.first;
        const VarPair var_11(name_1.first, name_1.first);
        double sigma_1 = std::sqrt(cov_matrix.at(var_11));
        const VarPair var_22(name_1.second, name_1.second);
        double sigma_2 = std::sqrt(cov_matrix.at(var_22));
        const VarPair var_12(name_1.first, name_1.second);
        correlation[var_12] = cov_matrix.at(var_12)/(sigma_1*sigma_2);
    }
    return correlation;
}

//Create covariance(correlation/mutual information) histo matrix
void CreateMatrixHistos(const MassVar& sample_vars, const MassVarPairMI& element, std::string type,
                        TDirectory* directory){
    std::string mass, class_sample;
    for(const auto mass_entry: sample_vars){
        int bin =  mass_entry.second.size();
        if (mass_entry.first == Bkg_mu && mass_entry.first != Signal_SM)
            class_sample =  "Background";
        else {
            class_sample =  "Signal";
            if (mass_entry.first == Signal_SM) mass = "_SM";
            else mass = "_mass"+std::to_string(mass_entry.first);
        }
        auto matrix = std::make_shared<TH2D>((type+"_"+class_sample+mass).c_str(),(type+"_"+class_sample+mass).c_str(),bin,0,bin,bin,0,bin);
        int i = 1;
        for(const auto& var_1 : mass_entry.second) {
            int j = 1;
            matrix->GetXaxis()->SetBinLabel(i, (var_1.first).c_str());
            matrix->GetYaxis()->SetBinLabel(i, (var_1.first).c_str());
            for(const auto& var_2 : mass_entry.second) {
                if (j >= i){
                    const VarPair var_12(var_1.first, var_2.first);
                    matrix->SetBinContent(i, j, element.at(mass_entry.first).at(var_12));
                    matrix->SetBinContent(j, i, element.at(mass_entry.first).at(var_12));
                }
                j++;
            }
            i++;
        }
        root_ext::WriteObject(*matrix, directory);
    }
}

//Create covariance(correlation/mutual information) histo matrix
void CreateMatrixHistosMassSelected(int massa, const VarPairMI& element, const std::set<std::string>& selected,
                                    std::string type, TDirectory* directory){
    std::string mass= "", class_sample = "";
    int bin =  selected.size();
    if (massa == Bkg_mu && massa != Signal_SM)
        class_sample =  "Background";
    else {
        class_sample =  "Signal";
        if (massa == Signal_SM) mass = "_SM";
        else mass = "_mass"+std::to_string(massa);
    }
    auto matrix = std::make_shared<TH2D>((type+"_"+class_sample+mass).c_str(),(type+"_"+class_sample+mass).c_str(),bin,0,bin,bin,0,bin);
    int i = 1;
    for(const auto& var_1 : selected) {
        int j = 1;
        matrix->GetXaxis()->SetBinLabel(i, (var_1).c_str());
        for(const auto& var_2 : selected) {
            VarPair var_12(var_1, var_2);
            if (!element.count(var_12)) var_12 = VarPair(var_2,var_1);
            matrix->GetYaxis()->SetBinLabel(i, (var_1).c_str());
            matrix->SetBinContent(i, j, (element.at(var_12)*100));
            j++;
        }
        i++;
    }
    root_ext::WriteObject(*matrix, directory);

}

std::shared_ptr<TGraph> CreatePlot(std::string title, std::string name, std::string x_axis, std::string y_axis){
    std::shared_ptr<TGraph> plot;

    plot = std::make_shared<TGraph>();
    plot->SetLineColor(kGreen+1);
    plot->SetLineWidth(1);
    plot->SetMarkerColor(1);
    plot->SetMarkerSize(1);
    plot->SetMarkerStyle(3);
    plot->SetTitle((title).c_str());
    plot->SetName((name).c_str());
    plot->GetHistogram()->GetXaxis()->SetTitle((x_axis).c_str());
    plot->GetHistogram()->GetYaxis()->SetTitle((y_axis).c_str());

    return plot;
}


using Name_NDDistance = std::map<Name_ND, double>;
using MassName_ND = std::map<int, Name_NDDistance>;
using VectorName_ND = std::deque<std::pair<Name_ND, double>>;
using MassVectorName_ND = std::map<int, VectorName_ND>;
using RangeSelected = std::map<int,std::set<std::string>>;

Name_NDDistance JensenDivergenceSB(const VarData& sample_signal, const VarData& sample_bkg,
                                   const VarBand& bandwidth_signal, const VarBand& bandwidth_bkg){
    Name_NDDistance  JSDivergenceSB;
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
        JSDivergenceND_future[Name_ND{entry_1.first}] = run::async(stat_estimators::JensenShannonDivergence_ND<double>, x, y, band_x, band_y);
        for (const auto& entry_2 : sample_signal){
            if ( entry_1.first == entry_2.first )
                continue;
            if (JSDivergenceND_future.count(Name_ND({entry_2.first,entry_1.first})) )
                continue;
            x.push_back(&sample_signal.at(entry_2.first));
            y.push_back(&sample_bkg.at(entry_2.first));
            band_x.push_back(bandwidth_signal.at(entry_2.first));
            band_y.push_back(bandwidth_bkg.at(entry_2.first));
            JSDivergenceND_future[Name_ND{entry_1.first,entry_2.first}] = run::async(stat_estimators::JensenShannonDivergence_ND<double>, x, y, band_x, band_y);
            x.erase(x.end() - 1);
            y.erase(y.end() - 1);
            band_x.erase(band_x.end() - 1);
            band_y.erase(band_y.end() - 1);
        }
    }
    for(auto& entry : JSDivergenceND_future) {
        JSDivergenceSB[entry.first] = entry.second.get();
    }
    return JSDivergenceSB;
}

std::set<std::string> CheckList(int mass,  MassVar sample, const Name_NDDistance& JSDivergenceND,
                                const VarPairMI& mutual_matrix_signal, const VarPairMI& mutual_matrix_bkg,
                                std::map<std::string, std::string>& eliminated, const MassVarBand& band,
                                size_t number_variables, TDirectory* directory){

    std::set<std::string> selected, not_corrected;
    VectorName_ND JSDivergence_vector(JSDivergenceND.begin(), JSDivergenceND.end());

    while(selected.size() < number_variables && JSDivergence_vector.size()) {
        std::sort(JSDivergence_vector.begin(), JSDivergence_vector.end(),
                  [](const std::pair<Name_ND, double>& el1, const std::pair<Name_ND, double>& el2){
            return el1.second > el2.second;
        });
        auto best_entry = JSDivergence_vector.front();

        for(const auto& name : best_entry.first.names) {
            if(selected.count(name))
                continue;
            const double JS_1d = JSDivergenceND.at(Name_ND({name}));
            for(auto& entry : JSDivergence_vector) {
                if(entry.first.names.count(name))
                    entry.second -= JS_1d;
            }
            if (not_corrected.count(name))
                continue;
            static constexpr double threashold = 0.5;
            for(const auto& other_entry : sample.at(mass)) {
                if(other_entry.first == name) continue;
                const Name_ND names({name, other_entry.first});
                const auto name_pair = names.ToPair();
                if(mutual_matrix_signal.at(name_pair) < threashold && mutual_matrix_bkg.at(name_pair) < threashold){
                    eliminated[other_entry.first] = "MI-"+name;
                    not_corrected.insert(other_entry.first);
                }
            }

            selected.insert(name);
        }
        VectorName_ND copy;
        for (const auto& entry : JSDivergence_vector){
            if(entry.second <= 0) continue;
            bool consider = true;
            for(const auto& other : entry.first) {
                consider = consider && !not_corrected.count(other);
            }
            if(consider)
                copy.push_back(entry);
        }
        JSDivergence_vector = std::move(copy);
    }

    auto directory_v = directory->GetDirectory("Selected_variables");
    if (directory_v == nullptr) {
        directory->mkdir("Selected_variables");
        directory_v = directory->GetDirectory("Selected_variables");
    }
    TDirectory* directory_m;
    if (mass==Signal_SM){
        directory_v->mkdir("Signal_SM");
        directory_m = directory_v->GetDirectory("Signal_SM");
    }
    else{
    directory_v->mkdir(("Signal_mass"+std::to_string(mass)).c_str());
    directory_m = directory_v->GetDirectory(("Signal_mass"+std::to_string(mass)).c_str());
    }
    for(const auto& selected_name : selected){
        std::shared_ptr<TGraph> plot;
        plot = CreatePlot(("JSD_"+selected_name).c_str(),("JSD_"+selected_name).c_str(),"mass","JSD");
        int i = 0;
        std::vector<const DataVector*> sample_1;
        DataVector band_1;
        sample_1.push_back(&sample.at(mass).at(selected_name));
        band_1.push_back(band.at(mass).at(selected_name));
        for (const auto& mass_entry: sample){
            if (mass_entry.first == Bkg_mu) continue;
            std::vector<const DataVector*> sample_2;
            DataVector band_2;
            sample_2.push_back(&sample.at(mass_entry.first).at(selected_name));
            band_2.push_back(band.at(mass_entry.first).at(selected_name));
            auto test = stat_estimators::JensenShannonDivergence_ND(sample_1, sample_2, band_1, band_2);
            plot->SetPoint(i,mass_entry.first,test);
            i++;
        }
        root_ext::WriteObject(*plot, directory_m);
    }
    std::cout<<"fine"<<std::endl;
    return selected;
}


std::set<std::string> CheckList_2(int min, int max, const std::map<int, double>& max_distance, MassVar sample,
                                  const MassName_ND& JSDivergenceND, const std::map<int, VectorName_ND>& vector,
                                  const MassVarPairMI& mutual_matrix, std::map<std::string, std::string>& eliminated,
                                  const MassVarBand& band, size_t number_variables, TDirectory* directory){

    int count = 0;
    std::set<std::string> selected, checked, not_corrected;
    std::map<int, VectorName_ND> JSDivergence_vector = vector;
    std::map<std::pair<int,int>, std::shared_ptr<TH1D>> histo_check;
    std::map<std::pair<int, int>, std::map<Name_ND,double>> JSDivergence_signal;
    std::ofstream ofc("Check2_JSD_bestentry"+std::to_string(min)+"-"+std::to_string(max)+".csv", std::ofstream::out);
    for (const auto& mass_1 : JSDivergenceND){
        if (mass_1.first < min || mass_1.first > max || mass_1.first ==  Signal_SM) continue;
        for (const auto& mass_2 : JSDivergenceND){
            if (mass_2.first < min || mass_2.first > max || mass_2.first ==  Signal_SM) continue;
            if (mass_2.first < mass_1.first) continue;
            std::pair<int,int> pair(mass_1.first, mass_2.first);
            histo_check[pair] = std::make_shared<TH1D>(("JSDrangecheck_"+std::to_string(mass_1.first)+"_"+std::to_string(mass_2.first)).c_str(), ("JensenShannonDivergence_"+std::to_string(mass_1.first)+"_"+std::to_string(mass_2.first)).c_str(), 30,0,1);
            histo_check[pair]->SetXTitle("JSD");
        }
    }

    auto directory_s = directory->GetDirectory("Selected");
    if (directory_s == nullptr) {
        directory->mkdir("Selected");
        directory_s = directory->GetDirectory("Selected");
    }
    auto directory_bs = directory->GetDirectory("Selected_Bkg");
    if (directory_bs == nullptr) {
        directory->mkdir("Selected_Bkg");
        directory_bs = directory->GetDirectory("Selected_Bkg");
    }
    while(selected.size() < number_variables && JSDivergence_vector.at(min).size()>1) {
        std::sort(JSDivergence_vector.at(min).begin(), JSDivergence_vector.at(min).end(), [](const std::pair<Name_ND, double>& el1, const std::pair<Name_ND, double>& el2){
            return el1.second > el2.second;
        });
        while(JSDivergence_vector.at(min).size() && JSDivergence_vector.at(min).back().second < 0){
            JSDivergence_vector.at(min).pop_back();
        }
        if(!JSDivergence_vector.at(min).size()) continue;
        VectorName_ND distance;
        for (const auto& entry: JSDivergence_vector.at(min)){
            double d = 0;
            for (const auto& mass: JSDivergenceND){
                if (mass.first ==  Signal_SM) continue;
                if (mass.first < min || mass.first > max) continue;
                d+= (max_distance.at(mass.first)-mass.second.at(entry.first)) / max_distance.at(mass.first);
            }
            distance.emplace_back(entry.first,d);
        }
        std::sort(distance.begin(), distance.end(), [](const std::pair<Name_ND, double>& el1, const std::pair<Name_ND, double>& el2){
        return el1.second < el2.second;
        });
        const auto& best_entry = distance.front();
        bool passed_be = true;

        for(const auto& name : best_entry.first.names) {
            eliminated[name] = "";
            checked.insert(name);
            if(selected.count(name)){

                continue;
            }
            for (const auto& mass : JSDivergenceND){
                if (mass.first ==  Signal_SM) continue;
                const double JS_1d = mass.second.at(Name_ND({name}));
                for(auto& entry : JSDivergence_vector.at(mass.first)) {
                    if(entry.first.names.count(name))
                        entry.second -= JS_1d;
                }
            }

            if (not_corrected.count(name)) continue;
            static const double threashold = 0.5;


            for (const auto& mass : JSDivergenceND){
                if (mass.first ==  Signal_SM) continue;
                for(const auto& other_entry : sample.at(mass.first)) {
                    if(other_entry.first == name) continue;
                    const Name_ND names({name, other_entry.first});
                    const auto name_pair = names.ToPair();
                    if(mutual_matrix.at(mass.first).at(name_pair) < threashold && mutual_matrix.at(Bkg_mu).at(name_pair) < threashold){
                        eliminated[other_entry.first] = "MI-"+name;
                        not_corrected.insert(other_entry.first);
                    }
                }
            }
            selected.insert(name);
        }
        std::cout<<"selezionati "<<selected.size()<<std::endl;

        VectorName_ND copy;
        for (const auto& entry : JSDivergence_vector.at(min)){
            if(entry.second <= 0) continue;
            bool consider = true;
            for(const auto& other : entry.first) {
                consider = consider && !not_corrected.count(other);
            }
            if(consider)
                copy.push_back(entry);
        }
        JSDivergence_vector.at(min) = std::move(copy);


        for(auto mass : JSDivergence_vector){
            int c = 0;
            for (auto entry : mass.second){
                if (entry.first.Count(best_entry)){
                    JSDivergence_vector.at(mass.first).erase(JSDivergence_vector.at(mass.first).begin() + c);
                }
                c++;
            }
        }


        if (passed_be == true) {
            count++;
            auto selected_name1 = *best_entry.first.names.begin();
            auto selected_name2 = *best_entry.first.names.rbegin();
            for (const auto& mass_1 : JSDivergenceND){
                if (mass_1.first < min || mass_1.first > max) continue;
                for (const auto& mass_2 : JSDivergenceND){
                    if (mass_2.first < min || mass_2.first > max) continue;
                    std::pair<int,int> pair(mass_1.first,mass_2.first);
                    if (mass_2.first < mass_1.first) continue;
                    std::vector<const DataVector*> sample_1, sample_2;
                    DataVector band_1, band_2;
                    sample_1.push_back(&sample.at(mass_1.first).at(selected_name1));
                    sample_2.push_back(&sample.at(mass_2.first).at(selected_name1));
                    band_1.push_back(band.at(mass_1.first).at(selected_name1));
                    band_2.push_back(band.at(mass_2.first).at(selected_name1));
                    if (selected_name1 == selected_name2){
                        std::cout<<"uguali"<<std::endl;
                        JSDivergence_signal[pair][Name_ND{selected_name1}] = stat_estimators::JensenShannonDivergence_ND(sample_1, sample_2, band_1, band_2);
                        histo_check [pair]->Fill(JSDivergence_signal.at(pair).at(Name_ND{selected_name1}));
                        ofc<<mass_1.first<<","<<mass_2.first<<","<<selected_name1<<","<<","<<JSDivergence_signal.at(pair).at(Name_ND{selected_name1});
                    }
                    else{
                        sample_1.push_back(&sample.at(mass_1.first).at(selected_name2));
                        sample_2.push_back(&sample.at(mass_2.first).at(selected_name2));
                        band_1.push_back(band.at(mass_1.first).at(selected_name2));
                        band_2.push_back(band.at(mass_2.first).at(selected_name2));
                        JSDivergence_signal[pair][Name_ND{selected_name1,selected_name2}] = stat_estimators::JensenShannonDivergence_ND(sample_1, sample_2, band_1, band_2);
                        histo_check[pair]->Fill(JSDivergence_signal.at(pair).at(Name_ND{selected_name1,selected_name2}));
                        ofc<<mass_1.first<<","<<mass_2.first<<","<<selected_name1<<","<<selected_name2<<","<<JSDivergence_signal.at(pair).at(Name_ND{selected_name1,selected_name2})<<std::endl;
                    }
                }
            }

            std::shared_ptr<TGraph> plot;
            plot = std::make_shared<TGraph>();
            plot->SetLineColor(kGreen+1);
            plot->SetLineWidth(1);
            plot->SetMarkerColor(1);
            plot->SetMarkerSize(1);
            plot->SetMarkerStyle(3);
            plot->SetTitle(("JSD_"+selected_name1+"_"+selected_name2+"_Range"+std::to_string(min)+"_"+std::to_string(max)).c_str());
            plot->SetName(("JSD_"+selected_name1+"_"+selected_name2+"_Range"+std::to_string(min)+"_"+std::to_string(max)).c_str());
            plot->GetHistogram()->GetXaxis()->SetTitle("mass");
            plot->GetHistogram()->GetYaxis()->SetTitle("JSD");
            int i = 0;
            std::vector<const DataVector*> sample_1;
            DataVector band_1;
            sample_1.push_back(&sample.at(min).at(selected_name1));
            sample_1.push_back(&sample.at(min).at(selected_name2));
            band_1.push_back(band.at(min).at(selected_name1));
            band_1.push_back(band.at(min).at(selected_name2));
            for (const auto& mass_entry: sample){
                if (mass_entry.first < min || mass_entry.first > max) continue;
                std::vector<const DataVector*> sample_2;
                DataVector band_2;
                sample_2.push_back(&sample.at(mass_entry.first).at(selected_name1));
                sample_2.push_back(&sample.at(mass_entry.first).at(selected_name2));
                band_2.push_back(band.at(mass_entry.first).at(selected_name1));
                band_2.push_back(band.at(mass_entry.first).at(selected_name2));
                auto test = stat_estimators::JensenShannonDivergence_ND(sample_1, sample_2, band_1, band_2);
                plot->SetPoint(i,mass_entry.first,test);
                i++;
            }
            root_ext::WriteObject(*plot, directory_s);


            std::shared_ptr<TGraph> plot_b;
            plot_b = std::make_shared<TGraph>();
            plot_b->SetLineColor(kGreen+1);
            plot_b->SetLineWidth(1);
            plot_b->SetMarkerColor(1);
            plot_b->SetMarkerSize(1);
            plot_b->SetMarkerStyle(3);
            plot_b->SetTitle(("JSD_"+selected_name1+"_"+selected_name2+"_Range"+std::to_string(min)+"_"+std::to_string(max)+"_SignalBKg").c_str());
            plot_b->SetName(("JSD_"+selected_name1+"_"+selected_name2+"_Range"+std::to_string(min)+"_"+std::to_string(max)+"_SignalBKg").c_str());
            plot_b->GetHistogram()->GetXaxis()->SetTitle("mass");
            plot_b->GetHistogram()->GetYaxis()->SetTitle("JSD");
            int i_b = 0;
            std::vector<const DataVector*> sample_1_b;
            DataVector band_1_b;
            sample_1_b.push_back(&sample.at(Bkg_mu).at(selected_name1));
            sample_1_b.push_back(&sample.at(Bkg_mu).at(selected_name2));
            band_1_b.push_back(band.at(Bkg_mu).at(selected_name1));
            band_1_b.push_back(band.at(Bkg_mu).at(selected_name2));
            for (const auto& mass_entry: sample){
                if (mass_entry.first < min || mass_entry.first > max) continue;
                std::vector<const DataVector*> sample_2;
                DataVector band_2;
                sample_2.push_back(&sample.at(mass_entry.first).at(selected_name1));
                sample_2.push_back(&sample.at(mass_entry.first).at(selected_name2));
                band_2.push_back(band.at(mass_entry.first).at(selected_name1));
                band_2.push_back(band.at(mass_entry.first).at(selected_name2));
                auto test = stat_estimators::JensenShannonDivergence_ND(sample_1_b, sample_2, band_1_b, band_2);
                plot_b->SetPoint(i_b,mass_entry.first,test);
                i_b++;
            }
            root_ext::WriteObject(*plot_b, directory_bs);
        }
    }

    std::map<std::pair<int,int>, std::shared_ptr<TH1D>> histo;
    for (const auto& pair_mass : JSDivergence_signal){
        histo[pair_mass.first] = std::make_shared<TH1D>(("JSDrange_"+std::to_string(pair_mass.first.first)+"_"+std::to_string(pair_mass.first.second)).c_str(), ("Jensen Shannon Divergence "+std::to_string(pair_mass.first.first)+"_"+std::to_string(pair_mass.first.second)).c_str(), 50,0,1);
        histo[pair_mass.first]->SetXTitle("JSD");
    }
    for (const auto& pair_mass : JSDivergence_signal){
        std::map<Name_ND, bool> inserted;
        for(const auto& selected_name1 : selected) {
            std::vector<const DataVector*> sample_1, sample_2;
            DataVector band_1, band_2;
            sample_1.push_back(&sample.at(pair_mass.first.first).at(selected_name1));
            sample_2.push_back(&sample.at(pair_mass.first.second).at(selected_name1));
            band_1.push_back(band.at(pair_mass.first.first).at(selected_name1));
            band_2.push_back(band.at(pair_mass.first.second).at(selected_name1));
            if(!inserted.count(Name_ND{selected_name1})){
                if (JSDivergence_signal.at(pair_mass.first).count(Name_ND{selected_name1})){
                    histo[pair_mass.first]->Fill(JSDivergence_signal.at(pair_mass.first).at(Name_ND{selected_name1}));
                    inserted[Name_ND{selected_name1}] = true;
                }
                else{
                    JSDivergence_signal[pair_mass.first][Name_ND{selected_name1}] = stat_estimators::JensenShannonDivergence_ND(sample_1, sample_2, band_1, band_2);
                    histo[pair_mass.first]->Fill(JSDivergence_signal.at(pair_mass.first).at(Name_ND{selected_name1}));
                    inserted[Name_ND{selected_name1}] = true;
                }
            }
            for(const auto& selected_name2 : selected) {
                if (selected_name1 == selected_name2) continue;
                if (inserted.count(Name_ND({selected_name2, selected_name1}))) continue;
                if (!inserted.count(Name_ND({selected_name1, selected_name2}))){
                    if (JSDivergence_signal.at(pair_mass.first).count(Name_ND{selected_name1,selected_name2})){
                        histo[pair_mass.first]->Fill(JSDivergence_signal.at(pair_mass.first).at(Name_ND{selected_name1,selected_name2}));
                        inserted[Name_ND{selected_name1,selected_name2}] = true;
                    }
                    else{
                        sample_1.push_back(&sample.at(pair_mass.first.first).at(selected_name2));
                        sample_2.push_back(&sample.at(pair_mass.first.second).at(selected_name2));
                        band_1.push_back(band.at(pair_mass.first.first).at(selected_name2));
                        band_2.push_back(band.at(pair_mass.first.second).at(selected_name2));
                        JSDivergence_signal[pair_mass.first][Name_ND{selected_name1,selected_name2}] = stat_estimators::JensenShannonDivergence_ND(sample_1, sample_2, band_1, band_2);
                        histo[pair_mass.first]->Fill(JSDivergence_signal.at(pair_mass.first).at(Name_ND{selected_name1,selected_name2}));
                        inserted[Name_ND{selected_name1,selected_name2}] = true;
                        sample_1.erase(sample_1.end() - 1);
                        sample_2.erase(sample_2.end() - 1);
                        band_1.erase(band_1.end() - 1);
                        band_2.erase(band_2.end() - 1);
                    }
                }
                else
                    continue;
            }
        }
    }
    std::cout<<count<<std::endl;
    auto directory_ss = directory->GetDirectory("Selected_pair_masses");
    if (directory_ss == nullptr) {
        directory->mkdir("Selected_pair_masses");
        directory_ss = directory->GetDirectory("Selected_pair_masses");
    }
    auto directory_ps = directory->GetDirectory("Pair_selected");
    if (directory_ps == nullptr) {
        directory->mkdir("Pair_selected");
        directory_ps = directory->GetDirectory("Pair_selected");
    }
    for (const auto& pair_mass : JSDivergence_signal){
        root_ext::WriteObject(*histo[pair_mass.first], directory_ps);
        root_ext::WriteObject(*histo_check[pair_mass.first], directory_ss);
    }

    auto directory_v = directory->GetDirectory("Variables");
    if (directory_v == nullptr) {
        directory->mkdir("Variables");
        directory_v = directory->GetDirectory("Variables");
    }
    for(const auto& selected_name : selected){
        std::shared_ptr<TGraph> plot = CreatePlot("JSD_"+selected_name+"_Range"+std::to_string(min)+"_"+std::to_string(max), "JSD_"+selected_name+"_Range"+std::to_string(min)+"_"+std::to_string(max) ,"mass", "JSD");
        int i = 0;
        std::vector<const DataVector*> sample_1;
        DataVector band_1;
        sample_1.push_back(&sample.at(min).at(selected_name));
        band_1.push_back(band.at(min).at(selected_name));
        for (const auto& mass_entry: sample){
            if (mass_entry.first < min || mass_entry.first > max) continue;
            std::vector<const DataVector*> sample_2;
            DataVector band_2;
            sample_2.push_back(&sample.at(mass_entry.first).at(selected_name));
            band_2.push_back(band.at(mass_entry.first).at(selected_name));
            auto test = stat_estimators::JensenShannonDivergence_ND(sample_1, sample_2, band_1, band_2);
            plot->SetPoint(i,mass_entry.first,test);
            i++;
        }
        root_ext::WriteObject(*plot, directory_v);
    }
   return selected;
}

void Plot(int massa, const std::set<std::string>& selected,const VarData& sample, TDirectory* directory){
    std::set<std::string> checked;
    std::string mass;
    if (massa == Bkg_mu) mass = "Bkg";
    else if (massa == Signal_SM) mass = "SM";
    else mass = std::to_string(massa);
    auto dir = directory->GetDirectory(mass.c_str());
    if (dir == nullptr){
        directory->mkdir(mass.c_str());
        dir = directory->GetDirectory(mass.c_str());
    }
    for (const auto & set_entry_1 :  selected){
        checked.insert(set_entry_1);
        auto histo = std::make_shared<TH1D>((set_entry_1+"_"+mass).c_str(), (set_entry_1+"_"+mass).c_str(), 50,0,0);
        histo->SetCanExtend(TH1::kAllAxes);
        histo->SetXTitle((set_entry_1).c_str());
        for(size_t i = 0; i < sample.at(set_entry_1).size(); i++){
                histo->Fill(sample.at(set_entry_1)[i]);
        }
        root_ext::WriteObject(*histo, dir);
    }
}


void Plot2D(int massa, const std::set<std::string>& selected,const VarData& sample, TDirectory* directory){
    std::set<std::string> checked;
    std::string mass;
    if (massa == Bkg_mu) mass = "Bkg";
    else if (massa == Signal_SM) mass = "SM";
    else mass = std::to_string(massa);
    auto dir = directory->GetDirectory(mass.c_str());
    if (dir == nullptr){
        directory->mkdir(mass.c_str());
        dir = directory->GetDirectory(mass.c_str());
    }
    for (const auto & set_entry_1 :  selected){
        checked.insert(set_entry_1);
        for (const auto & set_entry_2 :  selected){
            if (checked.count(set_entry_2)) continue;
            auto histo = std::make_shared<TH2D>((set_entry_1+"_"+set_entry_2+"_"+mass).c_str(), (set_entry_1+"_"+set_entry_2+"_"+mass).c_str(), 50,0,0,50,0,0);
            histo->SetCanExtend(TH1::kAllAxes);
            histo->SetXTitle((set_entry_1).c_str());
            histo->SetYTitle((set_entry_2).c_str());
            for(size_t i = 0; i < sample.at(set_entry_1).size(); i++){
                    histo->Fill(sample.at(set_entry_1)[i],sample.at(set_entry_2)[i]);
            }
            root_ext::WriteObject(*histo, dir);
        }
    }
}

void PlotJensenShannon(const MassName_ND& JSDivergenceND, TDirectory* directory){
    std::map<std::string, std::shared_ptr<TGraph>> plot;
    int i = 0;
    for (const auto& mass_entry: JSDivergenceND){
        if (mass_entry.first == Bkg_mu)
            continue;
        for (const auto& var : mass_entry.second){
            if (var.first.names.size()!=1)
                continue;
            const std::string name = *var.first.names.begin();
            if (plot.count(name) == 0) {
                plot[name] = CreatePlot("JSD_"+name+"_SignalBkg","JSD_"+name+"_SignalBkg","mass","JSD");
            }
            plot[name]->SetPoint(i,mass_entry.first,JSDivergenceND.at(mass_entry.first).at(Name_ND({name})));
        }
        i++;
    }
    for(const auto& var: plot){
        root_ext::WriteObject(*var.second, directory);
    }
}








using VecVariables =  std::vector<std::string>;
std::vector<VecVariables> VariablesSelection(const MassVar& sample, const MassVarBand& band,
                        const MassVarPairMI& mutual_matrix, MassSetSelected& selected, RangeSelected& selected2,
                        std::shared_ptr<TFile> outfile, size_t number_variables){

    MassName_ND JSDivergenceSB;
    std::cout<<"Jensen Shannon Signal Background";
    auto start = clock::now();
    for (const auto& mass_entry: sample){
        if (mass_entry.first == Bkg_mu)
            continue;
        JSDivergenceSB[mass_entry.first] = JensenDivergenceSB(mass_entry.second, sample.at(Bkg_mu), band.at(mass_entry.first), band.at(Bkg_mu));
    }
    std::cout<<" - Plot Jensen Shannon"<<std::endl;
    outfile->mkdir("JensenShannonDivergence");
    auto directory_js = outfile->GetDirectory("JensenShannonDivergence");
    directory_js->mkdir("Signal_Background");
    auto directory_sb = directory_js->GetDirectory("Signal_Background");
    directory_sb->mkdir("All_variables");
    auto directory_v = directory_sb->GetDirectory("All_variables");
    PlotJensenShannon(JSDivergenceSB, directory_v);
    auto stop = clock::now();
    std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;

    std::cout<<"Checklist1"<<std::endl;
    std::ofstream ofs("SelectedVariable.csv", std::ofstream::out);
    directory_sb->mkdir("1D");
    auto directory_1d = directory_sb->GetDirectory("1D");
    directory_js->mkdir("Signal_Signal");
    auto directory_ss = directory_js->GetDirectory("Signal_Signal");
    directory_ss->mkdir("Check_Mass");
    auto directory_cm = directory_ss->GetDirectory("Check_Mass");
    for (const auto& mass_entry : sample){
        std::cout<<mass_entry.first<<" ";
        std::cout.flush();
        if (mass_entry.first == Bkg_mu)
            continue;
        std::map<std::string, std::string> eliminated;
        selected[mass_entry.first] =  CheckList(mass_entry.first, sample,JSDivergenceSB.at(mass_entry.first), mutual_matrix.at(mass_entry.first), mutual_matrix.at(Bkg_mu), eliminated, band, number_variables, directory_cm);
        ofs << mass_entry.first << std::endl;
        for(auto& entry_selected: selected[mass_entry.first]){
                ofs << entry_selected<<",";
        }
        ofs << std::endl;
        std::string  mass = std::to_string(mass_entry.first);
        if (mass_entry.first == Signal_SM) mass = "SM";
        auto histo = std::make_shared<TH1D>(("JSD_Signal"+mass+"_Background").c_str(), ("JensenShannonDivergence_Signal"+mass+"_Background").c_str(), 50,0,1);
        histo->SetXTitle("JSD");
        std::ofstream of("InformationTable_mass"+(std::to_string(mass_entry.first))+".csv", std::ofstream::out);
        of<<"Var_1"<<","<<"Var_2"<<","<<"JSD_ND"<<","<<"JSD_12-(JSD_1+JSD_2)"<<","<<"ScaledMI_Signal_12" <<","<<"ScaledMI_Bkg_12"<<","<<"selected"<<","<<"eliminated by"<<std::endl;
        for(auto& entry : JSDivergenceSB.at(mass_entry.first)){
             histo->Fill(entry.second);
             if ((entry.first).names.size() == 1){
                 auto name_1 = ((entry.first).names).begin();
                 bool selected_1 = selected[mass_entry.first].count(*name_1);
                 of<<*name_1<<","<<"   "<<","<<JSDivergenceSB.at(mass_entry.first).at(Name_ND({*name_1}))<<","<<","<<","<<","<<selected_1;
                 if(eliminated.count(*name_1)) of<<","<<eliminated.at(*name_1);
                 of<<std::endl;
             }
             else  if ((entry.first).names.size() == 2){
                 auto name_1 = ((entry.first).names).begin();
                 auto name_2 = ((entry.first).names).rbegin();
                 bool selected_1 = selected[mass_entry.first].count(*name_1);
                 bool selected_2 = selected[mass_entry.first].count(*name_2);
                 of<<*name_1<<","<<*name_2<<","<<JSDivergenceSB.at(mass_entry.first).at(Name_ND({*name_1,*name_2}))<<","<<JSDivergenceSB.at(mass_entry.first).at(Name_ND({*name_1,*name_2}))-(JSDivergenceSB.at(mass_entry.first).at(Name_ND({*name_1}))+JSDivergenceSB.at(mass_entry.first).at(Name_ND({*name_2})))<<","<<mutual_matrix.at(mass_entry.first).at(VarPair(*name_1,*name_2))<<","<<mutual_matrix.at(Bkg_mu).at(VarPair(*name_1,*name_2))<<","<<selected_1<<","<<selected_2;
                 if(eliminated.count(*name_1)) of<<","<<eliminated.at(*name_1);
                 if(eliminated.count(*name_2)) of<<","<<eliminated.at(*name_2);
                 of<<std::endl;
             }
        }
        of.close();
        root_ext::WriteObject(*histo, directory_1d);
    }
    ofs.close();

    std::cout<<"intersection and union"<<std::endl;
    int bin = selected.size();
    auto matrix_intersection = std::make_shared<TH2D>("Intersection","Intersection of variables",bin, 0, bin, bin, 0, bin);
    matrix_intersection->SetXTitle("mass");
    matrix_intersection->SetYTitle("mass");
    int k = 1;
    std::vector<VecVariables> vec_union(3);
//    for (int i = 0; i <= 3; i++){
//        VecVariables unione;
//        vec_union.emplace_back(unione);
//    }
    for(const auto& mass_1: selected){
        int j = 1;
        std::string mass;
        if (mass_1.first != Bkg_mu && mass_1.first != Signal_SM) mass = std::to_string(mass_1.first);
        else if (mass_1.first == Signal_SM) mass = "SM";
        matrix_intersection->GetXaxis()->SetBinLabel(k, (mass).c_str());
        matrix_intersection->GetYaxis()->SetBinLabel(k, (mass).c_str());
        for(const auto& mass_2: selected){
            VecVariables intersection;
            std::set_intersection(mass_1.second.begin(), mass_1.second.end(), mass_2.second.begin(), mass_2.second.end(),
                                    std::back_inserter(intersection));
            matrix_intersection->SetBinContent(k, j, intersection.size());
            int i = 0;
            for (const auto& range: ranges){
                if (mass_1.first <= range.max && mass_1.first >= range.min && mass_2.first <= range.max && mass_1.first >= range.min){
                     std::set_union(mass_1.second.begin(), mass_1.second.end(), mass_2.second.begin(), mass_2.second.end(),
                                                        std::back_inserter(vec_union[i]));
                }
                i++;
            }
            j++;
        }
        k++;
    }
    root_ext::WriteObject(*matrix_intersection, outfile.get());
    outfile->mkdir("VariableDistribution1D");
    auto directory_vardistr1d = outfile->GetDirectory("VariableDistribution1D");
    directory_vardistr1d->mkdir("Mass_selected");
    auto directory_vardistrms = directory_vardistr1d->GetDirectory("Mass_selected");
    outfile->mkdir("VariableDistribution2D");
    auto directory_vardistr2d = outfile->GetDirectory("VariableDistribution2D");
    directory_vardistr2d->mkdir("Mass_selected");
    auto directory_vardistrms2d = directory_vardistr2d->GetDirectory("Mass_selected");
    std::cout<<"plot"<<std::endl;
    for (const auto& mass_entry : sample){
        if (mass_entry.first == Bkg_mu) {
            Plot(mass_entry.first, selected.at(250), sample.at(mass_entry.first), directory_vardistrms);
            Plot2D(mass_entry.first, selected.at(250), sample.at(mass_entry.first), directory_vardistrms2d);

        }
        else{
            Plot(mass_entry.first, selected.at(mass_entry.first), sample.at(mass_entry.first), directory_vardistrms);
            Plot2D(mass_entry.first, selected.at(mass_entry.first), sample.at(mass_entry.first), directory_vardistrms2d);
        }
    }

    for (auto& names : vec_union){
        std::sort(names.begin(),names.end() );
        names.erase( std::unique( names.begin(), names.end() ), names.end() );
    }
    std::cout<<"#variabli per range dopo unione: "<<vec_union[0].size()<<" "<<vec_union[1].size()<<" "<<vec_union[2].size()<<std::endl;

    directory_vardistr1d->mkdir("Range_united");
    auto directory_vardistrru = directory_vardistr1d->GetDirectory("Range_united");
    directory_vardistr2d->mkdir("Range_united");
    auto directory_vardistrru2d = directory_vardistr2d->GetDirectory("Range_united");
    std::cout<<"plot"<<std::endl;
//    int i = 0;
//    for (const auto& range: ranges){
//        std::set<std::string> sel;
//        for (const auto& name : vec_union[i]){
//            sel.insert(name);
//        }
//        for (const auto& mass_entry : sample){
//            if (mass_entry.first == Bkg_mu) {
//                Plot(mass_entry.first, sel, sample.at(mass_entry.first), directory_vardistrru);
//                Plot2D(mass_entry.first, sel, sample.at(mass_entry.first), directory_vardistrru2d);

//            }
//            if (mass_entry.first > range.max || mass_entry.first < range.min ) continue;
//            else {
//                Plot(mass_entry.first, sel, sample.at(mass_entry.first), directory_vardistrru);
//                Plot2D(mass_entry.first, sel, sample.at(mass_entry.first), directory_vardistrru2d);
//            }
//        }
//        i++;
//    }
    std::vector<VecVariables> vector_union;
    std::ofstream of("UnionVariables.csv", std::ofstream::out);
    int j = 0;
    for (const auto& range: ranges){
        vector_union.emplace_back(vec_union[j]);
        of<<"Massa"<<range.min<<"-"<<range.max<<std::endl;
        for(const auto& entry: vec_union[j]){
            of<<entry<<",";
        }
        of<<std::endl;
        j++;
    }
    of.close();

    std::cout<<"Checklist2"<<std::endl;
    std::map<std::string, std::string> eliminated2;
    std::map<int, VectorName_ND> JSDivergence_vector;
    std::map<int, double> max_distance;
    for (const auto& range: ranges){
        for (const auto& mass: JSDivergenceSB){
            if (mass.first < range.min || mass.first > range.max) continue;
            VectorName_ND  vector(mass.second.begin(), mass.second.end());
            JSDivergence_vector[mass.first] = vector;
            std::sort(JSDivergence_vector.at(mass.first).begin(), JSDivergence_vector.at(mass.first).end(), [](const std::pair<Name_ND, double>& el1, const std::pair<Name_ND, double>& el2){
                return el1.second > el2.second;
            });
            max_distance[mass.first] = JSDivergence_vector.at(mass.first).front().second;
        }
    }
    std::ofstream ofr("Selected_range.csv", std::ofstream::out);
    directory_ss->mkdir("Check_Range");
    auto directory_cr = directory_ss->GetDirectory("Check_Range");
    directory_vardistr1d->mkdir("Range");
    auto directory_vardistrr = directory_vardistr1d->GetDirectory("Range");
    directory_vardistr2d->mkdir("Range");
    auto directory_vardistrr2d = directory_vardistr2d->GetDirectory("Range");
    for (const auto& range: ranges){
        selected2[range.min] = CheckList_2(range.min, range.max, max_distance, sample, JSDivergenceSB, JSDivergence_vector, mutual_matrix, eliminated2, band, 20, directory_cr);
        std::cout<<"range: "<<range.min<<"-"<<range.max<<std::endl;
        ofr<<"range: "<<range.min<<"-"<<range.max<<std::endl;
        for (const auto& name : selected2[range.min]){
            ofr<<name<<",";
        }
        ofr<<std::endl;
        for (const auto& mass_entry : sample){
            if (mass_entry.first == Bkg_mu) {
                Plot(mass_entry.first, selected2[range.min], sample.at(mass_entry.first), directory_vardistrr);
                Plot2D(mass_entry.first, selected2[range.min], sample.at(mass_entry.first), directory_vardistrr2d);
            }
            if (mass_entry.first > range.max || mass_entry.first < range.min ) continue;
            else {
                Plot(mass_entry.first, selected2[range.min], sample.at(mass_entry.first), directory_vardistrr);
                Plot2D(mass_entry.first, selected2[range.min], sample.at(mass_entry.first), directory_vardistrr2d);
            }
        }
    }
    ofr.close();

    directory_ss->mkdir("Range_Bkg");
    auto directory_rb = directory_ss->GetDirectory("Range_Bkg");
    directory_rb->mkdir("Selected");
    auto directory_s = directory_rb->GetDirectory("Selected");

    for (int count_i = 0; count_i<3; count_i++){
        for (const auto& selected_name : selected2[ranges[count_i].min]){
            std::shared_ptr<TGraph> plot;
            plot = std::make_shared<TGraph>();
            plot->SetLineColor(kGreen+1);
            plot->SetLineWidth(1);
            plot->SetMarkerColor(1);
            plot->SetMarkerSize(1);
            plot->SetMarkerStyle(3);
            plot->SetTitle(("JSD_"+selected_name+"_SignalBKg").c_str());
            plot->SetName(("JSD_"+selected_name+"_SignalBKg").c_str());
            plot->GetHistogram()->GetXaxis()->SetTitle("mass");
            plot->GetHistogram()->GetYaxis()->SetTitle("JSD");
            int j = 0;
            for (const auto& range: ranges){
                std::vector<const DataVector*> sample_1, sample_2;
                DataVector band_1, band_2, entries;
                for (const auto& mass_entry: sample){
                    if (mass_entry.first < range.min || mass_entry.first > range.max) continue;
                    for (const auto& entry : sample.at(mass_entry.first).at(selected_name)){
                        entries.push_back(entry);
                    }
                }
                VarData vardata{{selected_name, entries}};
                VarBand band_s = OptimalBandwidth(vardata);

                sample_1.push_back(&entries);
                band_1.push_back(band_s.at(selected_name));
                sample_2.push_back(&sample.at(Bkg_mu).at(selected_name));
                band_2.push_back(band.at(Bkg_mu).at(selected_name));

                auto test = stat_estimators::JensenShannonDivergence_ND(sample_1, sample_2, band_1, band_2);
                plot->SetPoint(j,range.min,test);
                j++;
            }
            root_ext::WriteObject(*plot, directory_s);
        }
    }

    directory_rb->mkdir("Pair_selected");
    auto directory_ps = directory_rb->GetDirectory("Pair_selected");
    std::set<std::string> inserted;

    for (int count_i = 0; count_i<3; count_i++){
        for (const auto& selected_name1 : selected2[ranges[count_i].min]){
            inserted.insert(selected_name1);
            for (const auto& selected_name2 : selected2[ranges[count_i].min]){
                if (inserted.count(selected_name2)) continue;
                std::shared_ptr<TGraph> plot;
                plot = std::make_shared<TGraph>();
                plot->SetLineColor(kGreen+1);
                plot->SetLineWidth(1);
                plot->SetMarkerColor(1);
                plot->SetMarkerSize(1);
                plot->SetMarkerStyle(3);
                plot->SetTitle(("JSD_"+selected_name1+"_"+selected_name2+"_SignalBKg").c_str());
                plot->SetName(("JSD_"+selected_name1+"_"+selected_name2+"_SignalBKg").c_str());
                plot->GetHistogram()->GetXaxis()->SetTitle("mass");
                plot->GetHistogram()->GetYaxis()->SetTitle("JSD");
                int j = 0;
                for (const auto& range: ranges){
                    std::vector<const DataVector*> sample_1, sample_2;
                    DataVector band_1, band_2, entries1, entries2;
                    for (const auto& mass_entry: sample){
                        if (mass_entry.first < range.min || mass_entry.first > range.max) continue;
                        for (const auto& entry : sample.at(mass_entry.first).at(selected_name1)){
                            entries1.push_back(entry);
                        }
                        for (const auto& entry : sample.at(mass_entry.first).at(selected_name2)){
                            entries2.push_back(entry);
                        }
                    }
                    VarData vardata1{{selected_name1, entries1}};
                    VarBand band_s1 = OptimalBandwidth(vardata1);
                    VarData vardata2{{selected_name2, entries2}};
                    VarBand band_s2 = OptimalBandwidth(vardata2);

                    sample_1.push_back(&entries1);
                    sample_1.push_back(&entries2);
                    band_1.push_back(band_s1.at(selected_name1));
                    band_1.push_back(band_s2.at(selected_name2));
                    sample_2.push_back(&sample.at(Bkg_mu).at(selected_name1));
                    band_2.push_back(band.at(Bkg_mu).at(selected_name1));
                    sample_2.push_back(&sample.at(Bkg_mu).at(selected_name2));
                    band_2.push_back(band.at(Bkg_mu).at(selected_name2));
                    auto test = stat_estimators::JensenShannonDivergence_ND(sample_1, sample_2, band_1, band_2);
                    plot->SetPoint(j,range.min,test);
                    j++;
                }
                root_ext::WriteObject(*plot, directory_ps);
            }
        }
    }



    auto matrix_intersection2 = std::make_shared<TH2D>("Intersection_ranges","Intersection of variables",ranges.size(), 0, ranges.size(), ranges.size(), 0, ranges.size());
    matrix_intersection2->SetXTitle("range");
    matrix_intersection2->SetYTitle("range");
    int c = 1;
    for (const auto& range1: ranges){
        for (const auto& range2: ranges){
            if (range1.min != range2.min) continue;
            std::string label;
            label = std::to_string(range1.min) + "-" + std::to_string(range1.max);
            matrix_intersection2->GetXaxis()->SetBinLabel(c, (label).c_str());
            matrix_intersection2->GetYaxis()->SetBinLabel(c, (label).c_str());
            c++;
        }
    }
    int kk = 1;
    for (const auto& range1 : selected2){
        int jj = 1;
         for (const auto& range2 : selected2){
             VecVariables intersection;
             std::set_intersection(range1.second.begin(), range1.second.end(), range2.second.begin(), range2.second.end(),
                                     std::back_inserter(intersection));
             matrix_intersection2->SetBinContent(kk, jj, intersection.size());
             jj++;
         }
         kk++;
    }
    root_ext::WriteObject(*matrix_intersection2, outfile.get());

    return vector_union;
}


void KolmogorovSignalPlotSelected(const MassVar& sample_signal, const VecVariables& selected,
                          int min, int max, int medium, TDirectory* directory){
    std::map<std::string, std::shared_ptr<TGraph>> plot;
     int i = 0;
    for (const auto& mass_entry: sample_signal){
        if ( mass_entry.first == Bkg_mu )
            continue;
        if (mass_entry.first < min || mass_entry.first > max)
            continue;
        for (const auto& var : selected){
            std::vector<double> vector_signal = sample_signal.at(mass_entry.first).at(var);
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
                plot[var]->SetTitle(("KolmogorovSmirnov_"+var+"_Signal"+std::to_string(min)+"_"+std::to_string(max)).c_str());
                plot[var]->SetName(("Ks_"+var+"_Signal"+std::to_string(min)+"_"+std::to_string(max)).c_str());
                plot[var]->GetHistogram()->GetXaxis()->SetTitle("mass");
                plot[var]->GetHistogram()->GetYaxis()->SetTitle("KS Probability");
            }
            plot[var]->SetPoint(i,mass_entry.first,k);
        }
        i++;
    }

    for(const auto& var: selected){
        root_ext::WriteObject(*plot[var], directory);
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
            if (line.size()==0 || line.at(0)=='#')
                continue;
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

        run::ThreadPull threads(8);
        auto start = clock::now();

        std::map<std::string, MassVar> sample_vars;

        for (int i = 0; i<2; i++){
            std::string tree;
            if (i == 0) tree = "muTau";
            if (i == 1) tree = "eTau";
            for(const SampleEntry& entry:samples)
            {
                if ( i == 0 && entry.filename == "TT_ext3_eTau.root")  continue;
                if ( i == 1 && entry.filename == "TT_ext3_muTau.root") continue;
                auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
                EventTuple tuple(tree, input_file.get(), true, {} , GetEnabledBranches());
                std::cout << entry << " number of events: " << std::min(tuple.GetEntries(),args.number_events()) << std::endl;
                Long64_t current_entry = 0;
                Long64_t tot_entries = 0;
                while(tot_entries < std::min(tuple.GetEntries(),args.number_events())) {
                    tuple.GetEntry(current_entry);
                    const Event& event = tuple.data();
                    if (event.eventEnergyScale != 0 || (event.q_1+event.q_2) != 0 || event.jets_p4.size() < 2
                        || event.extraelec_veto == true || event.extramuon_veto == true || event.jets_p4[0].eta() > 2.4
                        || event.jets_p4[1].eta() > 2.4){
                        current_entry++;
                        continue;
                    }
                    LorentzVectorE_Float bb = event.jets_p4[0] + event.jets_p4[1];
                    double ellipse_cut = pow(event.SVfit_p4.mass()-116,2)/pow(35.,2) + pow(bb.mass()-111,2)/pow(45.,2);
                    if (ellipse_cut>1){
                        current_entry++;
                        continue;
                    }
                    vars.AddEvent(tree, event, entry.mass, entry.weight);
                    current_entry++;
                    tot_entries++;
                }
            }
            auto stop = clock::now();
            std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;
            sample_vars[tree] = vars.GetSampleVariables(tree);
        }

        std::cout<<"quanti tipi di fondo?"<<sample_vars.size()<<std::endl;

        std::cout << "n.variabili: " << sample_vars.at("muTau").at(Bkg_mu).size() << std::endl;
        std::cout << "n.masse segnale: " << sample_vars.at("muTau").size() - 1 << " + " << 1 << " fondo."<<std::endl;

        for (const auto& sample_bkg: sample_vars){

            outfile->mkdir((sample_bkg.first).c_str());
            auto directory_bkg = outfile->GetDirectory((sample_bkg.first).c_str());

            for (const auto& var: sample_bkg.second){
                auto map = var.second;
                if (var.first != Bkg_mu && var.first != Signal_SM)
                    std::cout << "massa: " << var.first << ", eventi segnale: " << map["pt_l1"].size() << std::endl;
                else if (var.first == Signal_SM) std::cout << "eventi SM: "<< map["pt_l1"].size() << std::endl;
                else if (var.first == Bkg_mu) std::cout << "eventi fondo: "<< map["pt_l1"].size() << std::endl;
            }

            MassVarBand bandwidth;
            MassVarPairMI mutual_matrix;
            std::cout << "bandwidth and mutual information" << std::endl;
            directory_bkg->mkdir("MutualInformation");
            auto directory_mi = directory_bkg->GetDirectory("MutualInformation");
//            start = clock::now();
            for (const auto& sample: sample_bkg.second){
                std::cout<<sample.first<<" bandwidth ";
                bandwidth[sample.first] = OptimalBandwidth(sample.second);
                std::cout<<" mutual ";
                mutual_matrix[sample.first] = Mutual(sample.second, bandwidth[sample.first]);
                if (sample.first == Bkg_mu || !mutual_matrix.count(Bkg_mu))
                    continue;
                else {
                    std::cout<<" mutual plot ";
                    MutualPlot(sample.first, mutual_matrix[sample.first], mutual_matrix[Bkg_mu], directory_mi);
                }
                std::cout<<std::endl;
            }
            std::cout <<"mutual histos" << std::endl;
            directory_mi->mkdir("Matrix");
            auto directory_s = directory_mi->GetDirectory("Matrix");
            CreateMatrixHistos(sample_bkg.second,mutual_matrix,"MI", directory_s);
//            stop = clock::now();
//            std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;

            start = clock::now();
            std::cout<<"Selection variables"<<std::endl;
            MassSetSelected selected;
            RangeSelected selected2;
            std::vector<VecVariables> union_selected = VariablesSelection(sample_bkg.second,bandwidth,mutual_matrix, selected, selected2, outfile, args.number_variables());
//            stop = clock::now();
//            std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;
            std::cout<<"Kolmogorov"<<std::endl;
            outfile->mkdir("Kolmogorov");
            auto directory_ks = outfile->GetDirectory("Kolmogorov");
            int count = 0;
            for (const auto& range : ranges){
                KolmogorovSignalPlotSelected(sample_bkg.second,union_selected[0], range.min , range.max , range.med, directory_ks);
                count++;
            }

            std::cout<<"Correlation mass selected"<<std::endl;
            outfile->mkdir("Correlation");
            auto directory_corr = outfile->GetDirectory("Correlation");
            directory_corr->mkdir("MassSelected");
            auto directory_ms = directory_corr->GetDirectory("MassSelected");
            for (const auto& sample: sample_bkg.second){
                if (sample.first == Bkg_mu) continue;
                 if (selected.at(sample.first).size() == 0) continue;
                VarPairEstCorr matrix_covariance = EstimateCovariance(selected.at(sample.first), sample_bkg.second.at(sample.first), UINT_FAST32_MAX);
                auto matrix_covariance_value = GetValue(matrix_covariance);
                auto matrix_corr = CovToCorr(matrix_covariance_value);
                CreateMatrixHistosMassSelected(sample.first, matrix_corr, selected.at(sample.first), "correlation", directory_ms);
            }
            std::cout<<"Correlation range union"<<std::endl;
            directory_corr->mkdir("Range_Union");
            auto directory_ru = directory_corr->GetDirectory("Range_Union");
            int i = 0;
            for (const auto& range: ranges){
                std::set<std::string> selected;
                for(const auto& var : union_selected[i]){
                    selected.insert(var);
                }
                if (!selected.size()) continue;
                VarData vardata;
                for(const auto& var : union_selected[i]){
                    for (const auto& sample: sample_bkg.second){
                        if (sample.first < range.min || sample.first > range.max) continue;
                        for (const auto& element: sample_bkg.second.at(sample.first).at(var)){
                            vardata[var].push_back(element);
                        }
                    }
                }
                VarPairEstCorr matrix_covariance = EstimateCovariance(selected, vardata, UINT_FAST32_MAX);
                auto matrix_covariance_value = GetValue(matrix_covariance);
                auto matrix_corr = CovToCorr(matrix_covariance_value);
                CreateMatrixHistosMassSelected(range.min, matrix_corr, selected, "correlation", directory_ru);
                i++;
            }
            std::cout<<"Correlation range"<<std::endl;
            directory_corr->mkdir("Range");
            auto directory_r = directory_corr->GetDirectory("Range");
            int j = 0;
            for (const auto& range: ranges){
                VarData vardata;
                for(const auto& var : selected2.at(range.min)){
                    for (const auto& sample: sample_bkg.second){
                        if (sample.first < range.min || sample.first > range.max) continue;
                        for (const auto& element: sample_bkg.second.at(sample.first).at(var)){
                            vardata[var].push_back(element);
                        }
                    }
                }
                VarPairEstCorr matrix_covariance = EstimateCovariance(selected2.at(range.min), vardata, UINT_FAST32_MAX);
                auto matrix_covariance_value = GetValue(matrix_covariance);
                auto matrix_corr = CovToCorr(matrix_covariance_value);
                CreateMatrixHistosMassSelected(range.min, matrix_corr, selected2.at(range.min), "correlation", directory_r);
                j++;
            }
        }

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
