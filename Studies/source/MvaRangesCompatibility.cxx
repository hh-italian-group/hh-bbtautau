/*! Study of correlation matrix and mutual information of BDT variables
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <fstream>
#include <random>
#include <algorithm>
#include <functional>
#include <future>
#include <initializer_list>
#include <TCanvas.h>
#include <TPad.h>
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
    REQ_ARG(unsigned long, number_threads);
    REQ_ARG(size_t, number_variables);
    OPT_ARG(Long64_t, odd, -1);

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

static constexpr int Bkg = -1;
static constexpr int Signal_SM = 0;
static constexpr double threashold_mi = 0.7;

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
using MassSetSelected =  std::map<int,std::set<std::string>>;
using VarPairEstCorr = std::map<VarPair,stat_estimators::EstimatedQuantity>;
using VarPairCorr =  std::map<VarPair,double>;

class MvaVariablesStudy : public MvaVariables {
private:
    std::map<std::string, double> variables;
    MassVar all_variables;

public:
    using MvaVariables::MvaVariables;

    virtual void SetValue(const std::string& name, double value) override
    {
        variables[name] = value;
    }

    virtual void AddEventVariables(bool istraining, int mass,  double weight) override
    {
        VarData& sample_vars = all_variables[mass];
        for(const auto& name_value : variables) {
            const std::string& name = name_value.first;
            const double value = name_value.second;
            sample_vars[name].push_back(value);
        }
    }

    const MassVar& GetSampleVariables() const
    {
        return all_variables;
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
    Name_ND(std::initializer_list<std::string> _name) : names(_name.begin(), _name.end()){}

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
    int min, max;
    Range() : min(0), max(0){}
    bool Contains(int v) const { return v >= min && v <= max; }

    Range(Double_t _min, Double_t _max ) : min(_min), max(_max) {}
};

Range low_mass(250, 280);
Range medium_mass(300, 400);
Range high_mass(450,900);
std::vector<Range> ranges{low_mass, medium_mass, high_mass};

//Create optimal bandwidth for each variable for a single value of mass
VarBand OptimalBandwidth(const VarData& sample){
    VarBand bandwidth;
    std::map<std::string, std::future<double>> bandwidth_future;
        for (const auto& var : sample){
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

//Create elements of mutual information matrix for a single value o f mass
VarPairMI Mutual(const VarData& sample, const VarBand& bandwidth){
    VarPairMI mutual_matrix;
    std::map<VarPair,std::future<double>> matrix_future;
    for (const auto& var_1: sample){
        for(const auto& var_2 : sample) {
            if (var_2.first <= var_1.first) continue;
            const VarPair var_12(var_1.first, var_2.first);
            matrix_future[var_12] = run::async(stat_estimators::ScaledMutualInformation<double>,
                                               std::cref(var_1.second), std::cref(var_2.second),
                                               bandwidth.at(var_1.first), bandwidth.at(var_2.first));
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
    auto histo2d = std::make_shared<TH2D>(("Signal_mass"+m+"_Background").c_str(),
                                          ("MutualInformation_Signal_mass"+m+"_Background").c_str(),50,0,1,50,0,1);
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
VarPairEstCorr EstimateCovariance(const VarData& sample_vars, uint_fast32_t seed){
    VarPairEstCorr covariance_matrix;
    for(const auto& var_1 : sample_vars) {
        for(const auto& var_2 : sample_vars) {
            if (var_2.first < var_1.first) continue;
            const VarPair var_12(var_1.first, var_2.first);
            covariance_matrix[var_12] = stat_estimators::EstimateWithErrorsByResampling(stat_estimators::Covariance<double>,
                                                                                        sample_vars.at(var_1.first), sample_vars.at(var_2.first),
                                                                                        true, true, 1000, 0.31731, seed);
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
        if (mass_entry.first == Bkg)
            class_sample =  "Background";
        else {
            class_sample =  "Signal";
            if (mass_entry.first == Signal_SM) mass = "_SM";
            else mass = "_mass"+std::to_string(mass_entry.first);
        }
        auto matrix = std::make_shared<TH2D>((type+"_"+class_sample+mass).c_str(),(type+"_"+class_sample+mass).c_str(),
                                             bin, 0, bin, bin, 0, bin);
        int i = 1;
        for(const auto& var_1 : mass_entry.second) {
            int j = 1;
            matrix->GetXaxis()->SetBinLabel(i, (var_1.first).c_str());
            matrix->GetYaxis()->SetBinLabel(i, (var_1.first).c_str());
            for(const auto& var_2 : mass_entry.second) {
                const VarPair var_12(var_1.first, var_2.first);
                if (element.at(mass_entry.first).count(var_12)){
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

std::shared_ptr<TGraph> CreatePlot(std::string title, std::string name, std::string x_axis, std::string y_axis){
    auto plot = std::make_shared<TGraph>();
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
        JSDivergenceND_future[Name_ND{entry_1.first}] = run::async(stat_estimators::JensenShannonDivergence_ND<double>,
                                                                   x, y, band_x, band_y);
        for (const auto& entry_2 : sample_signal){
            if ( entry_1.first == entry_2.first )
                continue;
            if (JSDivergenceND_future.count(Name_ND({entry_2.first,entry_1.first})) )
                continue;
            x.push_back(&sample_signal.at(entry_2.first));
            y.push_back(&sample_bkg.at(entry_2.first));
            band_x.push_back(bandwidth_signal.at(entry_2.first));
            band_y.push_back(bandwidth_bkg.at(entry_2.first));
            JSDivergenceND_future[Name_ND{entry_1.first,entry_2.first}] = run::async(stat_estimators::JensenShannonDivergence_ND<double>,
                                                                                     x, y, band_x, band_y);
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

            for(const auto& other_entry : sample.at(mass)) {
                if(other_entry.first == name) continue;
                const Name_ND names({name, other_entry.first});
                const auto name_pair = names.ToPair();
                if(mutual_matrix_signal.at(name_pair) < threashold_mi && mutual_matrix_bkg.at(name_pair) < threashold_mi){
                    eliminated[other_entry.first] = "MI-"+name;
                    not_corrected.insert(other_entry.first);
                }
            }
            selected.insert(name);
        }
        JSDivergence_vector.pop_front();
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
    auto directory_v = directory->GetDirectory("Plot_Selected_variables");
    if (directory_v == nullptr) {
        directory->mkdir("Plot_Selected_variables");
        directory_v = directory->GetDirectory("Plot_Selected_variables");
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

    std::map<std::string, std::map<int, std::future<double>>> JSDss_selected_future;
    std::map<std::string, std::map<int, double>> JSDss_selected;
    std::map<std::string, std::shared_ptr<TGraph>> plot;
    for(const auto& selected_name : selected){
        plot[selected_name] = CreatePlot(("JSD_"+selected_name).c_str(),("JSD_"+selected_name).c_str(),"mass","JSD");
        std::vector<const DataVector*> sample_1;
        DataVector band_1;
        sample_1.push_back(&sample.at(mass).at(selected_name));
        band_1.push_back(band.at(mass).at(selected_name));
        for (const auto& mass_entry: sample){
            if (mass_entry.first == Bkg) continue;
            std::vector<const DataVector*> sample_2;
            DataVector band_2;
            sample_2.push_back(&sample.at(mass_entry.first).at(selected_name));
            band_2.push_back(band.at(mass_entry.first).at(selected_name));
            JSDss_selected_future[selected_name][mass_entry.first] = run::async(stat_estimators::JensenShannonDivergence_ND<double>, sample_1, sample_2, band_1, band_2);
        }
    }
    for(auto& var : JSDss_selected_future) {
        int i = 0;
        for(auto& mass : var.second) {
            mass.second.wait();
            if(!mass.second.valid())
                throw exception("future not valid");
            JSDss_selected[var.first][mass.first] = mass.second.get();
            plot.at(var.first)->SetPoint(i,mass.first,JSDss_selected.at(var.first).at(mass.first));
            i++;
        }
        root_ext::WriteObject(*plot.at(var.first), directory_m);
    }
    auto directory_distribution = directory->GetDirectory("Distribution");
    if (directory_distribution == nullptr) {
        directory->mkdir("Distribution");
        directory_distribution = directory->GetDirectory("Distribution");
    }
    std::string  massa = std::to_string(mass);
    TDirectory* directory_mass;
    if (mass == Signal_SM){
        directory_distribution->mkdir("Signal_SM");
        directory_mass = directory_distribution->GetDirectory("Signal_SM");
    }
    else{
    directory_distribution->mkdir(("Signal_mass"+massa).c_str());
    directory_mass = directory_distribution->GetDirectory(("Signal_mass"+massa).c_str());
    }
    if (mass == Signal_SM) massa = "SM";

    std::map<int, std::shared_ptr<TH1D>> histo;
    for (const auto& mass_entry: sample){
        if (mass_entry.first == Bkg) continue;
        histo[mass_entry.first] = std::make_shared<TH1D>(("JSD_Signal"+massa+"_Signal"+std::to_string(mass_entry.first)).c_str(),
                                                         ("JensenShannonDivergence_Signal"+massa+"_"+std::to_string(mass_entry.first)).c_str(), 20,0,1);
        histo.at(mass_entry.first)->SetXTitle("JSD");

    }
    for(const auto& selected_name : selected){
        for (const auto& mass_entry: sample){
            if (mass_entry.first == Bkg) continue;
            histo.at(mass_entry.first)->Fill(JSDss_selected.at(selected_name).at(mass_entry.first));
        }
    }
    for (const auto& mass_entry: sample){
        if (mass_entry.first == Bkg) continue;
        root_ext::WriteObject(*histo[mass_entry.first], directory_mass);
    }
    return selected;
}

void PlotJensenShannon(const MassName_ND& JSDivergenceND, TDirectory* directory){
    std::map<std::string, std::shared_ptr<TGraph>> plot;
    int i = 0;
    for (const auto& mass_entry: JSDivergenceND){
        if (mass_entry.first == Bkg)
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

void KolmogorovSignalCompatibility(int mass, const MassVar& sample_var, const std::set<std::string>& selected, TDirectory* directory){
    auto directory_v = directory->GetDirectory("Plot_Selected_variables");
    if (directory_v == nullptr) {
        directory->mkdir("Plot_Selected_variables");
        directory_v = directory->GetDirectory("Plot_Selected_variables");
    }
    TDirectory* directory_m;
    if (mass == Signal_SM){
        directory_v->mkdir("Signal_SM");
        directory_m = directory_v->GetDirectory("Signal_SM");
    }
    else{
    directory_v->mkdir(("Signal_mass"+std::to_string(mass)).c_str());
    directory_m = directory_v->GetDirectory(("Signal_mass"+std::to_string(mass)).c_str());
    }
    auto directory_distribution = directory->GetDirectory("Distribution");
    if (directory_distribution == nullptr) {
        directory->mkdir("Distribution");
        directory_distribution = directory->GetDirectory("Distribution");
    }
    std::string  massa = std::to_string(mass);
    if (mass == Signal_SM) massa = "SM";
    TDirectory* directory_mass;
    if (mass == Signal_SM){
        directory_distribution->mkdir("Signal_SM");
        directory_mass = directory_distribution->GetDirectory("Signal_SM");
    }
    else{
    directory_distribution->mkdir(("Signal_mass"+massa).c_str());
    directory_mass = directory_distribution->GetDirectory(("Signal_mass"+massa).c_str());
    }
    std::map<std::string, std::shared_ptr<TGraph>> plot;
    std::map<int, std::shared_ptr<TH1D>> histo;
    std::map<int, std::map<std::string, double>> kolmogorov;
    int i = 0;
    for (const auto& mass_entry: sample_var){
        if ( mass_entry.first == Bkg )
            continue;
        for (const auto& var : selected){
            std::vector<double> vector_signal = sample_var.at(mass_entry.first).at(var);
            std::sort(vector_signal.begin(), vector_signal.end());\
            std::vector<double> vector_signal_2 = sample_var.at(mass).at(var);
            std::sort(vector_signal_2.begin(), vector_signal_2.end());
            Double_t* v_s = vector_signal.data(), *v_s_2 = vector_signal_2.data();
            kolmogorov[mass_entry.first][var] = TMath::KolmogorovTest(vector_signal.size(), v_s, vector_signal_2.size(), v_s_2, "");
            if (plot.count(var) == 0) {
                plot[var] = CreatePlot(("KolmogorovSmirnov_"+var).c_str(),
                                       ("KolmogorovSmirnov_"+var).c_str(),
                                       "mass","KS Probability" );
            }
            plot[var]->SetPoint(i,mass_entry.first,kolmogorov.at(mass_entry.first).at(var));
        }
        i++;
        histo[mass_entry.first] = std::make_shared<TH1D>(("KolmogorovSmirnov_Signal"+massa+"_Signal"+std::to_string(mass_entry.first)).c_str(),
                                                         ("KolmogorovSmirnov_Signal"+massa+"_"+std::to_string(mass_entry.first)).c_str(), 20,0,1);
        histo.at(mass_entry.first)->SetXTitle("KS");

    }

    for (const auto& mass_entry: sample_var){
        if (mass_entry.first == Bkg) continue;
        for(const auto& selected_name : selected){
            histo.at(mass_entry.first)->Fill(kolmogorov.at(mass_entry.first).at(selected_name));
        }
    }
    for(const auto& var: selected){
        root_ext::WriteObject(*plot[var], directory_m);
    }
    for (const auto& mass_entry: sample_var){
        if (mass_entry.first == Bkg) continue;
        root_ext::WriteObject(*histo[mass_entry.first], directory_mass);
    }
}

using VecVariables =  std::vector<std::string>;
std::vector<VecVariables> VariablesSelection(MassVar sample, const MassName_ND& JSDivergenceSB,
                                             const MassVarBand& bandwidth, const MassVarPairMI& mutual_matrix,
                                             MassSetSelected& mass_selected, std::shared_ptr<TFile> outfile, size_t number_variables){

    outfile->mkdir("Kolmogorov");
    auto directory_ks = outfile->GetDirectory("Kolmogorov");

    auto directory_jensenshannon = outfile->GetDirectory("JensenShannonDivergence");
    auto directory_jenshan_sgnlbkg = directory_jensenshannon->GetDirectory("Signal_Background");
    std::cout<<"Checklist1"<<std::endl;
    directory_jenshan_sgnlbkg->mkdir("Distribution");
    auto directory_jenshan_distrib = directory_jenshan_sgnlbkg->GetDirectory("Distribution");
    directory_jensenshannon->mkdir("Signal_Signal");
    auto directory_jenshan_sgnlsgnl = directory_jensenshannon->GetDirectory("Signal_Signal");
    directory_jenshan_sgnlsgnl->mkdir("Check_Mass");
    auto directory_jenshan_mass = directory_jenshan_sgnlsgnl->GetDirectory("Check_Mass");

    std::ofstream ofs("SelectedVariable.csv", std::ofstream::out);
    for (const auto& mass_entry : sample){
        ofs << mass_entry.first << ",";
        std::cout<<mass_entry.first<<" ";
        std::cout.flush();
        if (mass_entry.first == Bkg)
            continue;
        std::map<std::string, std::string> eliminated;
        mass_selected[mass_entry.first] =  CheckList(mass_entry.first, sample,JSDivergenceSB.at(mass_entry.first),
                                                     mutual_matrix.at(mass_entry.first), mutual_matrix.at(Bkg),
                                                     eliminated, bandwidth, number_variables, directory_jenshan_mass);

        KolmogorovSignalCompatibility(mass_entry.first, sample, mass_selected[mass_entry.first], directory_ks);
        ofs << mass_entry.first << std::endl;
        for(auto& entry_selected: mass_selected[mass_entry.first]){
                ofs << entry_selected<<",";
        }
        ofs << std::endl;
        std::string  mass = std::to_string(mass_entry.first);

        if (mass_entry.first == Signal_SM) mass = "SM";
        auto histo = std::make_shared<TH1D>(("JSD_Signal"+mass+"_Background").c_str(), ("JensenShannonDivergence_Signal"+mass+"_Background").c_str(), 50,0,1);
        histo->SetXTitle("JSD");
        std::ofstream of("InformationTable_mass"+(std::to_string(mass_entry.first))+".csv", std::ofstream::out);
        of<<"Var_1"<<","<<"Var_2"<<","<<"JSD_ND"<<","<<"JSD_12-(JSD_1+JSD_2)"<<","<<"ScaledMI_Signal_12" <<","<<"ScaledMI_Bkg_12"<<","<<"selected"<<","<<","<<"eliminated by"<<std::endl;
        for(auto& entry : JSDivergenceSB.at(mass_entry.first)){
             histo->Fill(entry.second);
             if ((entry.first).names.size() == 1){
                 auto name_1 = ((entry.first).names).begin();
                 bool selected_1 = mass_selected[mass_entry.first].count(*name_1);
                 of<<*name_1<<","<<"   "<<","<<JSDivergenceSB.at(mass_entry.first).at(Name_ND({*name_1}))<<","<<","<<","<<","<<selected_1;
                 if(eliminated.count(*name_1)) of<<","<<","<<eliminated.at(*name_1);
                 of<<std::endl;
             }
             else  if ((entry.first).names.size() == 2){
                 auto name_1 = ((entry.first).names).begin();
                 auto name_2 = ((entry.first).names).rbegin();
                 bool selected_1 = mass_selected[mass_entry.first].count(*name_1);
                 bool selected_2 = mass_selected[mass_entry.first].count(*name_2);
                 of<<*name_1<<","<<*name_2<<","<<JSDivergenceSB.at(mass_entry.first).at(Name_ND({*name_1,*name_2}))<<","<<JSDivergenceSB.at(mass_entry.first).at(Name_ND({*name_1,*name_2}))-(JSDivergenceSB.at(mass_entry.first).at(Name_ND({*name_1}))+JSDivergenceSB.at(mass_entry.first).at(Name_ND({*name_2})))<<","<<mutual_matrix.at(mass_entry.first).at(VarPair(*name_1,*name_2))<<","<<mutual_matrix.at(Bkg).at(VarPair(*name_1,*name_2))<<","<<selected_1<<","<<selected_2;
                 if(eliminated.count(*name_1)) of<<","<<eliminated.at(*name_1);
                 if(eliminated.count(*name_2)) of<<","<<eliminated.at(*name_2);
                 of<<std::endl;
             }
        }
        of.close();
        root_ext::WriteObject(*histo, directory_jenshan_distrib);
    }
    std::cout<<std::endl<<"Intersection and Union"<<std::endl;
    int bin = mass_selected.size();
    auto matrix_intersection = std::make_shared<TH2D>("Intersection","Intersection of variables",bin, 0, bin, bin, 0, bin);
    matrix_intersection->SetXTitle("mass");
    matrix_intersection->SetYTitle("mass");
    int k = 1;
    std::vector<VecVariables> vec_union(3);
    for(const auto& mass_1: mass_selected){
        int j = 1;
        std::string mass;
        if (mass_1.first == Signal_SM) mass = "SM";
        else mass = std::to_string(mass_1.first);
        matrix_intersection->GetXaxis()->SetBinLabel(k, (mass).c_str());
        matrix_intersection->GetYaxis()->SetBinLabel(k, (mass).c_str());
        for(const auto& mass_2: mass_selected){
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
    std::ofstream of("UnionVariables.csv", std::ofstream::out);
    std::vector<VecVariables> vector_union;
    int count_range = 0;
    for (auto& names : vec_union){
        std::sort(names.begin(),names.end() );
        names.erase( std::unique( names.begin(), names.end() ), names.end() );
        of<<"Range"<<ranges[count_range].min<<"-"<<ranges[count_range].max<<std::endl;
        for (const auto& entry: names){
            of<<entry<<",";
        }
        of<<std::endl;
        vector_union.emplace_back(names);
        count_range++;
    }
    of.close();
    std::cout<<"#variabli per range dopo unione: "<<vec_union[0].size()<<" "<<vec_union[1].size()<<" "<<vec_union[2].size()<<std::endl;
    return vector_union;
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

        auto start_tot = clock::now();
        run::ThreadPull threads(args.number_threads());

        MassVar sample_vars;
        auto start = clock::now();
        for(const SampleEntry& entry:samples)
        {
            if ( args.tree_name() == "muTau" && entry.filename == "TT_ext3_eTau.root")  continue;
            if ( args.tree_name() == "eTau" && entry.filename == "TT_ext3_muTau.root") continue;
            auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
            EventTuple tuple(args.tree_name(), input_file.get(), true, {} , GetEnabledBranches());
            Long64_t current_entry = 0, entries = 0, tot_entries = 0;
            while(entries < std::min(tuple.GetEntries(), args.number_events()) && current_entry < tuple.GetEntries()) {
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
                current_entry++;
                entries++;
                if (args.odd() == -1) {
                    tot_entries++;
                    vars.AddEvent(event, entry.mass, entry.weight);
                }
                else if (args.odd() == 0 && entries%2 == 0) {
                    tot_entries++;
                    vars.AddEvent(event, entry.mass, entry.weight);
                }
                else if (args.odd() == 1 && entries%2 != 0) {
                    tot_entries++;
                    vars.AddEvent(event, entry.mass, entry.weight);
                }
            }
            std::cout << entry << " number of events: " << tot_entries << std::endl;
        }
        auto stop = clock::now();
        std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;

        sample_vars = vars.GetSampleVariables();

        std::cout << "n.variabili: " << sample_vars.at(Bkg).size() << std::endl;
        std::cout << "n.masse segnale: " << sample_vars.size() - 1 << " + " << 1 << " fondo."<<std::endl;
        std::cout << "n.eventi 250: " << sample_vars.at(250).at("pt_l1").size() << " + " << sample_vars.at(Bkg).at("pt_l1").size() <<std::endl;

        MassVarPairMI mutual_matrix;
        MassVarBand bandwidth;
        MassVarPairMI correlation_matrix;

        std::cout << "Bandwidth and Mutual Information" << std::endl;
        outfile->mkdir("MutualInformation");
        auto directory_mutinf = outfile->GetDirectory("MutualInformation");
        start = clock::now();
        for (const auto& sample: sample_vars){
            std::cout<<"----"<<sample.first<<"----"<<std::endl;
            std::cout<<"covariance  ";
            auto matrix_covariance = EstimateCovariance(sample.second, UINT_FAST32_MAX);
            auto matrix_covariance_value = GetValue(matrix_covariance);
            correlation_matrix[sample.first] = CovToCorr(matrix_covariance_value);
            std::cout<<"bandwidth  ";
            bandwidth[sample.first] = OptimalBandwidth(sample.second);
            std::cout<<"mutual  ";
            mutual_matrix[sample.first] = Mutual(sample.second, bandwidth.at(sample.first));
            if (sample.first == Bkg || !mutual_matrix.count(Bkg)){
                std::cout<<std::endl;
                continue;
            }
            else {
                std::cout<<"mutual plot  ";
                MutualPlot(sample.first, mutual_matrix.at(sample.first), mutual_matrix.at(Bkg), directory_mutinf);
            }
            std::cout<<std::endl;
        }
        std::cout <<"Mutual histos" << std::endl;
        directory_mutinf->mkdir("Matrix");
        auto directory_mutinf_matrix = directory_mutinf->GetDirectory("Matrix");
        CreateMatrixHistos(sample_vars, mutual_matrix, "MI", directory_mutinf_matrix);

        std::cout <<"Correlation histos" << std::endl;
        outfile->mkdir("Correlation");
        auto directory_correlation = outfile->GetDirectory("Correlation");
        CreateMatrixHistos(sample_vars, correlation_matrix, "correlation", directory_correlation);
        stop = clock::now();
        std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;

        outfile->mkdir("MI_Correlation");
        auto directory_mi_corr = outfile->GetDirectory("MI_Correlation");
        for(const auto& sample : sample_vars){
            std::string mass= "";
            if (sample.first == Signal_SM) mass = "_SM";
            else if (sample.first == Bkg) mass = "_Background";
            else mass = "_mass"+std::to_string(sample.first);
            auto histo = std::make_shared<TH2D>(("MI_Correlation_"+mass).c_str(),("MI_Correlation_"+mass).c_str(),50,0,1,50,0,1);
            histo->SetXTitle("MI");
            histo->SetYTitle("Correlation");
            for (const auto pair: mutual_matrix.at(sample.first)){
                histo->Fill(pair.second, std::abs(correlation_matrix.at(sample.first).at(pair.first)));
            }
            root_ext::WriteObject(*histo, directory_mi_corr);
            }


        MassName_ND JSDivergenceSB;
        std::cout<<"Jensen Shannon Signal Background";
        start = clock::now();
        for (const auto& mass_entry: sample_vars){
            if (mass_entry.first == Bkg)
                continue;
            JSDivergenceSB[mass_entry.first] = JensenDivergenceSB(mass_entry.second, sample_vars.at(Bkg),
                                                                  bandwidth.at(mass_entry.first), bandwidth.at(Bkg));
        }
        std::cout<<" - Plot Jensen Shannon"<<std::endl;
        outfile->mkdir("JensenShannonDivergence");
        auto directory_jensenshannon = outfile->GetDirectory("JensenShannonDivergence");
        directory_jensenshannon->mkdir("Signal_Background");
        auto directory_jenshan_sgnlbkg = directory_jensenshannon->GetDirectory("Signal_Background");
        directory_jenshan_sgnlbkg->mkdir("All_variables");
        auto directory_jenshan_allvars = directory_jenshan_sgnlbkg->GetDirectory("All_variables");
        PlotJensenShannon(JSDivergenceSB, directory_jenshan_allvars);
        stop = clock::now();
        std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;

        std::cout<<"Selection variables"<<std::endl;
        MassSetSelected mass_selected;
        std::vector<VecVariables> union_selected = VariablesSelection(sample_vars, JSDivergenceSB, bandwidth, mutual_matrix,
                                                                      mass_selected, outfile, args.number_variables());



        auto stop_tot = clock::now();
        std::cout<<"secondi totali: "<<std::chrono::duration_cast<std::chrono::seconds>(stop_tot - start_tot).count()<<std::endl;
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

