/*Configuration for Mva studies*/

#include <future>
#include "AnalysisTools/Core/include/StatEstimators.h"
#include "AnalysisTools/Run/include/MultiThread.h"

namespace  analysis {

static constexpr int Bkg = -1;
static constexpr int Signal_SM = 0;

using DataVector = std::vector<double>;
using VarData = std::map<std::string, DataVector>;
using MassVar = std::map<int, VarData>;

struct SampleEntry{
  std::string filename;
  double weight;
  int mass;
  std::string channel;
  SampleEntry() : weight(-1){}
};

std::ostream& operator<<(std::ostream& os, const SampleEntry& entry)
{
    os << entry.filename << " " << entry.mass << " " << entry.weight << " " <<  entry.channel;
    return os;
}

std::istream& operator>>(std::istream& is, SampleEntry& entry)
{
    is >> entry.filename >> entry.mass >> entry.weight >>entry.channel ;
    return is;
}

using SampleEntryCollection = std::vector<SampleEntry>;
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

class MvaVariablesStudy : public MvaVariables {
private:
    std::map<std::string, double> variables;
    std::map<bool, MassVar> all_variables;

public:
    using MvaVariables::MvaVariables;

    virtual void SetValue(const std::string& name, double value) override
    {
        variables[name] = value;
    }

    virtual void AddEventVariables(bool istraining, int mass,  double weight) override
    {
        VarData& sample_vars = all_variables[istraining][mass];
        for(const auto& name_value : variables) {
            const std::string& name = name_value.first;
            const double value = name_value.second;
            sample_vars[name].push_back(value);
        }
    }

    const MassVar& GetSampleVariables(bool istraining = false) const
    {
        return all_variables.at(istraining);
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

    const_iterator begin() const { return names.begin(); }
    const_iterator end() const { return names.end(); }
    size_t size() const { return names.size(); }

    bool IsEqual(const Name_ND& other) const
    {
        if (names.size() == other.size()){
            return (*names.begin() == *other.names.begin() && *names.rbegin() == *other.names.rbegin());
        }
        else return false;
    }
};

using DataVector = std::vector<double>;
using VarData = std::map<std::string, DataVector>;
using MassVar = std::map<int, VarData>;

using VarNameSet = std::set<std::string>;
using MassVarNameSet =  std::map<int, VarNameSet>;
using NameElement =  std::map<Name_ND, double>;
using MassNameElement = std::map<int, NameElement>;
using VectorName_ND = std::deque<std::pair<Name_ND, double>>;
using MassVectorName_ND = std::map<int, VectorName_ND>;

//Calculate optimal bandwidth for each variable for a single value of mass
NameElement OptimalBandwidth(const VarData& sample){
    NameElement bandwidth;
    std::map<Name_ND, std::future<double>> bandwidth_future;
        for (const auto& var : sample){
            bandwidth_future[Name_ND{var.first}] = run::async(stat_estimators::OptimalBandwith<double>, std::cref(var.second), 0.01);
        }
        for(auto& var : bandwidth_future) {
            if(!var.second.valid())
                throw exception("future not valid");
            bandwidth[Name_ND{var.first}] = var.second.get();
        }
    return bandwidth;
}

//Create elements of mutual information matrix for a single value o f mass
NameElement Mutual(const VarData& sample, const NameElement& bandwidth){
    NameElement mutual_matrix;
    std::map<Name_ND, std::future<double>> matrix_future;
    for (const auto& var_1: sample){
        for(const auto& var_2 : sample) {
            if (var_2.first <= var_1.first) continue;
            matrix_future[Name_ND{var_1.first, var_2.first}] = run::async(stat_estimators::ScaledMutualInformation<double>, std::cref(var_1.second),
                                                                          std::cref(var_2.second), bandwidth.at(Name_ND{var_1.first}), bandwidth.at(Name_ND{var_2.first}));
        }
    }
    for(auto& var : matrix_future) {
        if(!var.second.valid())
            throw exception("future not valid");
        mutual_matrix[var.first] = var.second.get();
    }
    return mutual_matrix;
}

//Create histos of mutual information for signal and background
void MutualHisto(int mass, const NameElement& mutual_matrix_signal, const NameElement& mutual_matrix_bkg, TDirectory* directory){
    auto directory_1d = directory->GetDirectory("1D");
    if (directory_1d == nullptr) {
        directory->mkdir("1D");
        directory_1d = directory->GetDirectory("1D");
        auto histo = std::make_shared<TH1D>("Background","MutualInformation_Background",50,0,1);
        histo->SetXTitle("MI");
        for (const auto& entry : mutual_matrix_bkg){
            histo->Fill(entry.second);
        }
        root_ext::WriteObject(*histo, directory_1d);
    }
    auto directory_2d = directory->GetDirectory("2D");
    if (directory_2d == nullptr) {
        directory->mkdir("2D");
        directory_2d = directory->GetDirectory("2D");
    }
    std::string value;
    if ( mass ==  Signal_SM) value = "SM";
    else value = std::to_string(mass);
    auto histo2d = std::make_shared<TH2D>(("Signal_mass"+value+"_Background").c_str(),
                                          ("MutualInformation_Signal_mass"+value+"_Background").c_str(),50,0,1,50,0,1);
    histo2d->SetXTitle("MI Signal");
    histo2d->SetYTitle("MI Background");
    auto histo = std::make_shared<TH1D>(("Signal_mass"+value).c_str(),("MutualInformation_Signal_mass"+value).c_str(),50,0,1);
    histo->SetXTitle("MI");
    for (const auto& entry : mutual_matrix_signal){
        histo2d->Fill(entry.second, mutual_matrix_bkg.at(entry.first));
        histo->Fill(entry.second);
    }
    root_ext::WriteObject(*histo, directory_1d);
    root_ext::WriteObject(*histo2d, directory_2d);
}

//Estimate elements of covariance matrix for selected variables
NameElement Correlation(const VarData& sample_vars){
    NameElement corr_matrix;
    for(const auto& var_1 : sample_vars) {
        for(const auto& var_2 : sample_vars) {
            if (var_2.first < var_1.first) continue;
            corr_matrix[Name_ND{var_1.first, var_2.first}] = 100 * stat_estimators::Correlation(sample_vars.at(var_1.first), sample_vars.at(var_2.first));
        }
    }
    return corr_matrix;
}

//Create correlation/mutual information/JSD histo matrix
void CreateMatrixHistos(const MassVar sample, const MassNameElement& element, std::string type, TDirectory* directory){
    std::string mass, class_sample;
    for(const auto mass_entry: sample){
        if (!element.count(mass_entry.first)) continue;
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
                Name_ND var_12{var_1.first, var_2.first};
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
        plot->SetMarkerColor(8);
        plot->SetMarkerSize(1);
        plot->SetMarkerStyle(3);
    plot->SetTitle((title).c_str());
    plot->SetName((name).c_str());
    plot->GetHistogram()->GetXaxis()->SetTitle((x_axis).c_str());
    plot->GetHistogram()->GetYaxis()->SetTitle((y_axis).c_str());
    return plot;
}

NameElement JensenDivergenceSB(const VarData& sample_signal, const VarData& sample_bkg,
                                   const NameElement& bandwidth_signal, const NameElement& bandwidth_bkg){
    NameElement  JSDivergenceSB;
    std::map<Name_ND,std::future<double>> JSDivergenceND_future;
    for (const auto& entry_1 : sample_signal){
        std::vector<const DataVector*> x;
        std::vector<const DataVector*> y;
        DataVector band_x;
        DataVector band_y;
        x.push_back(&sample_signal.at(entry_1.first));
        y.push_back(&sample_bkg.at(entry_1.first));
        band_x.push_back(bandwidth_signal.at(Name_ND{entry_1.first}));
        band_y.push_back(bandwidth_bkg.at(Name_ND{entry_1.first}));
        JSDivergenceND_future[Name_ND{entry_1.first}] = run::async(stat_estimators::JensenShannonDivergence_ND<double>, x, y, band_x, band_y);
        for (const auto& entry_2 : sample_signal){
            if ( entry_1.first == entry_2.first )
                continue;
            if (JSDivergenceND_future.count(Name_ND({entry_2.first,entry_1.first})) )
                continue;
            x.push_back(&sample_signal.at(entry_2.first));
            y.push_back(&sample_bkg.at(entry_2.first));
            band_x.push_back(bandwidth_signal.at(Name_ND{entry_2.first}));
            band_y.push_back(bandwidth_bkg.at(Name_ND{entry_2.first}));
            JSDivergenceND_future[Name_ND{entry_1.first,entry_2.first}] = run::async(stat_estimators::JensenShannonDivergence_ND<double>,x, y,
                                                                                     band_x, band_y);
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

}
