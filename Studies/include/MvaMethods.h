/*! Definition of functions for calulating most used quantities, such as correlation, mutual information distance
 * and bandwidth.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "hh-bbtautau/Studies//include/MvaConfiguration.h"
#include "hh-bbtautau/Studies/include/MvaVariablesStudy.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include <TCanvas.h>
#include <future>
#include "AnalysisTools/Run/include/MultiThread.h"
#include "AnalysisTools/Core/include/StatEstimators.h"
#include "TGraphErrors.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"

namespace  analysis {
namespace mva_study{

using SetNamesVar = std::set<std::string>;
using SampleIdSetNamesVar = std::map<SampleId, std::set<std::string>>;
using SampleIdNameSet =  std::map<SampleId, Name_ND>;
using NameElement =  std::map<Name_ND, double>;
using NameElementFuture =  std::map<Name_ND, std::future<double>>;
using SampleIdNameElement = std::map<SampleId, NameElement>;
using VectorName_ND = std::deque<std::pair<Name_ND, double>>;
using SampleIdVectorName_ND = std::map<SampleId, VectorName_ND>;
using SamplePair = std::pair<SampleId,SampleId>;
using SamplePairNameNDElement = std::map<SamplePair, std::map<Name_ND, double>>;

template<typename Key, typename Value>
std::map<Key, Value> CollectFutures(std::map<Key, std::future<Value>>& futures)
{
    std::map<Key, Value> values;
    for(auto& var : futures) {
        values[var.first] = var.second.get();
    }
    return values;
}

//Read a .csv file and return or bandwidth of JSDivergence
NameElement Read_csvfile(const std::string& filecsv, const std::unordered_set<std::string>& disabled_vars)
{
    NameElement bandwidth;
    std::ifstream file(filecsv);
    std::string line;
    while(std::getline(file,line))
    {
        std::stringstream lineStream(line);
        std::string cell;
        double value = 0;
        Name_ND name;
        while(std::getline(lineStream,cell,','))
        {
            bool check = TryParse(cell,value);
            if (!check && !disabled_vars.count(cell)) name.insert(cell);
        }
        bandwidth[name] = value;
    }
    return bandwidth;
}

//Calculate optimal bandwidth for each variable for a single value of mass
inline NameElement OptimalBandwidth(const VarData& sample)
{
    NameElementFuture bandwidth_future;
    for (const auto& var : sample){
        bandwidth_future[Name_ND{var.first}] = run::async(stat_estimators::OptimalBandwith<double>, std::cref(var.second), 0.01);
    }
    return CollectFutures(bandwidth_future);
}

//Create elements of mutual information matrix for a single value o f mass
inline NameElement Mutual(const VarData& sample_vars, NameElement& bandwidth)
{
   NameElementFuture matrix_future;
   for(auto var_1 = sample_vars.begin(); var_1 != sample_vars.end(); ++var_1){
       for(auto var_2 = var_1; var_2 != sample_vars.end(); ++var_2) {
           if (bandwidth.at(Name_ND{var_1->first}) == 0) bandwidth[Name_ND{var_1->first}] = 0.0001;
           if (bandwidth.at(Name_ND{var_2->first}) == 0) bandwidth[Name_ND{var_2->first}] = 0.0001;
           matrix_future[Name_ND{var_1->first, var_2->first}] = run::async(stat_estimators::ScaledMutualInformation<double>, std::cref(var_1->second),
                                                                          std::cref(var_2->second), bandwidth.at(Name_ND{var_1->first}), bandwidth.at(Name_ND{var_2->first}));
        }
    }
    return CollectFutures(matrix_future);
}

//Estimate elements of correlation matrix for selected variables
inline NameElement Correlation(const VarData& sample_vars)
{
    NameElementFuture corr_matrix_future;
    for(auto var_1 = sample_vars.begin(); var_1 != sample_vars.end(); ++var_1){
        for(auto var_2 = var_1; var_2 != sample_vars.end(); ++var_2) {
            corr_matrix_future[Name_ND{var_1->first, var_2->first}] = run::async(stat_estimators::Correlation<double>, var_1->second, var_2->second, 1.);
        }
    }
    return CollectFutures(corr_matrix_future);
}

//Create histos of mutual information for signal and background
inline void MutualHisto(const SampleId& mass, const NameElement& mutual_matrix_signal, const NameElement& mutual_matrix_bkg, TDirectory* directory)
{
    auto directory_1d = directory->GetDirectory("1D");
    if (directory_1d == nullptr) {
        directory_1d = root_ext::GetDirectory(*directory, "1D");
        auto histo = std::make_shared<TH1D>("Background","MutualInformation_Background", 50, 0, 1);
        histo->SetXTitle("MID");
        for (const auto& entry : mutual_matrix_bkg){
            histo->Fill(entry.second);
        }
        root_ext::WriteObject(*histo, directory_1d);
    }

    auto directory_2d = root_ext::GetDirectory(*directory, "2D");
    std::string value = "Signal_"+ ToString(mass);
    auto histo2d = std::make_shared<TH2D>((value+"_Background").c_str(),
                                          ("MutualInformation_"+value+"_Background").c_str(),50,0,1,50,0,1);
    histo2d->SetXTitle("MI Signal");
    histo2d->SetYTitle("MI Background");
    auto histo = std::make_shared<TH1D>(value.c_str(),("MutualInformation_"+value).c_str(),50,0,1);
    histo->SetXTitle("MID");

    for (const auto& entry : mutual_matrix_signal){
        if (!mutual_matrix_bkg.count(entry.first)) continue;
        histo2d->Fill(entry.second, mutual_matrix_bkg.at(entry.first));
        histo->Fill(entry.second);
    }
    root_ext::WriteObject(*histo, directory_1d);
    root_ext::WriteObject(*histo2d, directory_2d);
}

//Create correlation/mutual information/JSD histo matrix
inline void CreateMatrixHistos(const SampleIdVarData& samples_mass, const SampleIdNameElement& element,
                               const std::string& type, TDirectory* directory, bool draw_diagonal = false)
{
    for(const auto& mass_entry: samples_mass){
        if (!element.count(mass_entry.first)) continue;
        int bin = static_cast<int>(mass_entry.second.size());
        std::string value = type+"_"+ ToString(mass_entry.first);
        auto matrix = std::make_shared<TH2D>(value.c_str(), value.c_str(),
                                             bin, 0, bin, bin, 0, bin);
        int i = 1;
        for(auto var_1 = mass_entry.second.begin(); var_1 != mass_entry.second.end(); ++var_1){
            if (!element.at(mass_entry.first).count(var_1->first)) continue;
            int j = i + !draw_diagonal;
            matrix->GetXaxis()->SetBinLabel(i, (var_1->first).c_str());
            matrix->GetYaxis()->SetBinLabel(i, (var_1->first).c_str());
            for(auto var_2 = std::next(var_1, !draw_diagonal); var_2 != mass_entry.second.end(); ++var_2) {
                Name_ND var_12{var_1->first, var_2->first};
                matrix->SetBinContent(i, j, element.at(mass_entry.first).at(var_12));
                matrix->SetBinContent(j, i, element.at(mass_entry.first).at(var_12));
                j++;
            }
            i++;
        }
        root_ext::WriteObject(*matrix, directory);
    }
}

inline void CreateMatrixHistosCompatibility(const SampleIdVarData& samples_mass, const std::map<SamplePair, double>& values,
                                            const std::string& name, TDirectory* directory)
{
    int k = 1;
    int bin = static_cast<int>(samples_mass.size());
    auto  histo_sgn_compatibility = std::make_shared<TH2D>(name.c_str(),name.c_str(), bin, 0, bin, bin, 0, bin);
    for(auto mass_entry_1 = samples_mass.begin(); mass_entry_1 != samples_mass.end(); ++mass_entry_1) {
        if ( mass_entry_1->first.IsBackground())
            continue;
        histo_sgn_compatibility->GetXaxis()->SetBinLabel(k, (ToString(mass_entry_1->first)).c_str());
        histo_sgn_compatibility->GetYaxis()->SetBinLabel(k, (ToString(mass_entry_1->first)).c_str());
        int j = k;
        for(auto mass_entry_2 = mass_entry_1; mass_entry_2 != samples_mass.end(); ++mass_entry_2) {
            if ( mass_entry_2->first.IsBackground())
                continue;
            SamplePair mass_pair(mass_entry_1->first, mass_entry_2->first);
            histo_sgn_compatibility->SetBinContent(k, j, values.at(mass_pair));
            histo_sgn_compatibility->SetBinContent(j, k, values.at(mass_pair));
            j++;
        }
        k++;
    }
    root_ext::WriteObject(*histo_sgn_compatibility, directory);
}

template<typename Graph = TGraph>
inline std::shared_ptr<Graph> CreatePlot(const std::string& title, const std::string& name, const std::string& x_axis, const std::string& y_axis)
{
    auto plot = std::make_shared<Graph>();
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

inline NameElement JensenDivergenceSamples(const VarData& sample_signal, const VarData& sample_bkg,
                               const NameElement& bandwidth_signal, const NameElement& bandwidth_bkg)
{
    NameElementFuture JSDivergenceND_future;
    for(auto entry_1 = sample_signal.begin(); entry_1 != sample_signal.end(); ++entry_1) {
        std::vector<const DataVector*> x;
        std::vector<const DataVector*> y;
        DataVector band_x;
        DataVector band_y;
        x.push_back(&sample_signal.at(entry_1->first));
        y.push_back(&sample_bkg.at(entry_1->first));
        band_x.push_back(bandwidth_signal.at(Name_ND{entry_1->first}));
        band_y.push_back(bandwidth_bkg.at(Name_ND{entry_1->first}));        
        JSDivergenceND_future[Name_ND{entry_1->first}] = run::async(stat_estimators::JensenShannonDivergence_ND<double>, x, y, band_x, band_y);
        for(auto entry_2 = std::next(entry_1); entry_2 != sample_signal.end(); ++entry_2) {
            x.push_back(&sample_signal.at(entry_2->first));
            y.push_back(&sample_bkg.at(entry_2->first));
            band_x.push_back(bandwidth_signal.at(Name_ND{entry_2->first}));
            band_y.push_back(bandwidth_bkg.at(Name_ND{entry_2->first}));
            JSDivergenceND_future[Name_ND{entry_1->first, entry_2->first}] = run::async(stat_estimators::JensenShannonDivergence_ND<double>,x, y,
                                                                                     band_x, band_y);
            x.erase(x.end() - 1);
            y.erase(y.end() - 1);
            band_x.erase(band_x.end() - 1);
            band_y.erase(band_y.end() - 1);
        }
    }
    return CollectFutures(JSDivergenceND_future);
}

struct TimeReporter {
public:
    using clock = std::chrono::system_clock;
    TimeReporter() : start(clock::now()), start_tot(clock::now()) {}

    void TimeReport(bool tot = false)
    {

        auto stop = clock::now();
        auto delta = tot ? stop - start_tot : stop - start;
        std::cout<<"secondi";
        if(tot)
            std::cout << " totali";
        std::cout << ": " << std::chrono::duration_cast<std::chrono::seconds>(delta).count()<<std::endl;
        if(!tot)
            start = stop;
    }

private:
    clock::time_point start, start_tot;
};

inline VectorName_ND CopySelectedVariables(const VectorName_ND& JSDivergence_vector, const Name_ND& best_entry, const SetNamesVar& not_corrected)
{
    VectorName_ND copy;
    for (const auto& entry : JSDivergence_vector){
        if(entry.second <= 0) continue;
        bool consider = entry.first != best_entry;
        for(const auto& other : entry.first) {
            consider = consider && !not_corrected.count(other);
        }
        if(consider)
            copy.push_back(entry);
    }
    return copy;
}


class BDTData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    using Entry = root_ext::AnalyzerDataEntry<TH1D>;
    using Hist = Entry::Hist;

    TH1D_ENTRY(bdt_out, 101, -1.01, 1.01)
    TH1D_ENTRY(difference, 200, -1., 1.)
};

const SampleId mass_tot = SampleId::MassTot();
const SampleId bkg = SampleId::Bkg();
const std::string tot = "full";
const size_t test_train = 10; //test and training together in the evaluation of the method
const int spin_tot = 3;
const std::string all_channel = "all_channel";
const int SM_spin = 1;
const int bkg_spin = -1;

struct ChannelSampleIdSpin{
    std::string channel;
    SampleId sample_id;
    int spin;
    ChannelSampleIdSpin() : channel("all_channel"), sample_id(mass_tot), spin(spin_tot) {}
    ChannelSampleIdSpin(std::string _channel, SampleId _sample_id, int _spin) : channel(_channel), sample_id(_sample_id), spin(_spin) {}

    bool operator<(const ChannelSampleIdSpin& x) const
    {
        if (channel != x.channel) return channel < x.channel;
        if (sample_id != x.sample_id) return sample_id < x.sample_id;
        return spin < x.spin;
    }

    bool operator ==(const ChannelSampleIdSpin& x) const
    {
        return x.channel==channel && x.sample_id == sample_id && x.spin == spin ;
    }

    bool operator !=(const ChannelSampleIdSpin& x) const {return !(x==*this);}

    bool IsAllChannel() const { return channel == all_channel; }
    bool IsAllSpin() const { return spin == spin_tot; }

};

const ChannelSampleIdSpin AllSgn{all_channel, mass_tot, spin_tot};
const ChannelSampleIdSpin AllBkg{all_channel, bkg, spin_tot};


template<typename TestFn>
inline std::map<ChannelSampleIdSpin,double> DistributionCompatibilityTest(const std::string& which_test, const std::map<ChannelSampleIdSpin, std::map<size_t, std::vector<double>>>& evaluation,
                                                 BDTData::Entry& outputBDT, TDirectory* directory, bool ver, const TestFn& testFn)
{
    std::map<ChannelSampleIdSpin,double>  test;
    std::shared_ptr<TH1D> histo;
    histo = std::make_shared<TH1D>(which_test.c_str(), which_test.c_str(), 50, 0, 1.01);
    histo->SetXTitle(which_test.c_str());
    for (const auto& id : evaluation){
        test[id.first] = testFn(outputBDT(id.first.channel, id.first.sample_id, id.first.spin, 0),
                               outputBDT(id.first.channel, id.first.sample_id, id.first.spin, 1));
        histo->Fill(test[id.first]);
        if (ver)
            std::cout<<id.first.channel<<"  "<<id.first.sample_id.sampleType<<id.first.sample_id.mass<<"  "<<id.first.spin<<"   "<<test[id.first]<<std::endl;
    }
    histo->Scale(1/histo->Integral());
    root_ext::WriteObject(*histo, directory);
    return test;
}

inline std::map<ChannelSampleIdSpin,double> KolmogorovTest(const std::map<ChannelSampleIdSpin, std::map<size_t, std::vector<double>>>& evaluation,
                                            BDTData::Entry& outputBDT, TDirectory* directory, bool ver)
{
    auto test_fn = [](const TH1D& h1, const TH1D& h2) { return h1.KolmogorovTest(&h2, "X"); };
    return DistributionCompatibilityTest("KS", evaluation, outputBDT, directory, ver, test_fn);
}

inline std::map<ChannelSampleIdSpin,double>  ChiSquareTest(const std::map<ChannelSampleIdSpin, std::map<size_t, std::vector<double>>>& evaluation,
                                            BDTData::Entry& outputBDT, TDirectory* directory, bool ver)
{
    auto test_fn = [](const TH1D& h1, const TH1D& h2) { return h1.Chi2Test(&h2, "WW"); };
    return DistributionCompatibilityTest("chi^2 probability", evaluation, outputBDT, directory, ver, test_fn);
}

struct OptimalSignificance {
    double cut;
    PhysicalValue significance;
    OptimalSignificance() {}
    OptimalSignificance(double _cut, PhysicalValue _significance) : cut(_cut), significance(_significance) {}
};

inline std::vector<OptimalSignificance> Calculate_CutSignificance(const ChannelSampleIdSpin& id_sgn_mass, const ChannelSampleIdSpin& id_bkg_mass,
                                                                       const std::string& title, BDTData::Entry& outputBDT,
                                                                       TDirectory* directory){
     std::vector<OptimalSignificance> cuts;
     auto int_S_tot = Integral(outputBDT(id_sgn_mass.channel, id_sgn_mass.sample_id, id_sgn_mass.spin, tot), true);
     auto int_B_tot = Integral(outputBDT(id_bkg_mass.channel, id_bkg_mass.sample_id, id_bkg_mass.spin, tot), true);
     auto relative = int_S_tot/std::sqrt(int_B_tot);
     auto histo_sign = CreatePlot<TGraphErrors>(title, title,"output BDT","S/(sqrt(B))");
     int nbin = outputBDT(id_sgn_mass.channel, id_sgn_mass.sample_id, id_sgn_mass.spin, tot).GetNbinsX();
     for(int i = 0; i<=nbin; ++i){
         auto output = outputBDT(id_sgn_mass.channel, id_sgn_mass.sample_id, id_sgn_mass.spin, tot).GetBinCenter(i);
         auto int_S = Integral(outputBDT(id_sgn_mass.channel, id_sgn_mass.sample_id, id_sgn_mass.spin,tot), i, nbin+1);
         auto int_B = Integral(outputBDT(id_bkg_mass.channel, id_bkg_mass.sample_id, id_bkg_mass.spin, tot), i, nbin+1);
         PhysicalValue significance;
         static const PhysicalValue sign_inf(-1,0);
         if ( int_B.GetValue() != 0)
             significance = (int_S/std::sqrt(int_B))/relative;
         else
             significance = sign_inf;
         cuts.emplace_back(output, significance);
         histo_sign->SetPoint(i, output, significance.GetValue());
         histo_sign->SetPointError(i, 0, significance.GetFullError());
     }
     std::sort(cuts.begin(), cuts.end(), [](auto el1, auto el2){
         return el1.significance.GetValue() > el2.significance.GetValue();
     });
     root_ext::WriteObject(*histo_sign, directory);
     return cuts;
}

inline std::map<ChannelSampleIdSpin, OptimalSignificance> EstimateSignificativity(const std::string& channel, const int& spin, const std::vector<int>& mass_range,
                                                                                              BDTData::Entry& outputBDT, TDirectory* directory, bool total)

{
    std::map<ChannelSampleIdSpin, OptimalSignificance> sign;

    if (spin == SM_spin) {
        SampleId sgn_SM(SampleType::Sgn_NonRes, 0);
        SampleId bkg(SampleType::Bkg_TTbar, 0);
        ChannelSampleIdSpin id_sgn_SM{channel, sgn_SM, spin};
        ChannelSampleIdSpin id_bkg{channel, bkg, bkg_spin};
        auto cuts = Calculate_CutSignificance(id_sgn_SM, id_bkg, "SM"+channel+"_"+std::to_string(spin)+"_significance", outputBDT, directory);
        sign[id_sgn_SM] = cuts.front();

        ChannelSampleIdSpin id_sgn_SM_allchannel{all_channel, sgn_SM, spin};
        ChannelSampleIdSpin id_bkg_allchannel{all_channel, bkg, bkg_spin};
        cuts = Calculate_CutSignificance(id_sgn_SM_allchannel, id_bkg_allchannel, "SM_allchannel_"+std::to_string(spin)+"_significance", outputBDT, directory);
        sign[id_sgn_SM_allchannel] = cuts.front();
    }
    else{
        auto mass_sign = CreatePlot<TGraphErrors>("Mass_significance_"+channel+"_"+std::to_string(spin),"Mass_significance_"+channel+"_"+std::to_string(spin),"mass","S/(sqrt(B))");
        auto mass_sign_allchannel = CreatePlot<TGraphErrors>("Mass_significance_allchannel_"+std::to_string(spin),"Mass_significance_allchannel_"+std::to_string(spin),"mass","S/(sqrt(B))");
        int index = 0;
        for(const auto& mass : mass_range){
            SampleId sgn_mass(SampleType::Sgn_Res, mass);
            SampleId bkg_mass(SampleType::Bkg_TTbar, mass);
            ChannelSampleIdSpin id_sgn_mass{channel, sgn_mass, spin};
            ChannelSampleIdSpin id_bkg_mass{channel, bkg_mass, bkg_spin};
            auto cuts = Calculate_CutSignificance(id_sgn_mass, id_bkg_mass, std::to_string(mass)+channel+"_"+std::to_string(spin)+"_significance", outputBDT, directory);
            sign[id_sgn_mass] = cuts.front();
            mass_sign->SetPoint(index, mass, sign[id_sgn_mass].significance.GetValue());
            mass_sign->SetPointError(index, 0, sign[id_sgn_mass].significance.GetStatisticalError());
            ++index;

            ChannelSampleIdSpin id_sgn_mass_allchannel{all_channel, sgn_mass, spin};
            ChannelSampleIdSpin id_bkg_mass_allchannel{all_channel, bkg_mass, bkg_spin};
            cuts = Calculate_CutSignificance(id_sgn_mass_allchannel, id_bkg_mass_allchannel, std::to_string(mass)+"_allchannel_"+std::to_string(spin)+"_significance", outputBDT, directory);
            sign[id_sgn_mass] = cuts.front();
            mass_sign_allchannel->SetPoint(index, mass, sign[id_sgn_mass].significance.GetValue());
            mass_sign_allchannel->SetPointError(index, 0, sign[id_sgn_mass].significance.GetStatisticalError());
            ++index;
        }
        if (total){
            ChannelSampleIdSpin id_masstot{channel, mass_tot, spin};
            ChannelSampleIdSpin id_masstot_allchannel{all_channel, mass_tot, spin};
            ChannelSampleIdSpin id_bkg{channel,bkg, bkg_spin};
            ChannelSampleIdSpin id_bkg_allchannel{all_channel, bkg, bkg_spin};

            auto cuts = Calculate_CutSignificance(id_masstot, id_bkg, "Significance", outputBDT, directory);
            sign[id_masstot] = cuts.front();
            cuts = Calculate_CutSignificance(id_masstot_allchannel, id_bkg_allchannel, "Significance", outputBDT, directory);
            sign[id_masstot_allchannel] = cuts.front();
        }
        root_ext::WriteObject(*mass_sign, directory);
        root_ext::WriteObject(*mass_sign_allchannel, directory);
    }
    return sign;
}

}
}
