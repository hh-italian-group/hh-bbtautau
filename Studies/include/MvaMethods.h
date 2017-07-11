/*! Definition of functions for calulating most used quantities, such as correlation, mutual information distance
 * and bandwidth.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "hh-bbtautau/Studies//include/MvaConfiguration.h"
#include "hh-bbtautau/Analysis/include/MvaVariablesStudy.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include <TCanvas.h>
#include <future>
#include "AnalysisTools/Run/include/MultiThread.h"
#include "AnalysisTools/Core/include/StatEstimators.h"
#include "TGraphErrors.h"

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


template<typename Key, typename Value>
std::map<Key, Value> CollectFutures(std::map<Key, std::future<Value>>& futures)
{
    std::map<Key, Value> values;
    for(auto& var : futures) {
        values[var.first] = var.second.get();
    }
    return values;
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
inline NameElement Mutual(const VarData& sample_vars, const NameElement& bandwidth)
{
   NameElementFuture matrix_future;
   for(auto var_1 = sample_vars.begin(); var_1 != sample_vars.end(); ++var_1){
       for(auto var_2 = var_1; var_2 != sample_vars.end(); ++var_2) {
            matrix_future[Name_ND{var_1->first, var_2->first}] = run::async(stat_estimators::ScaledMutualInformation<double>, std::cref(var_1->second),
                                                                          std::cref(var_2->second), bandwidth.at(Name_ND{var_1->first}), bandwidth.at(Name_ND{var_2->first}));
        }
    }
    return CollectFutures(matrix_future);
}

//Estimate elements of covariance matrix for selected variables
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
        histo->SetXTitle("1-MI");
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
    histo->SetXTitle("1-MI");
    for (const auto& entry : mutual_matrix_signal){
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
    int bin = static_cast<int>(samples_mass.size()) - 1;
    auto  histo_sgn_compatibility = std::make_shared<TH2D>(name.c_str(),name.c_str(), bin, 0, bin, bin, 0, bin);
    for(auto mass_entry_1 = samples_mass.begin(); mass_entry_1 != samples_mass.end(); ++mass_entry_1) {
        if ( mass_entry_1->first.IsBackground())
            continue;
        histo_sgn_compatibility->GetXaxis()->SetBinLabel(k, (ToString(mass_entry_1->first)).c_str());
        histo_sgn_compatibility->GetYaxis()->SetBinLabel(k, (ToString(mass_entry_1->first)).c_str());
        int j = k;
        for(auto mass_entry_2 = std::next(mass_entry_1); mass_entry_2 != samples_mass.end(); ++mass_entry_2) {
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

inline std::shared_ptr<TGraph> CreatePlot(const std::string& title, const std::string& name, const std::string& x_axis, const std::string& y_axis)
{
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


inline std::shared_ptr<TGraphErrors> CreatePlotErrors(const std::string& title, const std::string& name, const std::string& x_axis, const std::string& y_axis)
{
    auto plot = std::make_shared<TGraphErrors>();
    plot->SetLineColor(kBlue);
    plot->SetLineWidth(1);
    plot->SetMarkerColor(1);
    plot->SetMarkerSize(1);
    plot->SetMarkerStyle(0);
    plot->SetTitle((title).c_str());
    plot->SetName((name).c_str());
    plot->GetHistogram()->GetXaxis()->SetTitle((x_axis).c_str());
    plot->GetHistogram()->GetYaxis()->SetTitle((y_axis).c_str());
    return plot;
}

inline NameElement JensenDivergenceSB(const VarData& sample_signal, const VarData& sample_bkg,
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

inline std::map<SampleId,double> Kolmogorov(const std::map<SampleId, std::map<size_t, std::vector<double>>>& evaluation, TDirectory* directory)
{
    std::map<SampleId,double> kolmogorov;
    std::shared_ptr<TH1D> histo_kolmogorov;
    histo_kolmogorov = std::make_shared<TH1D>("kolmogorov", "kolmogorov", 50, 0, 1.01);
    histo_kolmogorov->SetXTitle("KS");
    for (const auto& sample : evaluation){
        std::map<size_t, std::vector<double>> ks_vector;
        for (auto tvt : sample.second){
            auto type = tvt.first;
            std::sort(tvt.second.begin(), tvt.second.end());
            ks_vector[type] = std::move(tvt.second);
        }
        double ks = TMath::KolmogorovTest(static_cast<int>(ks_vector.at(0).size()), ks_vector.at(0).data(),
                                          static_cast<int>(ks_vector.at(1).size()), ks_vector.at(1).data(), "");
        std::cout<<sample.first.sampleType<<" "<<sample.first.mass<<"    "<<ks<<std::endl;
        kolmogorov[sample.first] = ks;
        histo_kolmogorov->Fill(ks);
    }
    root_ext::WriteObject(*histo_kolmogorov, directory);
    return kolmogorov;
}


const SampleId mass_tot = SampleId::MassTot();
const SampleId bkg = SampleId::Bkg();
constexpr int nbin = 202;
constexpr double bin_min = -1.01;
constexpr double bin_max = 1.01;
constexpr size_t tot = 10;

std::map<int, std::pair<double, PhysicalValue>> EstimateSignificativity(const std::vector<int>& mass_range,
                                                                        const std::map<SampleId, std::map<size_t,std::shared_ptr<TH1D>>>& outputBDT,
                                                                        TDirectory* directory, bool total)
{
    std::map<int, std::pair<double, PhysicalValue>> sign;
    auto mass_sign = CreatePlotErrors("Mass_significance","Mass_significance","mass","S/(sqrt(B))");

    int index = 0;
    for(const auto& mass : mass_range){
        std::vector<std::pair<double,PhysicalValue>> cuts;
        SampleId sgn_mass(SampleType::Sgn_Res, mass);
        SampleId bkg_mass(SampleType::Bkg_TTbar, mass);
        auto histo_sign = CreatePlotErrors((std::to_string(mass)+"_significance").c_str(),(std::to_string(mass)+"_significance").c_str(),"output BDT","S/(sqrt(B))");
        auto int_S_tot = Integral(*outputBDT.at(sgn_mass).at(tot).get(), true);
        auto int_B_tot = Integral(*outputBDT.at(bkg_mass).at(tot).get(), true);
        auto relative = int_S_tot/std::sqrt(int_B_tot);
        for(int i = 0; i<=nbin; ++i){
            auto output = outputBDT.at(sgn_mass).at(tot)->GetBinCenter(i);
            auto int_S = Integral(*outputBDT.at(sgn_mass).at(tot).get(), i, nbin+1);
            auto int_B = Integral(*outputBDT.at(bkg_mass).at(tot).get(), i, nbin+1);
            PhysicalValue significance;
            if ( int_S.GetValue() != 0 && int_B.GetValue() != 0) {
                significance = (int_S/std::sqrt(int_B))/relative;
            }
            cuts.emplace_back(output, significance);
            histo_sign->SetPoint(i, output, significance.GetValue());
            histo_sign->SetPointError(i, 0, significance.GetFullError());
        }
        std::sort(cuts.begin(), cuts.end(), [](auto el1, auto el2){
            return el1.second.GetValue() > el2.second.GetValue();
        });
        sign[mass] = cuts.front();
        mass_sign->SetPoint(index, mass, sign[mass].second.GetValue());
        mass_sign->SetPointError(index, 0, sign[mass].second.GetStatisticalError());
        std::cout<<index<<" "<<mass<<"  "<<ToString(sign[mass].second)<<std::endl;
        root_ext::WriteObject(*histo_sign, directory);
        ++index;
    }
    root_ext::WriteObject(*mass_sign, directory);


    if (total){
        std::vector<std::pair<double,PhysicalValue>> cuts;
        auto histo_sign = CreatePlotErrors("Significance","Significance","output BDT","S/(sqrt(B))");
        auto int_S_tot = Integral(*outputBDT.at(mass_tot).at(tot).get(), true);
        auto int_B_tot = Integral(*outputBDT.at(bkg).at(tot).get(), true);
        auto relative = int_S_tot/std::sqrt(int_B_tot);
        for(int i = 0; i<=nbin; ++i){
            auto output = outputBDT.at(mass_tot).at(tot)->GetBinCenter(i);
            auto int_S = Integral(*outputBDT.at(mass_tot).at(tot).get(), i, nbin+1);
            auto int_B = Integral(*outputBDT.at(bkg).at(tot).get(), i, nbin+1);
            PhysicalValue significance;
            if ( int_S.GetValue() != 0 && int_B.GetValue() != 0) {
                significance = (int_S/std::sqrt(int_B))/relative;
            }
            cuts.emplace_back(output, significance);
            histo_sign->SetPoint(i, output, significance.GetValue());
            histo_sign->SetPointError(i, 0, significance.GetFullError());
        }
        std::sort(cuts.begin(), cuts.end(), [](auto el1, auto el2){
            return el1.second.GetValue() > el2.second.GetValue();
        });
        sign[mass_tot.mass] = cuts.front();
        mass_sign->SetPoint(index, mass_tot.mass, sign[mass_tot.mass].second.GetValue());
        mass_sign->SetPointError(index, 0, sign[mass_tot.mass].second.GetStatisticalError());
        std::cout<<index<<" "<<mass_tot<<"  "<<ToString(sign[mass_tot.mass].second)<<std::endl;
        root_ext::WriteObject(*histo_sign, directory);
    }
    return sign;
}


}
}
