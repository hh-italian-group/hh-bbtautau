/*! Definition of functions for calulating most used quantities, such as correlation, mutual information distance
 * and bandwidth.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/MvaConfiguration.h"
#include "AnalysisTools/Core/include/RootExt.h"

#pragma once

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


inline NameElement CollectFutures(NameElementFuture& futures)
{
    NameElement values;
    for(auto& var : futures) {
        values[var.first] = var.second.get();
    }
    return values;
}

//Calculate optimal bandwidth for each variable for a single value of mass
NameElement OptimalBandwidth(const VarData& sample){
    NameElementFuture bandwidth_future;
    for (const auto& var : sample){
        bandwidth_future[Name_ND{var.first}] = run::async(stat_estimators::OptimalBandwith<double>, std::cref(var.second), 0.01);
    }
    return CollectFutures(bandwidth_future);
}

//Create elements of mutual information matrix for a single value o f mass
NameElement Mutual(const VarData& sample_vars, const NameElement& bandwidth){
   NameElementFuture matrix_future;
   for(auto var_1 = sample_vars.begin(); var_1 != sample_vars.end(); ++var_1){
       for(auto var_2 = std::next(var_1); var_2 != sample_vars.end(); ++var_2) {
            if (var_2->first <= var_1->first) continue;
            matrix_future[Name_ND{var_1->first, var_2->first}] = run::async(stat_estimators::ScaledMutualInformation<double>, std::cref(var_1->second),
                                                                          std::cref(var_2->second), bandwidth.at(Name_ND{var_1->first}), bandwidth.at(Name_ND{var_2->first}));
        }
    }
    return CollectFutures(matrix_future);
}

//Estimate elements of covariance matrix for selected variables
NameElement Correlation(const VarData& sample_vars){
    NameElementFuture corr_matrix_future;
    for(auto var_1 = sample_vars.begin(); var_1 != sample_vars.end(); ++var_1){
        for(auto var_2 = var_1; var_2 != sample_vars.end(); ++var_2) {
            corr_matrix_future[Name_ND{var_1->first, var_2->first}] = run::async(stat_estimators::Correlation<double>, var_1->second, var_2->second, 1.);
        }
    }
    return CollectFutures(corr_matrix_future);
}

//Create histos of mutual information for signal and background
void MutualHisto(const SampleId& mass, const NameElement& mutual_matrix_signal, const NameElement& mutual_matrix_bkg, TDirectory* directory){
    auto directory_1d = directory->GetDirectory("1D");
    if (directory_1d == nullptr) {
        directory_1d = root_ext::GetDirectory(*directory, "1D");
        auto histo = std::make_shared<TH1D>("Background","MutualInformation_Background", 50, 0, 1);
        histo->SetXTitle("MI");
        for (const auto& entry : mutual_matrix_bkg){
            histo->Fill(entry.second);
        }
        root_ext::WriteObject(*histo, directory_1d);
    }
    auto directory_2d = root_ext::GetDirectory(*directory, "2D");
    std::string value = "_"+ ToString(mass);
    auto histo2d = std::make_shared<TH2D>(("Signal_"+value+"_Background").c_str(),
                                          ("MutualInformation_Signal"+value+"_Background").c_str(),50,0,1,50,0,1);
    histo2d->SetXTitle("MI Signal");
    histo2d->SetYTitle("MI Background");
    auto histo = std::make_shared<TH1D>(("Signal"+value).c_str(),("MutualInformation_Signal"+value).c_str(),50,0,1);
    histo->SetXTitle("MI");
    for (const auto& entry : mutual_matrix_signal){
        histo2d->Fill(entry.second, mutual_matrix_bkg.at(entry.first));
        histo->Fill(entry.second);
    }
    root_ext::WriteObject(*histo, directory_1d);
    root_ext::WriteObject(*histo2d, directory_2d);
}

//Create correlation/mutual information/JSD histo matrix
void CreateMatrixHistos(const SampleIdVarData sample_vars, const SampleIdNameElement& element, const std::string& type, TDirectory* directory){
    std::string mass;
    for(const auto& mass_entry: sample_vars){
        if (!element.count(mass_entry.first)) continue;
        int bin =  mass_entry.second.size();
        auto matrix = std::make_shared<TH2D>((type+"_"+ ToString(mass_entry.first)).c_str(),(type+"_"+ToString(mass_entry.first)).c_str(),
                                             bin, 0, bin, bin, 0, bin);
        int i = 1;
        for(auto var_1 = mass_entry.second.begin(); var_1 != mass_entry.second.end(); ++var_1){
            int j = i;
            matrix->GetXaxis()->SetBinLabel(i, (var_1->first).c_str());
            matrix->GetYaxis()->SetBinLabel(i, (var_1->first).c_str());
            for(auto var_2 = var_1; var_2 != mass_entry.second.end(); ++var_2) {
                Name_ND var_12{var_1->first, var_2->first};
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

std::shared_ptr<TGraph> CreatePlot(const std::string& title, const std::string& name, const std::string& x_axis, const std::string& y_axis){
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

NameElement JensenDivergenceSB(const VarData& sample_signal, const VarData& sample_bkg,
                               const NameElement& bandwidth_signal, const NameElement& bandwidth_bkg){
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
}
}
