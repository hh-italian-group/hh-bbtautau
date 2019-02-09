/*! Definition of functions for calulating most used quantities, such as correlation, mutual information distance
 * and bandwidth.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <future>
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "hh-bbtautau/Studies/include/MvaVariablesStudy.h"

namespace  analysis {
namespace mva_study {

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
NameElement Read_csvfile(const std::string& filecsv, const std::unordered_set<std::string>& disabled_vars);

//Calculate optimal bandwidth for each variable for a single value of mass
NameElement OptimalBandwidth(const VarData& sample);

//Create elements of mutual information matrix for a single value o f mass
NameElement Mutual(const VarData& sample_vars, NameElement& bandwidth);

//Estimate elements of correlation matrix for selected variables
NameElement Correlation(const VarData& sample_vars);

//Create histos of mutual information for signal and background
void MutualHisto(const SampleId& mass, const NameElement& mutual_matrix_signal, const NameElement& mutual_matrix_bkg,
                 TDirectory* directory);

//Create correlation/mutual information/JSD histo matrix
void CreateMatrixHistos(const SampleIdVarData& samples_mass, const SampleIdNameElement& element,
                        const std::string& type, TDirectory* directory, bool draw_diagonal = false);

void CreateMatrixHistosCompatibility(const SampleIdVarData& samples_mass, const std::map<SamplePair, double>& values,
                                     const std::string& name, TDirectory* directory);

template<typename Graph>
std::shared_ptr<Graph> CreatePlot(const std::string& title, const std::string& name, const std::string& x_axis,
                                  const std::string& y_axis)
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

NameElement JensenDivergenceSamples(const VarData& sample_signal, const VarData& sample_bkg,
                                    const NameElement& bandwidth_signal, const NameElement& bandwidth_bkg);

struct TimeReporter {
public:
    using clock = std::chrono::system_clock;
    TimeReporter();

    void TimeReport(bool tot = false);

private:
    clock::time_point start, start_tot;
};

VectorName_ND CopySelectedVariables(const VectorName_ND& JSDivergence_vector, const Name_ND& best_entry,
                                    const SetNamesVar& not_corrected);

class BDTData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    using Entry = root_ext::AnalyzerDataEntry<TH1D>;
    using Hist = Entry::Hist;

    TH1D_ENTRY(bdt_out, 101, -1.01, 1.01)
    TH1D_ENTRY(difference, 200, -1., 1.)
};

extern const SampleId mass_tot;
extern const SampleId bkg;
extern const std::string tot;
constexpr size_t test_train = 10; //test and training together in the evaluation of the method
constexpr int spin_tot = 3;
extern const std::string all_channel;
constexpr int SM_spin = 1;
constexpr int bkg_spin = -1;

struct ChannelSampleIdSpin{
    std::string channel;
    SampleId sample_id;
    int spin;
    ChannelSampleIdSpin();
    ChannelSampleIdSpin(std::string _channel, SampleId _sample_id, int _spin);
    bool operator<(const ChannelSampleIdSpin& x) const;
    bool operator ==(const ChannelSampleIdSpin& x) const;
    bool operator !=(const ChannelSampleIdSpin& x) const;
    bool IsAllChannel() const;
    bool IsAllSpin() const;
};

extern const ChannelSampleIdSpin AllSgn;
extern const ChannelSampleIdSpin AllBkg;

template<typename TestFn>
std::map<ChannelSampleIdSpin,double> DistributionCompatibilityTest(const std::string& which_test,
                                                                   const std::map<ChannelSampleIdSpin,
                                                                   std::map<size_t, std::vector<double>>>& evaluation,
                                                                   BDTData::Entry& outputBDT, TDirectory* directory,
                                                                   bool ver, const TestFn& testFn)
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

std::map<ChannelSampleIdSpin,double> KolmogorovTest(const std::map<ChannelSampleIdSpin,
                                                    std::map<size_t, std::vector<double>>>& evaluation,
                                                    BDTData::Entry& outputBDT, TDirectory* directory, bool ver);

std::map<ChannelSampleIdSpin,double> ChiSquareTest(const std::map<ChannelSampleIdSpin,
                                                   std::map<size_t, std::vector<double>>>& evaluation,
                                                   BDTData::Entry& outputBDT, TDirectory* directory, bool ver);

struct OptimalSignificance {
    double cut;
    PhysicalValue significance;
};

std::vector<OptimalSignificance> Calculate_CutSignificance(const ChannelSampleIdSpin& id_sgn_mass,
                                                           const ChannelSampleIdSpin& id_bkg_mass,
                                                           const std::string& title, BDTData::Entry& outputBDT,
                                                           TDirectory* directory);

std::map<ChannelSampleIdSpin, OptimalSignificance> EstimateSignificativity(const std::string& channel, const int& spin,
                                                                           const std::vector<int>& mass_range,
                                                                           BDTData::Entry& outputBDT,
                                                                           TDirectory* directory, bool total);

}
}
