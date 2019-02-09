/*! Definition of histogram containers for flat tree analyzers.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnaTuple.h"
#include "SampleDescriptor.h"

namespace analysis {

class EventAnalyzerData : public root_ext::AnalyzerData {
public:
    using Entry = root_ext::AnalyzerDataEntry<TH1D>;
    using EntryPtr = std::shared_ptr<Entry>;
    using ValueType = float;
    using EntrySource = std::pair<EntryPtr, const ValueType*>;
    using HistContainer = std::map<std::string, EntrySource>;
    using HistDesc = PropertyConfigReader::Item;
    using HistDescCollection = PropertyConfigReader::ItemCollection;
    using SampleUnc = ModellingUncertainty::SampleUnc;

    EventAnalyzerData(std::shared_ptr<TFile> outputFile, const std::string& directoryName, Channel channel,
                      const EventAnalyzerDataId& dataId, const bbtautau::AnaTuple* anaTuple,
                      const std::set<std::string>& histogram_names, const HistDescCollection& descriptors,
                      const SampleUnc& sample_unc);
    void Fill(double weight);

private:
    static const HistDesc& FindDescriptor(const std::string& h_name, Channel channel, const EventAnalyzerDataId& dataId,
                                          const HistDescCollection& descriptors);

private:
    HistContainer histograms;
};

} // namespace analysis
