/*! Definition of the stacked prots producer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Print/include/DrawOptions.h"
#include "EventAnalyzerDataCollection.h"
#include "SampleDescriptor.h"

namespace analysis {

class StackedPlotsProducer {
public:
    using AnaData = ::analysis::EventAnalyzerData;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection;
    using Sample = ::analysis::SampleDescriptorBase;
    using SampleCollection = std::vector<const Sample*>;
    using Hist = TH1D;
    using HistPtr = std::shared_ptr<root_ext::SmartHistogram<Hist>>;
    using PlotConfigReader = ::analysis::PropertyConfigReader;
    using PlotConfig = PlotConfigReader::ItemCollection;
    using PageOptions = ::root_ext::draw_options::Page;

    static SampleCollection CreateOrderedSampleCollection(const std::vector<std::string>& draw_sequence,
                                                          const SampleDescriptorCollection& samples,
                                                          const CombinedSampleDescriptorCollection& combined_samples,
                                                          const std::vector<std::string>& signals,
                                                          const std::vector<std::string>& data, Channel channel);

    StackedPlotsProducer(AnaDataCollection& _anaDataCollection, const SampleCollection& _samples,
                         const std::string& plot_cfg_name, const std::string& page_opt_name,
                         const std::set<std::string>& _histogramNames = {});

    Channel ChannelId() const;
    const std::string& ChannelNameLatex() const;

    void PrintStackedPlots(const std::string& outputFileNamePrefix, const EventRegion& eventRegion,
                           const EventCategorySet& eventCategories,
                           const EventSubCategorySet& eventSubCategories, const std::set<std::string>& signals,
                           const Sample* total_bkg = nullptr);

private:
    HistPtr GetHistogram(const EventAnalyzerDataId& metaId, const std::string& sample_name,
                         const std::string& hist_name) const;

private:
    AnaDataCollection* anaDataCollection;
    SampleCollection samples;
    std::set<std::string> histogramNames;
    PlotConfig plot_cfg;
    PageOptions page_opt;
};

} // namespace analysis
