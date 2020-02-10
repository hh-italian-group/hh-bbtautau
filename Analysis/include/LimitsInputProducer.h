/*! Definition of the limits input producer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "EventAnalyzerDataCollection.h"
#include "SampleDescriptor.h"

namespace analysis {

class LimitsInputProducer {
public:
    using AnaData = ::analysis::EventAnalyzerData;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection;
    using Sample = ::analysis::SampleDescriptorBase;
    using SampleWP = Sample::Point;
    using SampleWPCollection = std::map<std::string, SampleWP>;
    using Hist = TH1D;
    using HistPtr = std::shared_ptr<root_ext::SmartHistogram<Hist>>;

    static std::string FullDataCardName(const std::string& datacard_name, UncertaintySource unc_source,
                                        UncertaintyScale unc_scale);
    static std::string EventRegionSuffix(EventRegion region);

    template<typename ...Args>
    LimitsInputProducer(AnaDataCollection& _anaDataCollection, Args&&... sample_descriptors) :
        anaDataCollection(&_anaDataCollection)
    {
        CollectWorkingPoints(std::forward<Args>(sample_descriptors)...);
    }


    void Produce(const std::string& outputFileNamePrefix, const std::string& setup_name,
                 const std::map<EventCategory, std::string>& eventCategories, EventSubCategory eventSubCategory,
                 const std::set<UncertaintySource>& uncertaintySources, const EventRegionSet& eventRegions,
                 const std::map<SelectionCut, std::string>& sel_aliases);

private:
    void CollectWorkingPoints() {}

    template<typename SampleCollection, typename ...Args>
    void CollectWorkingPoints(const SampleCollection& sample_descriptors, Args&&... other_sample_descriptors)
    {
        for(const auto& sample : sample_descriptors) {
            for(const SampleWP& wp : sample.second.working_points) {
                if(wp.datacard_name.size())
                    sampleWorkingPoints[wp.full_name] = wp;
            }
        }
        CollectWorkingPoints(std::forward<Args>(other_sample_descriptors)...);
    }

    bool CanHaveEmptyHistogram(const EventAnalyzerDataId& id, bool& print_warning) const;

private:
    AnaDataCollection* anaDataCollection;
    SampleWPCollection sampleWorkingPoints;
};

} // namespace analysis
