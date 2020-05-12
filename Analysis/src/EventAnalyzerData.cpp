/*! Definition of histogram containers for flat tree analyzers.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "hh-bbtautau/Analysis/include/EventAnalyzerData.h"

namespace analysis {

EventAnalyzerData::EventAnalyzerData(std::shared_ptr<TFile> outputFile, const std::string& directoryName,
                                     Channel channel, const EventAnalyzerDataId& dataId,
                                     const std::set<std::string>& histogram_names,
                                     const HistDescCollection& descriptors, const SampleUnc& sample_unc,
                                     bool readMode) :
    AnalyzerData(outputFile, directoryName, readMode)
{
    for(const auto& h_name : histogram_names) {
        const auto& desc = FindDescriptor(h_name, channel, dataId, descriptors);
        auto entry_ptr = std::make_shared<Entry>(h_name, this, desc);
        (*entry_ptr)().SetSystematicUncertainty(sample_unc.unc);
        (*entry_ptr)().SetPostfitScaleFactor(sample_unc.sf);
        histograms[h_name] = entry_ptr;
    }
}

EventAnalyzerData::Entry& EventAnalyzerData::GetHistogram(const std::string& hist_name) const
{
    const auto iter = histograms.find(hist_name);
    if(iter == histograms.end())
        throw exception("Histogram with name '%1%' not found.") % hist_name;
    return *iter->second;
}

const EventAnalyzerData::HistDesc& EventAnalyzerData::FindDescriptor(const std::string& h_name, Channel channel,
                                                                     const EventAnalyzerDataId& dataId,
                                                                     const HistDescCollection& descriptors)
{
    const std::vector<std::string> desc_name_candidates = {
        boost::str(boost::format("%1%/%2%/%3%/%4%") % h_name % channel % dataId.Get<EventCategory>()
                                                    % dataId.Get<EventSubCategory>()),
        boost::str(boost::format("%1%/%2%/*/%3%") % h_name % channel % dataId.Get<EventSubCategory>()),
        boost::str(boost::format("%1%/*/*/%2%") % h_name % dataId.Get<EventSubCategory>()),
        boost::str(boost::format("%1%/%2%/%3%") % h_name % channel % dataId.Get<EventCategory>()),
        boost::str(boost::format("%1%/*/%2%") % h_name % dataId.Get<EventCategory>()),
        boost::str(boost::format("%1%/%2%") % h_name % channel),
        h_name
    };
    for(const auto& desc_name : desc_name_candidates) {
        auto iter = descriptors.find(desc_name);
        if(iter != descriptors.end())
            return iter->second;
    }
    throw exception("Descriptor for histogram '%1%' not found.") % h_name;
}

} // namespace analysis
