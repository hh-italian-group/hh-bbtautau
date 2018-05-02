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
                      const SampleUnc& sample_unc) :
        AnalyzerData(outputFile, directoryName, anaTuple == nullptr)
    {
        for(const auto& h_name : histogram_names) {
            const auto& desc = FindDescriptor(h_name, channel, dataId, descriptors);
            auto entry_ptr = std::make_shared<Entry>(h_name, this, desc);
            const auto branch_ptr = anaTuple ? &anaTuple->get<ValueType>(h_name) : nullptr;
            (*entry_ptr)().SetSystematicUncertainty(sample_unc.unc);
            (*entry_ptr)().SetPostfitScaleFactor(sample_unc.sf);
            histograms[h_name] = EntrySource(entry_ptr, branch_ptr);
        }
    }

    void Fill(double weight)
    {
        if(ReadMode())
            throw exception("Fill is not available in read mode.");
        for(auto& entry : histograms)
            (*entry.second.first)().Fill(*entry.second.second, weight);
    }

private:
    static const HistDesc& FindDescriptor(const std::string& h_name, Channel channel, const EventAnalyzerDataId& dataId,
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

private:
    HistContainer histograms;
};

} // namespace analysis
