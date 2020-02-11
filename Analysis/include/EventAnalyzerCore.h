/*! Definition of core functionality for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Run/include/program_main.h"
#include "hh-bbtautau/Analysis/include/EventAnalyzerDataCollection.h"
#include "h-tautau/JetTools/include/BTagger.h"

namespace analysis {

struct CoreAnalyzerArguments {
    REQ_ARG(std::string, sources);
    REQ_ARG(std::string, setup);
    OPT_ARG(std::string, mva_sources, "");
    OPT_ARG(std::string, mva_setup, "");
    OPT_ARG(std::string, working_path, "./");
    OPT_ARG(unsigned, n_threads, 1);

    CoreAnalyzerArguments() {}
    CoreAnalyzerArguments(const CoreAnalyzerArguments&) = default;
    virtual ~CoreAnalyzerArguments() {}
};

class EventAnalyzerCore {
public:
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection;

    EventAnalyzerCore(const CoreAnalyzerArguments& args, Channel _channel);

    EventAnalyzerCore(const EventAnalyzerCore&) = delete;
    EventAnalyzerCore& operator=(const EventAnalyzerCore&) = delete;
    virtual ~EventAnalyzerCore() {}

    const std::string& ChannelNameLatex() const;
    std::string FullPath(const std::string& path) const;
    static bool FixNegativeContributions(TH1D& histogram, std::string& debug_info, std::string& negative_bins_info);
    void ProcessCombinedSamples(AnaDataCollection& anaDataCollection, const EventSubCategory& subCategory,
                                const std::vector<std::string>& sample_names);
    void AddSampleToCombined(AnaDataCollection& anaDataCollection, const EventSubCategory& subCategory,
                             CombinedSampleDescriptor& sample, SampleDescriptor& sub_sample);

private:
    void CreateEventSubCategoriesToProcess();
    void CreateMvaSelectionAliases();

    template<typename SampleCollection>
    std::vector<std::string> FilterInactiveSamples(const SampleCollection& samples,
                                                   const std::vector<std::string>& sample_names) const
    {
        std::vector<std::string> filtered;
        for(const auto& name : sample_names) {
            if(!samples.count(name))
                throw exception("Sample '%1%' not found.") % name;
            const auto& sample = samples.at(name);
            if(sample.channels.size() && !sample.channels.count(channelId)) continue;
            filtered.push_back(name);
        }
        return filtered;
    }

    void RemoveUnusedSamples();

    template<typename SampleCollection>
    void RemoveUnusedSamples(SampleCollection& samples,
                             const std::vector< std::vector<std::string>*> selected_sample_names) const
    {
        std::set<std::string> used_samples;
        for(auto selected_collection : selected_sample_names) {
            *selected_collection = FilterInactiveSamples(samples, *selected_collection);
            used_samples.insert(selected_collection->begin(), selected_collection->end());
        }

        const std::set<std::string> all_sample_names = tools::collect_map_keys(samples);
        for(const auto& sample_name : all_sample_names) {
            if(!used_samples.count(sample_name))
                samples.erase(sample_name);
        }
    }

protected:
    AnalyzerSetup ana_setup;
    boost::optional<MvaReaderSetup> mva_setup;
    SampleDescriptorCollection sample_descriptors;
    CombinedSampleDescriptorCollection cmb_sample_descriptors;
    EventSubCategorySet sub_categories_to_process;
    Channel channelId;
    std::map<SelectionCut, std::string> mva_sel_aliases;
    std::string working_path;
    std::shared_ptr<BTagger> bTagger;
    SignalObjectSelector signalObjectSelector;
};

} // namespace analysis
