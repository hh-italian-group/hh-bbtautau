/*! Definition of BaseEventAnalyzer class, the base class for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "DYModel.h"
#include "EventAnalyzerCore.h"
#include "MvaReader.h"
#include "NonResModel.h"
#include "SyncTupleHTT.h"
#include <boost/regex.hpp>

namespace analysis {

struct AnalyzerArguments : CoreAnalyzerArguments {
    REQ_ARG(std::string, input);
    REQ_ARG(std::string, output);
    OPT_ARG(std::string, output_sync,"sync.root");
};

struct SyncDescriptor {
    std::shared_ptr<htt_sync::SyncTuple> sync_tree;
    std::shared_ptr<boost::regex> regex_pattern;

    SyncDescriptor(std::shared_ptr<htt_sync::SyncTuple> _sync_tree,
                   std::shared_ptr<boost::regex> _regex_pattern) :
    sync_tree(_sync_tree), regex_pattern(_regex_pattern) {}

};

class BaseEventAnalyzer : public EventAnalyzerCore {
public:
    using Event = ntuple::Event;

    BaseEventAnalyzer(const AnalyzerArguments& _args, Channel channel);
    void Run();

protected:
    EventCategorySet DetermineEventCategories(EventInfoBase& event);
    virtual EventRegion DetermineEventRegion(EventInfoBase& event, EventCategory eventCategory) = 0;

    void InitializeMvaReader();
    virtual EventSubCategory DetermineEventSubCategory(EventInfoBase& event, const EventCategory& category,
                                                       std::map<SelectionCut, double>& mva_scores);
    void ProcessSamples(const std::vector<std::string>& sample_names, const std::string& sample_set_name);

    void ProcessDataSource(const SampleDescriptor& sample, const SampleDescriptor::Point& sample_wp,
                           std::shared_ptr<ntuple::EventTuple> tuple, const ntuple::ProdSummary& prod_summary);

    virtual void ProcessSpecialEvent(const SampleDescriptor& sample, const SampleDescriptor::Point& /*sample_wp*/,
                                     const EventAnalyzerDataId& anaDataId, EventInfoBase& event, double weight,
                                     double shape_weight, bbtautau::AnaTupleWriter::DataIdMap& dataIds);

protected:
    AnalyzerArguments args;
    bbtautau::AnaTupleWriter anaTupleWriter;
    mva_study::MvaReader mva_reader;
    std::shared_ptr<TFile> outputFile_sync;
    // std::map<EventAnalyzerDataId, std::shared_ptr<htt_sync::SyncTuple>> syncTuple_map;
    std::vector<SyncDescriptor> sync_descriptors;
    std::map<std::string,std::shared_ptr<DYModel>> dymod;
    std::shared_ptr<NonResModel> nonResModel;
    const std::vector<std::string> trigger_patterns;
};

} // namespace analysis
