/*! Definition of BaseEventAnalyzer class, the base class for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "h-tautau/McCorrections/include/LeptonWeights.h"
#include "h-tautau/McCorrections/include/TauIdWeight.h"
#include "h-tautau/McCorrections/include/BTagWeight.h"
#include "hh-bbtautau/McCorrections/include/EventWeights_HH.h"
#include "h-tautau/Analysis/include/SignalObjectSelector.h"
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

    SyncDescriptor(const std::string& desc_str, std::shared_ptr<TFile> outputFile_sync);

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
    std::vector<SyncDescriptor> sync_descriptors;
    std::map<std::string,std::shared_ptr<DYModelBase>> dymod;
    std::shared_ptr<NonResModel> nonResModel;
    const std::vector<std::string> trigger_patterns;
    std::shared_ptr<mc_corrections::EventWeights_HH> eventWeights_HH;
    SignalObjectSelector signalObjectSelector;
};

} // namespace analysis
