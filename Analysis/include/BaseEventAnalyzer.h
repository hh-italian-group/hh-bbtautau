/*! Definition of BaseEventAnalyzer class, the base class for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "h-tautau/McCorrections/include/LeptonWeights.h"
#include "h-tautau/McCorrections/include/BTagWeight.h"
#include "hh-bbtautau/McCorrections/include/EventWeights_HH.h"
#include "h-tautau/Analysis/include/SignalObjectSelector.h"
#include "h-tautau/McCorrections/include/GenEventWeight.h"
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
    OPT_ARG(std::string, output_sync, "sync.root");
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
    void SetEventCategoryFlags(EventInfo& event, bbtautau::AnaTupleWriter::CategoriesFlags& category_flags) const;

    virtual EventRegion DetermineEventRegion(EventInfo& event, EventCategory eventCategory) = 0;

    void ProcessSamples(const std::vector<std::string>& sample_names, const std::string& sample_set_name);

    void ProcessDataSource(const SampleDescriptor& sample, const SampleDescriptor::Point& sample_wp,
                           std::shared_ptr<ntuple::EventTuple> tuple, const ntuple::ProdSummary& prod_summary);

    virtual void ProcessSpecialEvent(const SampleDescriptor& sample, const SampleDescriptor::Point& /*sample_wp*/,
                                     const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight,
                                     double shape_weight, bbtautau::AnaTupleWriter::DataIdMap& dataIds,
                                     double cross_section,
                                     std::map<UncertaintySource,std::map<UncertaintyScale,float>>& uncs_weight_map);

    bool SetRegionIsoRange(const LepCandidate& cand, EventRegion& region) const;

protected:
    AnalyzerArguments args;
    bbtautau::AnaTupleWriter anaTupleWriter;
    std::shared_ptr<TFile> outputFile_sync;
    std::vector<SyncDescriptor> sync_descriptors;
    std::map<std::string,std::shared_ptr<DYModelBase>> dymod;
    std::shared_ptr<NonResModel> nonResModel;
    const std::vector<std::string> trigger_patterns;
    std::shared_ptr<mc_corrections::EventWeights_HH> eventWeights_HH;
    bbtautau::AnaTupleWriter::BTagWeights btag_weights;
};

} // namespace analysis
