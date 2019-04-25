/*! Produce synchronization tree.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "h-tautau/Analysis/include/EventLoader.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/Btag_2017.h"
#include "h-tautau/McCorrections/include/EventWeights.h"
#include "hh-bbtautau/Analysis/include/SampleDescriptorConfigEntryReader.h"
#include "hh-bbtautau/Analysis/include/SyncTupleHTT.h"

struct Arguments {
    REQ_ARG(std::string, mode);
    REQ_ARG(std::string, input_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(std::string, period);
    REQ_ARG(std::string, output_file);
    OPT_ARG(std::string, mva_setup, "");
    OPT_ARG(bool, fill_tau_es_vars, false);
    OPT_ARG(bool, fill_jet_es_vars, false);
    OPT_ARG(std::string, jet_unc_source, "");
    OPT_ARG(std::string, jet_uncertainty, "");
};

namespace analysis {

enum class SyncMode { HTT, HH };

ENUM_NAMES(SyncMode) = {
    { SyncMode::HTT, "htt" },
    { SyncMode::HH, "hh" }
};

class SyncTreeProducer {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;
    using SyncEvent = htt_sync::SyncEvent;
    using SyncTuple = htt_sync::SyncTuple;

    static constexpr float default_value = std::numeric_limits<float>::lowest();
    static constexpr int default_int_value = std::numeric_limits<int>::lowest();

    SyncTreeProducer(const Arguments& _args) : args(_args), eventWeights(analysis::EnumNameMap<analysis::Period>::GetDefault().Parse(args.period()), JetOrdering::DeepCSV, DiscriminatorWP::Medium, true)
    {
        std::istringstream ss_mode(args.mode());
        ss_mode >> syncMode;
        run_period = analysis::EnumNameMap<analysis::Period>::GetDefault().Parse(args.period());

        if(args.mva_setup().size()) {
            ConfigReader config_reader;

            MvaReaderSetupCollection mva_setup_collection;
            MvaReaderSetupEntryReader mva_entry_reader(mva_setup_collection);
            config_reader.AddEntryReader("MVA", mva_entry_reader, true);
            config_reader.ReadConfig(args.mva_setup());

            std::vector<MvaReaderSetup> mva_setups;
            for(const auto& mva_setup_element : mva_setup_collection) {
                mva_setups.push_back(mva_setup_element.second);
            }
            mva_setup = mva_setups.size() == 1 ? mva_setups.front() : MvaReaderSetup::Join(mva_setups);

            mva_reader = std::make_shared<analysis::mva_study::MvaReader>();
            InitializeMvaReader();
        }


    }

    void Run()
    {
        std::cout << boost::format("Processing input file '%1%' into output file '%2%' using %3% mode.\n")
                   % args.input_file() % args.output_file() % args.mode();

        auto originalFile = root_ext::OpenRootFile(args.input_file());
        auto outputFile = root_ext::CreateRootFile(args.output_file());
        auto originalTuple = ntuple::CreateEventTuple(args.tree_name(),originalFile.get(),true,ntuple::TreeState::Full);
        SyncTuple sync(args.tree_name(), outputFile.get(), false);
        auto summaryTuple = ntuple::CreateSummaryTuple("summary", originalFile.get(), true, ntuple::TreeState::Full);
        summaryTuple->GetEntry(0);
        SummaryInfo summaryInfo(summaryTuple->data(), args.jet_unc_source());
        EventIdentifier current_id = EventIdentifier::Undef_event();
        std::map<EventEnergyScale, ntuple::Event> events;
        for(const auto& event : *originalTuple) {
            EventIdentifier event_id(event);
            if(event_id != current_id) {
                if(!events.empty()) {
                    FillSyncTuple(sync, events, summaryInfo);
                    events.clear();
                }
                current_id = event_id;
            }

            auto central_evt_iter = events.find(EventEnergyScale::Central);
            const ntuple::Event* central_evt = central_evt_iter != events.end() ? &central_evt_iter->second : nullptr;
            Event full_event = event;
            ntuple::EventLoader::Load(full_event, central_evt);
            const auto es = static_cast<EventEnergyScale>(event.eventEnergyScale);
            events[es] = full_event;
        }

        if(!events.empty()){
            FillSyncTuple(sync, events, summaryInfo);
        }

        sync.Write();

    }

private:

    void InitializeMvaReader()
    {
        using MvaKey = mva_study::MvaReader::MvaKey;
        if(!mva_setup.is_initialized())
            throw analysis::exception("Mva setup is not initialized.");
        for(const auto& method : mva_setup->trainings) {
            const auto& name = method.first;
            const auto& file = method.second;
            const auto& vars = mva_setup->variables.at(name);
            const auto& masses = mva_setup->masses.at(name);
            const auto& spins = mva_setup->spins.at(name);
            const bool legacy = mva_setup->legacy.count(name);
            const bool legacy_lm = legacy && mva_setup->legacy.at(name) == "lm";
            const size_t n_wp = masses.size();
            for(size_t n = 0; n < n_wp; ++n) {
                const MvaKey key{name, static_cast<int>(masses.at(n)), spins.at(n)};
                mva_reader->Add(key, file, vars, legacy, legacy_lm);
            }
        }
    }

    void FillSyncTuple(SyncTuple& sync, const std::map<EventEnergyScale, ntuple::Event>& events,const SummaryInfo& summaryInfo) const
    {
        static const std::map<Channel, std::vector<std::string>> triggerPaths = {
            { Channel::ETau, { "HLT_Ele32_WPTight_Gsf_v", "HLT_Ele35_WPTight_Gsf_v", "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v" } },
            { Channel::MuTau, { "HLT_IsoMu24_v", "HLT_IsoMu27_v", "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v" } },
            { Channel::TauTau, { "HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v",
                "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v",
                "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v"} },
            { Channel::MuMu, { "HLT_IsoMu24_v", "HLT_IsoMu27_v" } },
        };
        const Channel channel = Parse<Channel>(args.tree_name());
        std::map<EventEnergyScale, std::shared_ptr<EventInfoBase>> event_infos;
        for(const auto& entry : events) {
            const auto es = entry.first;
            const auto& event = entry.second;

            if(!args.fill_tau_es_vars() && (es == EventEnergyScale::TauUp || es == EventEnergyScale::TauDown)) continue;
            if((!args.fill_jet_es_vars() || !args.jet_uncertainty().empty())
                    && (es == EventEnergyScale::JetUp || es == EventEnergyScale::JetDown)) continue;
            if(syncMode == SyncMode::HH && (event.extraelec_veto || event.extramuon_veto)) continue;

            JetOrdering jet_ordering = run_period == Period::Run2017 ? JetOrdering::DeepCSV : JetOrdering::CSV;
            boost::optional<EventInfoBase> event_info_base = CreateEventInfo(event,&summaryInfo,analysis::TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017,run_period,jet_ordering);
            if(!event_info_base.is_initialized()) continue;
            if(syncMode == SyncMode::HH && !event_info_base->HasBjetPair()) continue;
            if(syncMode == SyncMode::HH && !event_info_base->GetTriggerResults().AnyAcceptAndMatch(triggerPaths.at(channel))) continue;
            if(syncMode == SyncMode::HTT && !event_info_base->GetTriggerResults().AnyAccept(triggerPaths.at(channel))) continue;


            /*
            static const std::vector<std::string> trigger_patterns = {
                "HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg_v"
            };
            analysis::EventInfoBase::JetCollection jets_vbf;
            analysis::EventInfoBase::JetPair vbf_jet_pair;
            jets_vbf = event_info->SelectJets(30, 5, std::numeric_limits<double>::lowest(),analysis::JetOrdering::Pt,
                                              event.GetSelectedBjetIndicesSet());
            vbf_jet_pair = event_info->SelectVBFJetPair(jets_vbf);
            if(vbf_jet_pair.first >= (*event_info)->jets_p4.size()
                    || vbf_jet_pair.second >= (*event_info)->jets_p4.size())
                continue;
            std::vector<ULong64_t> jet_trigger_match = {
                (*event_info)->jets_triggerFilterMatch.at(vbf_jet_pair.first),
                (*event_info)->jets_triggerFilterMatch.at(vbf_jet_pair.second)
            };
            if(!event_info->GetTriggerResults().AnyAcceptAndMatchEx(trigger_patterns, jet_trigger_match))
                continue;
            */
            std::shared_ptr<analysis::EventInfoBase> event_info = std::make_shared<analysis::EventInfoBase>(*event_info_base);
            event_infos[entry.first] = event_info;
        }

        if(!event_infos.count(EventEnergyScale::Central)) return;

        if(!args.jet_uncertainty().empty()) {
            event_infos[EventEnergyScale::JetUp] = event_infos[EventEnergyScale::Central]->ApplyShift(Parse<UncertaintySource>(args.jet_uncertainty()), UncertaintyScale::Up);
            event_infos[EventEnergyScale::JetDown] = event_infos[EventEnergyScale::Central]->ApplyShift(Parse<UncertaintySource>(args.jet_uncertainty()), UncertaintyScale::Down);
        }

        htt_sync::FillSyncTuple(*event_infos[EventEnergyScale::Central], sync, run_period, false, 1,
                                mva_reader.get(),
                                event_infos[EventEnergyScale::TauUp].get(),
                                event_infos[EventEnergyScale::TauDown].get(),
                                event_infos[EventEnergyScale::JetUp].get(),
                                event_infos[EventEnergyScale::JetDown].get());
    }

private:
    Arguments args;
    SyncMode syncMode;
    analysis::Period run_period;
    mc_corrections::EventWeights eventWeights;
    boost::optional<MvaReaderSetup> mva_setup;
    std::shared_ptr<analysis::mva_study::MvaReader> mva_reader;
};

} // namespace analysis

PROGRAM_MAIN(analysis::SyncTreeProducer, Arguments)
