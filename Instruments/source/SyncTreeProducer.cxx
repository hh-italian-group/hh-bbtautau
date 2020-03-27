/*! Produce synchronization tree.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/McCorrections/include/EventWeights.h"
#include "hh-bbtautau/Analysis/include/SampleDescriptorConfigEntryReader.h"
#include "hh-bbtautau/Analysis/include/SyncTupleHTT.h"
#include "h-tautau/Analysis/include/SignalObjectSelector.h"
#include "AnalysisTools/Core/include/EventIdentifier.h"

struct Arguments {
    REQ_ARG(std::string, mode);
    REQ_ARG(std::string, input_file);
    REQ_ARG(std::string, channel);
    REQ_ARG(std::string, period);
    REQ_ARG(std::string, trigger_cfg);
    REQ_ARG(std::string, output_file);
    REQ_ARG(bool, apply_trigger_vbf);
    REQ_ARG(bool, isData);
    OPT_ARG(bool, use_svFit, false);
    OPT_ARG(std::string, tree_name, "");
    OPT_ARG(std::string, mva_setup, "");
    OPT_ARG(bool, fill_tau_es_vars, false);
    OPT_ARG(bool, fill_jet_es_vars, false);
    OPT_ARG(std::string, jet_unc_source, "");
    OPT_ARG(std::string, jet_uncertainty, "");
    OPT_ARG(std::string, event_id, "");
    OPT_ARG(bool, debug, false);
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
    // static constexpr int default_int_value = std::numeric_limits<int>::lowest();

    SyncTreeProducer(const Arguments& _args) : args(_args), syncMode(Parse<SyncMode>(args.mode())),
                                                            run_period(Parse<analysis::Period>(args.period())),
                                                            channel(Parse<Channel>(args.channel())),
                                                            signalObjectSelector(ConvertMode(syncMode))
                                                            // eventWeights(Parse<analysis::Period>(args.period()), JetOrdering::DeepCSV, DiscriminatorWP::Medium, true),

    {
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
        EventCandidate::InitializeUncertainties(run_period, false, ".",
                                                TauIdDiscriminator::byDeepTau2017v2p1VSjet);

    }

    void Run()
    {
        std::cout << boost::format("Processing input file '%1%' into output file '%2%' using %3% mode.\n")
                   % args.input_file() % args.output_file() % args.mode();

        std::map<std::string,std::pair<std::shared_ptr<ntuple::EventTuple>,Long64_t>> map_event;
        const std::string tree_name = args.tree_name().empty() ? args.channel() : args.tree_name();

        auto originalFile = root_ext::OpenRootFile(args.input_file());
        auto outputFile = root_ext::CreateRootFile(args.output_file());
        auto originalTuple = ntuple::CreateEventTuple(tree_name, originalFile.get(), true, ntuple::TreeState::Full);
        const Long64_t n_entries = originalTuple->GetEntries();

        SyncTuple sync(args.channel(), outputFile.get(), false);
        auto summaryTuple = ntuple::CreateSummaryTuple("summary", originalFile.get(), true, ntuple::TreeState::Full);
        summaryTuple->GetEntry(0);
        SummaryInfo summaryInfo(summaryTuple->data(), channel, args.trigger_cfg());
        std::cout << "n_entries " << n_entries << '\n';

        boost::optional<EventIdentifier> selected_event_id;
        if(!args.event_id().empty())
            selected_event_id = EventIdentifier(args.event_id());
        debug = selected_event_id.is_initialized() || args.debug();
        for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) {
            originalTuple->GetEntry(current_entry);
            ntuple::Event event = (*originalTuple).data();
            event.isData = args.isData();
            if(static_cast<Channel>(event.channelId) != channel) continue;

            if(debug) {
                const EventIdentifier eventId(event.run, event.lumi, event.evt);
                if(selected_event_id.is_initialized() && eventId != *selected_event_id) continue;
                std::cout << "Event: " << eventId << std::endl;
            }

            FillSyncTuple(sync, event, summaryInfo, debug);
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

    SignalMode ConvertMode(SyncMode syncMode)
    {
        std::map<SyncMode,SignalMode> signalMode_map ={
            {SyncMode::HTT, SignalMode::HTT},
            {SyncMode::HH, SignalMode::HH}
        };
        return signalMode_map.at(syncMode);
    }

    void FillSyncTuple(SyncTuple& sync, const ntuple::Event& event, const SummaryInfo& summaryInfo, bool debug) const
    {
        static const std::map<std::pair<Period, Channel>, std::vector<std::string>> triggerPaths = {
            { { Period::Run2016, Channel::ETau }, { "HLT_Ele25_eta2p1_WPTight_Gsf_v" } },
            { { Period::Run2016, Channel::MuTau }, { "HLT_IsoMu22_v", "HLT_IsoMu22_eta2p1_v", "HLT_IsoTkMu22_v",
                                "HLT_IsoTkMu22_eta2p1_v", "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v",
                                "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v" } },
            { { Period::Run2016, Channel::TauTau }, { "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v",
                                "HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v" } },
            { { Period::Run2016, Channel::MuMu }, { "HLT_IsoMu22_v" } },
            { { Period::Run2017, Channel::ETau }, { "HLT_Ele32_WPTight_Gsf_v", "HLT_Ele35_WPTight_Gsf_v",
                                "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v",
                                "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v" } },
            { { Period::Run2017, Channel::MuTau }, { "HLT_IsoMu24_v", "HLT_IsoMu27_v",
                                "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v" } },
            { { Period::Run2017, Channel::TauTau }, { "HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v",
                                "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v",
                                "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v" } },
            { { Period::Run2017, Channel::MuMu }, { "HLT_IsoMu24_v", "HLT_IsoMu27_v" } },
            { { Period::Run2018, Channel::ETau }, { "HLT_Ele32_WPTight_Gsf_v", "HLT_Ele35_WPTight_Gsf_v",
                                "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v",
                                "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v" } },
            { { Period::Run2018, Channel::MuTau }, { "HLT_IsoMu24_v", "HLT_IsoMu27_v",
                                "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v",
                                "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v" } },
            { { Period::Run2018, Channel::TauTau }, { "HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v",
                                 "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v",
                                 "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v",
                                 "HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v" } },
            { { Period::Run2018, Channel::MuMu }, { "HLT_IsoMu24_v", "HLT_IsoMu27_v" } },
        };

        static const std::map<std::pair<Period, Channel>, std::vector<std::string>> trigger_patterns_vbf = {
            { { Period::Run2017, Channel::TauTau }, {"HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg_v" } },
            { { Period::Run2018, Channel::TauTau }, { "HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_v",
                                 "HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_v" } },
        };
        static const JetOrdering jet_ordering = JetOrdering::DeepFlavour;
        const auto trig_key = std::make_pair(run_period, channel);

        if(debug)
            std::cout << "Creating event info..." << std::endl;
        auto event_info = CreateEventInfo(event, signalObjectSelector, &summaryInfo, run_period, jet_ordering, true,
                                          UncertaintySource::None, UncertaintyScale::Central, debug);
        if(debug)
            std::cout << "Event info is_initialized = " << event_info.is_initialized() << std::endl;
        if(!event_info.is_initialized()) return;
        if(debug) {
            for(size_t leg_id = 1; leg_id <= 2; ++leg_id) {
                const auto& leg = event_info->GetLeg(leg_id);
                std::cout << "Leg " << leg_id << ": type=" << leg->leg_type();
                if(!args.isData())
                    std::cout << ", gen_match=" << leg->gen_match();
                std::cout << ", charge=" << leg->charge() << ", dz=" << leg->dz() << "\n"
                          << "\t" << LorentzVectorToString(leg.GetMomentum());
                if(!args.isData())
                    std::cout << ", uncorrected " << LorentzVectorToString(leg->p4());
                std::cout << "\n";
                if(leg->leg_type() == LegType::tau) {
                    std::cout << "\tdecayMode=" << leg->decayMode() << "\n";
                }
            }
            if(event_info->HasBjetPair()) {
                for(size_t leg_id = 1; leg_id <= 2; ++leg_id) {
                    const auto& leg = event_info->GetBJet(leg_id);
                    std::cout << "b jet " << leg_id << ":\n"
                              << "\t" << LorentzVectorToString(leg.GetMomentum());
                    if(!args.isData())
                        std::cout << ", uncorrected " << LorentzVectorToString(leg->p4());
                    std::cout << "\n";
                }
            } else {
                std::cout << "No b-jet pair was selected.\n";
            }
            if(event_info->HasVBFjetPair()) {
                for(size_t leg_id = 1; leg_id <= 2; ++leg_id) {
                    const auto& leg = event_info->GetVBFJet(leg_id);
                    std::cout << "VBF jet " << leg_id << ":\n"
                              << "\t" << LorentzVectorToString(leg.GetMomentum());
                    if(!args.isData())
                        std::cout << ", uncorrected " << LorentzVectorToString(leg->p4());
                    std::cout << "\n";
                }
            } else {
                std::cout << "No VBF pair was selected.\n";
            }
            const auto& trig_list = triggerPaths.at(trig_key);
            for(const auto& trig : trig_list) {
                std::cout << trig << ": accept=" << event_info->GetTriggerResults().Accept(trig)
                          << " match=" << event_info->GetTriggerResults().MatchEx(trig,
                                                                        event_info->GetFirstLeg().GetMomentum().pt(),
                                                                        event_info->GetSecondLeg().GetMomentum().pt())
                          << std::endl;
            }
            if(args.apply_trigger_vbf() && trigger_patterns_vbf.count(trig_key)) {
                boost::optional<JetCandidate> first_vbf_jet, second_vbf_jet;
                boost::optional<std::vector<boost::multiprecision::uint256_t>> jet_trigger_match;
                if(event_info->HasVBFjetPair()) {
                    first_vbf_jet = event_info->GetVBFJet(1);
                    second_vbf_jet = event_info->GetVBFJet(2);
                    jet_trigger_match = std::vector<boost::multiprecision::uint256_t>({
                        (*first_vbf_jet)->triggerFilterMatch(), (*second_vbf_jet)->triggerFilterMatch(),
                    });
                }
                const auto& vbf_trig_list = trigger_patterns_vbf.at(trig_key);
                for(const auto& trig : vbf_trig_list) {
                    const bool match = event_info->HasVBFjetPair()
                                     ? event_info->GetTriggerResults().MatchEx(trig,
                                                           event_info->GetFirstLeg().GetMomentum().pt(),
                                                           event_info->GetSecondLeg().GetMomentum().pt(),
                                                           *jet_trigger_match)
                                     : false;
                    std::cout << trig << ": accept=" << event_info->GetTriggerResults().Accept(trig)
                              << " match=" << match << std::endl;
                }
            }
            std::cout << "PassLeptonVetoSelection = " << signalObjectSelector.PassLeptonVetoSelection(event) << "\n";
            std::cout << "Other leptons:\n";
            for(unsigned n = 0; n < event.other_lepton_p4.size(); ++n){
                if(static_cast<LegType>(event.other_lepton_type.at(n)) == LegType::e){
                    const DiscriminatorIdResults eleId_iso(event.other_lepton_eleId_iso.at(n));
                    const DiscriminatorIdResults eleId_noIso(event.other_lepton_eleId_noIso.at(n));
                    std::cout << "\tele " << LorentzVectorToString(event.other_lepton_p4.at(n))
                              << ", pass Medium iso id = " << eleId_iso.Passed(DiscriminatorWP::Medium)
                              << ", pass Medium noIso id = " << eleId_noIso.Passed(DiscriminatorWP::Medium)
                              << ", pfRelIso04 = " << event.other_lepton_iso.at(n) << "\n";
                }
                if(static_cast<LegType>(event.other_lepton_type.at(n)) == LegType::mu){
                    analysis::DiscriminatorIdResults muonId(event.other_lepton_muonId.at(n));
                    std::cout << "\tmuon " << LorentzVectorToString(event.other_lepton_p4.at(n))
                              << ", pass Medium id = " << muonId.Passed(DiscriminatorWP::Medium)
                              << ", pass Tight id = " << muonId.Passed(DiscriminatorWP::Tight)
                              << ", pfRelIso04 = " << event.other_lepton_iso.at(n) << "\n";
                }
            }
            std::cout << "PassMETfilters = " << signalObjectSelector.PassMETfilters(event,run_period,args.isData())
                      << "\n";
            std::cout << "MET " << LorentzVectorToString(event_info->GetMET().GetMomentum(), LVectorRepr::PxPyPtPhi);
            if(!args.isData())
                std::cout << ", uncorrected " << LorentzVectorToString(event_info->GetMET()->p4(),
                                                                       LVectorRepr::PxPyPtPhi);
            std::cout << "\n";
            if(!args.isData()) {
                auto& evt_cand = event_info->GetEventCandidate();
                std::cout << "Lepton corrections:\n";
                LorentzVectorXYZ total_delta(0, 0, 0, 0);
                for(const auto& lep : evt_cand.GetLeptons()) {
                    const auto& delta = lep.GetMomentum() - lep->p4();
                    total_delta += delta;
                    std::cout << "\t" << lep->leg_type() << ", ";
                    if(lep->leg_type() == LegType::tau) {
                        std::cout << "decayMode=" << lep->decayMode() << ", gen_match=" << lep->gen_match() << ", ";
                    }
                    std::cout << LorentzVectorToString(lep->p4(), LVectorRepr::PxPyPzE)
                              << " -> " << LorentzVectorToString(lep.GetMomentum(), LVectorRepr::PxPyPzE, false)
                              << ", shift = " << LorentzVectorToString(delta, LVectorRepr::PxPyPzE, false) << "\n";
                }
                std::cout << "Total lepton shift: " << LorentzVectorToString(total_delta, LVectorRepr::PxPyPzE) << "\n";
            }

            const auto& sv_fit = event_info->GetSVFitResults(true, 1);
            std::cout << "SVfit: " << LorentzVectorToString(sv_fit.momentum, LVectorRepr::PtEtaPhiM) << "\n";
            const auto& kin_fit = event_info->GetKinFitResults(true, 1);
            std::cout << "KinFit: convergence=" << kin_fit.convergence << ", mass=" << kin_fit.mass
                      << ", chi2=" << kin_fit.chi2 << "\n";
        }
        bool pass_trigger = event_info->GetTriggerResults().AnyAcceptAndMatchEx(triggerPaths.at(trig_key),
                event_info->GetFirstLeg().GetMomentum().pt(), event_info->GetSecondLeg().GetMomentum().pt());
        if(!pass_trigger && args.apply_trigger_vbf() && trigger_patterns_vbf.count(trig_key)
                && event_info->HasVBFjetPair()) {
            const auto& first_vbf_jet = event_info->GetVBFJet(1);
            const auto& second_vbf_jet = event_info->GetVBFJet(2);

            std::vector<boost::multiprecision::uint256_t> jet_trigger_match = {
                first_vbf_jet->triggerFilterMatch(),
                second_vbf_jet->triggerFilterMatch()
            };
            pass_trigger = event_info->GetTriggerResults().AnyAcceptAndMatchEx(trigger_patterns_vbf.at(trig_key),
                    event_info->GetFirstLeg().GetMomentum().pt(), event_info->GetSecondLeg().GetMomentum().pt(),
                    jet_trigger_match);
        }
        if(!pass_trigger) return;
        if(syncMode == SyncMode::HH && !event_info->HasBjetPair()) return;
        if(syncMode == SyncMode::HH && !signalObjectSelector.PassLeptonVetoSelection(event)) return;
        if(syncMode == SyncMode::HH && !signalObjectSelector.PassMETfilters(event,run_period,args.isData())) return;
        htt_sync::FillSyncTuple(*event_info, sync, run_period, args.use_svFit(), 1,
                                mva_reader.get(), nullptr, nullptr, nullptr, nullptr);
    }

private:
    Arguments args;
    SyncMode syncMode;
    analysis::Period run_period;
    Channel channel;
    // mc_corrections::EventWeights eventWeights;
    boost::optional<MvaReaderSetup> mva_setup;
    std::shared_ptr<analysis::mva_study::MvaReader> mva_reader;
    SignalObjectSelector signalObjectSelector;
    bool debug{false};
};

} // namespace analysis

PROGRAM_MAIN(analysis::SyncTreeProducer, Arguments)
