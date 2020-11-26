/*! Definition of BaseEventAnalyzer class, the base class for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */


#include "hh-bbtautau/Analysis/include/BaseEventAnalyzer.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "AnalysisTools/Run/include/MultiThread.h"
#include "h-tautau/McCorrections/include/JetPuIdWeights.h"
#include "h-tautau/McCorrections/include/TopPtWeight.h"

namespace analysis {

SyncDescriptor::SyncDescriptor(const std::string& desc_str, std::shared_ptr<TFile> outputFile_sync)
{
    auto tree_regexes = SplitValueList(desc_str, false, ":");
    if(tree_regexes.size() != 2)
        throw exception("The Number of parameters is %1%, only 2 are allowed") %  tree_regexes.size();
    sync_tree = std::make_shared<htt_sync::SyncTuple>(tree_regexes.at(0),outputFile_sync.get(),false);
    regex_pattern = std::make_shared<boost::regex>(tree_regexes.at(1));
}

BaseEventAnalyzer::BaseEventAnalyzer(const AnalyzerArguments& _args, Channel channel) :
    EventAnalyzerCore(_args, channel, true), args(_args), anaTupleWriter(args.output(), channel, ana_setup.use_kinFit,
                      ana_setup.use_svFit,ana_setup.allow_calc_svFit), trigger_patterns(ana_setup.trigger.at(channel))

{
    std::cout << "Initializing uncertainties... " << std::flush;
    EventCandidate::InitializeUncertainties(ana_setup.period, false, args.working_path(),
                                            signalObjectSelector.GetTauVSjetDiscriminator().first);
    std::cout << "done." << std::endl;
    if(ana_setup.syncDataIds.size()) {
        outputFile_sync = root_ext::CreateRootFile(args.output_sync());
        for(unsigned n = 0; n < ana_setup.syncDataIds.size(); ++n){
            sync_descriptors.emplace_back(ana_setup.syncDataIds.at(n),outputFile_sync);
        }
    }
    for(auto& x: sample_descriptors) {
        SampleDescriptor& sample = x.second;
        if(sample.sampleType == SampleType::DY){
            if(sample.fit_method == DYFitModel::Htt)
                dymod[sample.name] = std::make_shared<DYModel_HTT>(sample, args.working_path());
            else
                dymod[sample.name] = std::make_shared<DYModel>(sample, args.working_path());
        }
    }
    std::cout << "Creating event weights... " << std::endl;
    eventWeights_HH = std::make_shared<mc_corrections::EventWeights_HH>(ana_setup.period, *bTagger);
}

void BaseEventAnalyzer::Run()
{
    run::ThreadPull threads(args.n_threads());
    ProcessSamples(ana_setup.signals, "signal");
    ProcessSamples(ana_setup.data, "data");
    ProcessSamples(ana_setup.backgrounds, "background");
    std::cout << "Saving output file..." << std::endl;
    for (size_t n = 0; n < sync_descriptors.size(); ++n) {
        auto sync_tree = sync_descriptors.at(n).sync_tree;
        sync_tree->Write();
    }
}

void BaseEventAnalyzer::SetEventCategoryFlags(EventInfo& event,
        bbtautau::AnaTupleWriter::CategoriesFlags& category_flags) const
{

    static const std::map<DiscriminatorWP, size_t> btag_working_points = {{DiscriminatorWP::Loose, 0},
                                                                          {DiscriminatorWP::Medium, 0},
                                                                          {DiscriminatorWP::Tight, 0}};
    // const bool is_boosted =  false;
    category_flags.fat_jet_cand = SignalObjectSelector::SelectFatJet(event.GetEventCandidate(),
                                                                     event.GetSelectedSignalJets());
    category_flags.is_boosted = category_flags.fat_jet_cand != nullptr;
    category_flags.is_vbf = false;
    if(event.HasVBFjetPair()) {
        const auto vbf_jet_1 = event.GetVBFJet(1).GetMomentum();
        const auto vbf_jet_2 = event.GetVBFJet(2).GetMomentum();
        const auto deta_jj =  std::abs(vbf_jet_1.Eta() - vbf_jet_2.Eta());
        const auto m_jj = (vbf_jet_1 + vbf_jet_2).M();

        category_flags.is_vbf = (m_jj > cuts::hh_bbtautau_Run2::VBF::mass_jj
                                && deta_jj > cuts::hh_bbtautau_Run2::VBF::deltaeta_jj);
    }

    const auto& all_jets = event.GetCentralJets();
    category_flags.num_jets = all_jets.size();
    category_flags.num_btag = btag_working_points;
    category_flags.nbtag = btag_working_points;
    if(event.HasBjetPair()) {
        for(size_t bjet_index = 1; bjet_index <= 2; ++bjet_index) {
            const auto& jet = event.GetBJet(bjet_index);
            for(const auto& btag_wp : btag_working_points) {
                if(bTagger->Pass(*jet, btag_wp.first))
                    ++category_flags.num_btag[btag_wp.first];
            }
        }
    }
    for(const auto& jet : all_jets) {
        for(const auto& btag_wp : btag_working_points) {
            if(bTagger->Pass(**jet, btag_wp.first))
                ++category_flags.nbtag[btag_wp.first];
        }
    }
}

void BaseEventAnalyzer::ProcessSamples(const std::vector<std::string>& sample_names, const std::string& sample_set_name)
{
    std::cout << "Processing " << sample_set_name << " samples... " << std::endl;
    for(size_t sample_index = 0; sample_index < sample_names.size(); ++sample_index) {
        const std::string& sample_name = sample_names.at(sample_index);
        if(!sample_descriptors.count(sample_name))
            throw exception("Sample '%1%' not found while processing.") % sample_name;
        SampleDescriptor& sample = sample_descriptors.at(sample_name);
        if(sample.sampleType == SampleType::QCD || (sample.channels.size() && !sample.channels.count(channelId)))
            continue;
        std::cout << '\t' << sample.name << std::endl;

        std::set<std::string> processed_files;
        for(const auto& sample_wp : sample.working_points) {
            if(!sample_wp.file_path.size() || processed_files.count(sample_wp.file_path)) continue;
            auto file = root_ext::OpenRootFile(tools::FullPath({args.input(), sample_wp.file_path}));
            auto tuple = ntuple::CreateEventTuple(ToString(channelId), file.get(), true,
                                                  ntuple::TreeState::Skimmed);
            auto summary_tuple = ntuple::CreateSummaryTuple("summary", file.get(), true,
                                                            ntuple::TreeState::Skimmed);
            const auto prod_summary = ntuple::MergeSummaryTuple(*summary_tuple);
            if(sample.sampleType == SampleType::ggHH_NonRes) {
                std::cout << "\t\tpreparing NonResModel... ";
                std::cout.flush();
                if(!crossSectionProvider)
                    throw exception("path to the cross section config should be specified.");
                nonResModel = std::make_shared<NonResModel>(*eventWeights_HH, sample, file, *crossSectionProvider);
                std::cout << "done." << std::endl;
            }
            ProcessDataSource(sample, sample_wp, tuple, prod_summary);

            processed_files.insert(sample_wp.file_path);
        }
    }
}

void BaseEventAnalyzer::ProcessDataSource(const SampleDescriptor& sample, const SampleDescriptor::Point& sample_wp,
                                          std::shared_ptr<ntuple::EventTuple> tuple,
                                          const ntuple::ProdSummary& prod_summary)
{
    std::vector<std::string> vbf_triggers;
    std::map<UncertaintySource, std::map<UncertaintyScale, float>>  uncs_weight_map;
    std::vector<DiscriminatorWP> btag_wps = { DiscriminatorWP::Loose, DiscriminatorWP::Medium,
                                              DiscriminatorWP::Tight };

    if(ana_setup.trigger_vbf.count(channelId))
        vbf_triggers = ana_setup.trigger_vbf.at(channelId);
    const auto summary = std::make_shared<SummaryInfo>(prod_summary, channelId, ana_setup.trigger_path,
                                                       trigger_patterns, vbf_triggers);
    std::set<UncertaintySource> unc_sources = { UncertaintySource::None };
    if(sample.sampleType != SampleType::Data)
        unc_sources = ana_setup.unc_sources;
    const auto unc_variations = EnumerateUncVariations(unc_sources);
    const bool is_data = sample.sampleType == SampleType::Data;
    for(const auto& tupleEvent : *tuple) {
        if(!signalObjectSelector.PassMETfilters(tupleEvent,ana_setup.period,is_data)) continue;
        if(!signalObjectSelector.PassLeptonVetoSelection(tupleEvent)) continue;

        for(const auto& [unc_source, unc_scale] : unc_variations) {
            auto event = EventInfo::Create(tupleEvent, signalObjectSelector, *bTagger, DiscriminatorWP::Medium,
                                           summary, unc_source, unc_scale);
            if(!event || !event->HasBjetPair()) continue;
            const bool pass_normal_trigger = event->PassNormalTriggers();
            const bool pass_vbf_trigger = event->PassVbfTriggers();
            const bool pass_trigger = pass_normal_trigger || pass_vbf_trigger;
            const bool pass_only_vbf_trigger = !pass_normal_trigger && pass_vbf_trigger;
            if(!pass_trigger) continue;
            // std::map<std::pair<DiscriminatorWP, bool>, std::map<UncertaintyScale, float>> btag_weights;
            bbtautau::AnaTupleWriter::DataIdMap dataIds;
            std::map<size_t, bool> sync_event_selected;

            double lepton_id_iso_weight = 1., trigger_weight = 1., prescale_weight = 1., l1_prefiring_weight = 1.,
                   jet_pu_id_weight = 1,  btag_weight = 1.;

            if(sample.sampleType != SampleType::Data) {
                auto lepton_weight_provider = eventWeights_HH->GetProviderT<mc_corrections::LeptonWeights>(
                        mc_corrections::WeightType::LeptonTrigIdIso);
                lepton_id_iso_weight = lepton_weight_provider->GetIdIsoWeight(*event,
                        signalObjectSelector.GetTauVSeDiscriminator(channelId).second,
                        signalObjectSelector.GetTauVSmuDiscriminator(channelId).second,
                        signalObjectSelector.GetTauVSjetDiscriminator().second,
                        unc_source, unc_scale);
                trigger_weight = lepton_weight_provider->GetTriggerWeight(*event,
                        signalObjectSelector.GetTauVSjetDiscriminator().second, unc_source, unc_scale);

                prescale_weight = lepton_weight_provider->GetTriggerPrescaleWeight(*event);

                if(ana_setup.period == Period::Run2016 || ana_setup.period == Period::Run2017)
                    l1_prefiring_weight = (*event)->l1_prefiring_weight;

                uncs_weight_map[UncertaintySource::L1_prefiring][UncertaintyScale::Central] =
                        ana_setup.period == Period::Run2016 || ana_setup.period == Period::Run2017
                        ? static_cast<float>((*event)->l1_prefiring_weight) : 1;
                uncs_weight_map[UncertaintySource::L1_prefiring][UncertaintyScale::Up] =
                        ana_setup.period == Period::Run2016 || ana_setup.period == Period::Run2017
                        ? static_cast<float>((*event)->l1_prefiring_weight_up) : 1;
                uncs_weight_map[UncertaintySource::L1_prefiring][UncertaintyScale::Down] =
                        ana_setup.period == Period::Run2016 || ana_setup.period == Period::Run2017
                        ? static_cast<float>((*event)->l1_prefiring_weight_down) : 1;

                auto jet_pu_id_weight_provided = eventWeights_HH->GetProviderT<mc_corrections::JetPuIdWeights>(
                        mc_corrections::WeightType::JetPuIdWeights);
                jet_pu_id_weight = jet_pu_id_weight_provided->GetWeight(*event);

                auto top_pt_weight_provided = eventWeights_HH->GetProviderT<mc_corrections::TopPtWeight>(
                        mc_corrections::WeightType::TopPt);

                auto btag_weight_provider = eventWeights_HH->GetProviderT<mc_corrections::BTagWeight>(
                    mc_corrections::WeightType::BTag);


                    static const std::vector<UncertaintySource> btag_sources = {
                        UncertaintySource::btag_lf, UncertaintySource::btag_hf, UncertaintySource::btag_hfstats1,
                        UncertaintySource::btag_hfstats2, UncertaintySource::btag_lfstats1,
                        UncertaintySource::btag_lfstats2, UncertaintySource::btag_cferr1, UncertaintySource::btag_cferr2
                     };


                for(const bool iter_fit : {false, true}) {
                    for(const auto wp : btag_wps){
                        const std::pair<DiscriminatorWP, bool> key(wp, iter_fit);
                        btag_weights.weights[wp][UncertaintyScale::Central] = static_cast<float>(
                                btag_weight_provider->Get(*event, wp, iter_fit, unc_source, unc_scale, false));
                        if(unc_source == UncertaintySource::None && !iter_fit) {
                            for(const auto scale : {UncertaintyScale::Up, UncertaintyScale::Down})
                                btag_weights.weights[wp][scale] = static_cast<float>(btag_weight_provider->Get(*event,
                                    wp, iter_fit, UncertaintySource::Eff_b, scale, false));
                        }
                        else if(unc_source == UncertaintySource::None && iter_fit) {
                            for(const auto scale : {UncertaintyScale::Up, UncertaintyScale::Down}) {
                                for(const bool apply_JES : {false, true}) {
                                    if(apply_JES) {
                                        for(const auto unc_jes : {UncertaintySource::JetFull_Total,
                                                                  UncertaintySource::JetReduced_Total})
                                            btag_weights.weights_iter_jes[wp][scale] = static_cast<float>(
                                                btag_weight_provider->Get(*event, wp, iter_fit, unc_jes, scale, true));
                                    } else {
                                        for(UncertaintySource unc : btag_sources)
                                            btag_weights.weights_iter[wp][scale] = static_cast<float>(
                                            btag_weight_provider->Get(*event, wp, iter_fit, unc, scale, false));
                                    }
                                }
                            }
                        }
                    }
                }

                if(unc_source == UncertaintySource::None) {
                    static const std::vector<UncertaintySource> uncs_trigger_weight = {
                        UncertaintySource::EleTriggerUnc, UncertaintySource::MuonTriggerUnc,
                        UncertaintySource::TauTriggerUnc_DM0, UncertaintySource::TauTriggerUnc_DM1,
                        UncertaintySource::TauTriggerUnc_DM10, UncertaintySource::TauTriggerUnc_DM11,
                        UncertaintySource::VBFTriggerUnc_jets, UncertaintySource::VBFTauTriggerUnc_DM0,
                        UncertaintySource::VBFTauTriggerUnc_DM1, UncertaintySource::VBFTauTriggerUnc_3prong
                    };
                    static const std::vector<UncertaintySource> uncs_tau_weight = { UncertaintySource::TauVSjetSF_DM0,
                        UncertaintySource::TauVSjetSF_DM1, UncertaintySource::TauVSjetSF_3prong,
                        UncertaintySource::TauVSjetSF_pt20to25, UncertaintySource::TauVSjetSF_pt25to30,
                        UncertaintySource::TauVSjetSF_pt30to35, UncertaintySource::TauVSjetSF_pt35to40,
                        UncertaintySource::TauVSjetSF_ptgt40, UncertaintySource::TauVSeSF_barrel,
                        UncertaintySource::TauVSeSF_endcap, UncertaintySource::TauVSmuSF_etaLt0p4,
                        UncertaintySource::TauVSmuSF_eta0p4to0p8, UncertaintySource::TauVSmuSF_eta0p8to1p2,
                        UncertaintySource::TauVSmuSF_eta1p2to1p7, UncertaintySource::TauVSmuSF_etaGt1p7,
                        UncertaintySource::EleIdIsoUnc, UncertaintySource::MuonIdIsoUnc,
                        UncertaintySource::TauCustomSF_DM0, UncertaintySource::TauCustomSF_DM1,
                        UncertaintySource::TauCustomSF_DM10, UncertaintySource::TauCustomSF_DM11 };

                    for (UncertaintySource unc : uncs_trigger_weight) {
                        for(UncertaintyScale unc_scale_eval : GetAllUncertaintyScales()) {
                            const UncertaintySource unc_eval = unc_scale_eval == UncertaintyScale::Central ?
                                                               UncertaintySource::None : unc;
                            uncs_weight_map[unc][unc_scale_eval] = static_cast<float>(
                                lepton_weight_provider->GetTriggerWeight(*event,
                                    signalObjectSelector.GetTauVSjetDiscriminator().second, unc_eval, unc_scale_eval));
                        }
                    }
                    for (UncertaintySource unc : uncs_tau_weight) {
                        for(UncertaintyScale unc_scale_eval : GetAllUncertaintyScales()) {
                            const UncertaintySource unc_eval = unc_scale_eval == UncertaintyScale::Central ?
                                                               UncertaintySource::None : unc;
                            uncs_weight_map[unc][unc_scale_eval] = static_cast<float>(
                                lepton_weight_provider->GetIdIsoWeight(*event,
                                signalObjectSelector.GetTauVSeDiscriminator(channelId).second,
                                signalObjectSelector.GetTauVSmuDiscriminator(channelId).second,
                                signalObjectSelector.GetTauVSjetDiscriminator().second, unc_eval, unc_scale_eval));
                        }
                    }
                    uncs_weight_map[UncertaintySource::TopPt][UncertaintyScale::Central] =
                            static_cast<float>(1. / (*summary)->totalShapeWeight);
                    if(sample.apply_top_pt_unc)
                        uncs_weight_map[UncertaintySource::TopPt][UncertaintyScale::Up] = static_cast<float>(
                               top_pt_weight_provided->Get(*event) / (*summary)->totalShapeWeight_withTopPt);

                    if(sample.sampleType != SampleType::ggHH_NonRes){
                        uncs_weight_map[UncertaintySource::PileUp][UncertaintyScale::Central] = static_cast<float>(
                               tupleEvent.weight_pu / (*summary)->totalShapeWeight);
                        uncs_weight_map[UncertaintySource::PileUp][UncertaintyScale::Up] = static_cast<float>(
                               tupleEvent.weight_pu_up / (*summary)->totalShapeWeight_withPileUp_Up);
                        uncs_weight_map[UncertaintySource::PileUp][UncertaintyScale::Down] = static_cast<float>(
                               tupleEvent.weight_pu_down / (*summary)->totalShapeWeight_withPileUp_Down);
                    }

                    for(auto unc_jet_pu_id : {UncertaintySource::PileUpJetId_eff, UncertaintySource::PileUpJetId_mistag}){
                        for(UncertaintyScale unc_scale_eval : GetAllUncertaintyScales()) {
                                const UncertaintySource unc_eval = unc_scale_eval == UncertaintyScale::Central ?
                                                                   UncertaintySource::None : unc_jet_pu_id;
                                uncs_weight_map[unc_jet_pu_id][unc_scale_eval] = static_cast<float>(
                                    jet_pu_id_weight_provided->GetWeight(*event, unc_eval, unc_scale_eval));
                        }
                    }
                }
            }

            bbtautau::AnaTupleWriter::CategoriesFlags category_flags;
            category_flags.pass_vbf_trigger = pass_vbf_trigger;
            category_flags.pass_only_vbf_trigger = pass_only_vbf_trigger;
            SetEventCategoryFlags(*event, category_flags);
            static const EventCategory eventCategory = EventCategory::Parse("2j");
            const EventRegion eventRegion = DetermineEventRegion(*event, eventCategory);
            for(const auto& region : ana_setup.regions) {
                if(!eventRegion.Implies(region)) continue;
                static const EventSubCategory subCategory = EventSubCategory::NoCuts();
                const EventAnalyzerDataId anaDataId(eventCategory, subCategory, region,
                                                    unc_source, unc_scale, sample_wp.full_name);
                double shape_weight = 1.;
                if(sample.sampleType == SampleType::Data) {
                    dataIds[anaDataId] = 1. ;
                } else {

                    if(eventCategory.HasBtagConstraint())
                         btag_weight = btag_weights.weights.at(eventCategory.BtagWP()).at(UncertaintyScale::Central);

                    double cross_section = (*summary)->cross_section > 0 ? (*summary)->cross_section :
                                                                            sample.cross_section;

                    auto gen_weight_provider = eventWeights_HH->GetProviderT<mc_corrections::GenEventWeight>(
                            mc_corrections::WeightType::GenEventWeight);
                    shape_weight = cross_section * gen_weight_provider->Get(*event);
                    const double weight = (*event)->weight_total * cross_section *  ana_setup.int_lumi
                                           * lepton_id_iso_weight * trigger_weight * prescale_weight
                                           * l1_prefiring_weight * btag_weight * jet_pu_id_weight
                                           / (*summary)->totalShapeWeight;
                    if(sample.sampleType == SampleType::MC) {
                        dataIds[anaDataId] = weight;
                    } else
                        ProcessSpecialEvent(sample, sample_wp, anaDataId, *event, weight,
                                            (*summary)->totalShapeWeight, dataIds, cross_section, uncs_weight_map);
                }

                for(size_t n = 0; n < sync_descriptors.size(); ++n) {
                    if(sync_event_selected[n]) continue;
                    const auto& regex_pattern = sync_descriptors.at(n).regex_pattern;
                    for(auto& dataId : dataIds) {
                        if(!boost::regex_match(dataId.first.GetName(), *regex_pattern)) continue;

                        htt_sync::FillSyncTuple(*event, *sync_descriptors.at(n).sync_tree, ana_setup.period,
                                                ana_setup.use_svFit, dataId.second,
                                                lepton_id_iso_weight, trigger_weight, btag_weight,
                                                shape_weight, jet_pu_id_weight, prescale_weight);
                        sync_event_selected[n] = true;
                        break;
                    }
                }
            }
            //dataId
            anaTupleWriter.AddEvent(*event, dataIds, category_flags, btag_weights, uncs_weight_map);
        }
    }
}

void BaseEventAnalyzer::ProcessSpecialEvent(const SampleDescriptor& sample,
                                            const SampleDescriptor::Point& /*sample_wp*/,
                                            const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight,
                                            double shape_weight, bbtautau::AnaTupleWriter::DataIdMap& dataIds,
                                            double cross_section,
                                            std::map<UncertaintySource,std::map<UncertaintyScale,float>>& uncs_weight_map)
{
    if(sample.sampleType == SampleType::DY){
        dymod.at(sample.name)->ProcessEvent(anaDataId,event,weight,dataIds);
    } else if(sample.sampleType == SampleType::ggHH_NonRes) {
        nonResModel->ProcessEvent(anaDataId, event, weight, shape_weight, dataIds, cross_section, uncs_weight_map);
    } else
        throw exception("Unsupported special event type '%1%'.") % sample.sampleType;
}

bool BaseEventAnalyzer::SetRegionIsoRange(const LepCandidate& cand, EventRegion& region) const
{
    if(cand->leg_type() == LegType::tau) {
        const auto& [tau_discr, wp_max] = signalObjectSelector.GetTauVSjetDiscriminator();
        const auto wp_min = signalObjectSelector.GetTauVSjetSidebandWPRange().first;
        for(int wp_index = static_cast<int>(wp_max); wp_index >= static_cast<int>(wp_min); --wp_index) {
            const DiscriminatorWP wp = static_cast<DiscriminatorWP>(wp_index);
            if(cand->Passed(tau_discr, wp)) {
                region.SetLowerIso(wp);
                if(wp != wp_max){
                    region.SetUpperIso(static_cast<DiscriminatorWP>(wp_index + 1));}
                break;
            }
        }
    } else if(cand->leg_type() == LegType::mu) {
        static const std::map<DiscriminatorWP, double> working_points = {
            { DiscriminatorWP::VVLoose, 2.0 },
            { DiscriminatorWP::Loose, 0.3 },
            { DiscriminatorWP::Medium, ::cuts::hh_bbtautau_Run2::MuTau::muonID::pfRelIso04 }
        };
        if(cand.GetIsolation() > working_points.at(DiscriminatorWP::Medium)
                && cand.GetIsolation() < working_points.at(DiscriminatorWP::Loose)) return false;
        for(auto wp = working_points.rbegin(); wp != working_points.rend(); ++wp) {
            if(cand.GetIsolation() < wp->second) {
                region.SetLowerIso(wp->first);
                if(wp != working_points.rbegin())
                    region.SetUpperIso((--wp)->first);
                break;
            }
        }
    } else if(cand->leg_type() == LegType::e) {
        static const std::map<DiscriminatorWP, double> working_points = {
            { DiscriminatorWP::VVLoose, 2.0 },
            { DiscriminatorWP::Loose, 0.3 },
            { DiscriminatorWP::Medium, ::cuts::H_tautau_Run2::ETau::electronID::pfRelIso04 }
        };
        if(cand.GetIsolation() > working_points.at(DiscriminatorWP::Medium)
                && cand.GetIsolation() < working_points.at(DiscriminatorWP::Loose)) return false;
        for(auto wp = working_points.rbegin(); wp != working_points.rend(); ++wp) {
            if(cand.GetIsolation() < wp->second) {
                region.SetLowerIso(wp->first);
                if(wp != working_points.rbegin())
                    region.SetUpperIso((--wp)->first);
                break;
            }
        }
    } else {
        throw exception("BaseEventAnalyzer::SetRegionIsoRange: leg type not supported.");
    }
    return region.HasLowerIso();
}

} // namespace analysis
