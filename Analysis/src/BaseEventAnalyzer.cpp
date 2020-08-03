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
    EventCandidate::InitializeUncertainties(ana_setup.period, false, args.working_path(),
                                            signalObjectSelector.GetTauVSjetDiscriminator().first);
    InitializeMvaReader();
    if(ana_setup.syncDataIds.size()){
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

std::pair<EventCategorySet,
          bbtautau::AnaTupleWriter::CategoriesFlags> BaseEventAnalyzer::DetermineEventCategories(EventInfo& event,
                                                                                                 bool pass_vbf_trigger)
{

    static const std::map<DiscriminatorWP, size_t> btag_working_points = {{DiscriminatorWP::Loose, 0},
                                                                          {DiscriminatorWP::Medium, 0},
                                                                          {DiscriminatorWP::Tight, 0}};
    EventCategorySet categories;
    bbtautau::AnaTupleWriter::CategoriesFlags categories_flags;

    // const bool is_boosted =  false;
    auto fatJet = SignalObjectSelector::SelectFatJet(event.GetEventCandidate(), event.GetSelectedSignalJets());
    const bool is_boosted = fatJet != nullptr;
    bool is_VBF = false;
    boost::optional<DiscriminatorWP> vbf_tag;
    std::set<size_t> jets_to_exclude;
    if(event.HasVBFjetPair()){
        const auto vbf_jet_1 = event.GetVBFJet(1).GetMomentum();
        const auto vbf_jet_2 = event.GetVBFJet(2).GetMomentum();
        const auto deta_jj =  std::abs(vbf_jet_1.Eta() - vbf_jet_2.Eta());
        const auto m_jj = (vbf_jet_1 + vbf_jet_2).M();

        is_VBF = (m_jj > cuts::hh_bbtautau_Run2::VBF::mass_jj
                    && deta_jj > cuts::hh_bbtautau_Run2::VBF::deltaeta_jj);
        if(is_VBF) {
            const bool is_tight = m_jj > cuts::hh_bbtautau_Run2::VBF::mass_jj_tight && pass_vbf_trigger;
            vbf_tag = is_tight ? DiscriminatorWP::Tight : DiscriminatorWP::Loose;
        }
        jets_to_exclude.insert(event.GetVBFJet(1)->jet_index());
        jets_to_exclude.insert(event.GetVBFJet(2)->jet_index());
    }

    const auto& all_jets = event.GetCentralJets();
    auto bjet_counts = btag_working_points;
    if(event.HasBjetPair()) {
        for(size_t bjet_index = 1; bjet_index <= 2; ++bjet_index) {
            const auto& jet = event.GetBJet(bjet_index);
            for(const auto& btag_wp : btag_working_points) {
                if(bTagger->Pass(*jet, btag_wp.first))
                    ++bjet_counts[btag_wp.first];
            }
        }
    }
    categories_flags.num_jets = all_jets.size();
    categories_flags.num_btag_loose = bjet_counts.at(DiscriminatorWP::Loose);
    categories_flags.num_btag_medium = bjet_counts.at(DiscriminatorWP::Medium);
    categories_flags.num_btag_tight = bjet_counts.at(DiscriminatorWP::Tight);
    categories_flags.is_vbf = is_VBF;
    categories_flags.is_boosted = is_boosted;
    categories_flags.fat_jet_cand = fatJet;

    for(const auto& category : ana_setup.categories_base) {
        if(category.Contains(all_jets.size(), bjet_counts, is_VBF, is_boosted, vbf_tag))
            categories.insert(category);
    }
    return std::make_pair(categories, categories_flags);
}

void BaseEventAnalyzer::InitializeMvaReader()
{
    using MvaKey = mva_study::MvaReader::MvaKey;
    if(!mva_setup.is_initialized()) return;
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
            mva_reader.Add(key, FullPath(file), vars, legacy, legacy_lm);
        }
    }
}

EventSubCategory BaseEventAnalyzer::DetermineEventSubCategory(EventInfo& event, const EventCategory& category,
                                                              std::map<SelectionCut, double>& mva_scores)
{
    using namespace cuts::hh_bbtautau_Run2::hh_tag;
    using MvaKey = mva_study::MvaReader::MvaKey;

    EventSubCategory sub_category;
    if(event.HasBjetPair()) {
        const double mbb = event.GetHiggsBB().GetMomentum().mass();

        if(category.HasBoostConstraint() && category.IsBoosted()) {
            if(ana_setup.use_svFit) {
                const auto htt = event.GetHiggsTTMomentum(true, ana_setup.allow_calc_svFit);
                const bool isInsideBoostedCut = htt && IsInsideBoostedMassWindow(htt->mass(), mbb);
                sub_category.SetCutResult(SelectionCut::mh, isInsideBoostedCut);
            }
        } else {
            if(!ana_setup.use_svFit && ana_setup.massWindowParams.count(SelectionCut::mh))
                throw exception("Category mh inconsistent with the false requirement of SVfit.");
            if(ana_setup.massWindowParams.count(SelectionCut::mh)) {
                const bool cut_result = ana_setup.use_svFit
                    && event.GetSVFitResults(ana_setup.allow_calc_svFit).has_valid_momentum
                    && ana_setup.massWindowParams.at(SelectionCut::mh).IsInside(
                            event.GetHiggsTTMomentum(true, ana_setup.allow_calc_svFit)->mass(), mbb);
                sub_category.SetCutResult(SelectionCut::mh, cut_result);
            }
            if(ana_setup.massWindowParams.count(SelectionCut::mhVis))
                sub_category.SetCutResult(SelectionCut::mhVis,ana_setup.massWindowParams.at(SelectionCut::mhVis)
                        .IsInside(event.GetHiggsTTMomentum(false)->mass(), mbb));

            if(ana_setup.massWindowParams.count(SelectionCut::mhMET))
                sub_category.SetCutResult(SelectionCut::mhMET,ana_setup.massWindowParams.at(SelectionCut::mhMET)
                        .IsInside((*event.GetHiggsTTMomentum(false) + event.GetMET().GetMomentum()).mass(), mbb));

        }
        if(ana_setup.use_kinFit)
            sub_category.SetCutResult(SelectionCut::KinematicFitConverged,
                                      event.GetKinFitResults(ana_setup.allow_calc_svFit).HasValidMass());
    }
    if(mva_setup.is_initialized()) {

        std::map<MvaKey, std::future<double>> scores;
        for(const auto& mva_sel : mva_setup->selections) {
            const auto& params = mva_sel.second;
            const MvaKey key{params.name, static_cast<int>(params.mass), params.spin};
            if(!scores.count(key)) {
                auto eval = std::bind(&mva_study::MvaReader::Evaluate, &mva_reader, key, &event);
                scores[key] = run::async(eval);
            }
        }
        for(const auto& mva_sel : mva_setup->selections) {
            const auto& params = mva_sel.second;
            const MvaKey key{params.name, static_cast<int>(params.mass),params.spin};
            const double score = scores.at(key).get();
            const bool pass = score > params.cut;
            sub_category.SetCutResult(mva_sel.first, pass);
            mva_scores[mva_sel.first] = score;
        }
    }

    return sub_category;
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
                nonResModel = std::make_shared<NonResModel>(ana_setup.period, sample, file, *crossSectionProvider);
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

            if(!event) continue;
            const bool pass_normal_trigger = event->PassNormalTriggers();
            const bool pass_vbf_trigger = event->PassVbfTriggers();
            const bool pass_trigger = pass_normal_trigger || pass_vbf_trigger;
            if(!pass_trigger) continue;
            std::map<DiscriminatorWP, std::map<UncertaintyScale, float>> btag_weights;
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

                auto jet_pu_id_weight_provided = eventWeights_HH->GetProviderT<mc_corrections::JetPuIdWeights>(
                        mc_corrections::WeightType::JetPuIdWeights);
                jet_pu_id_weight = jet_pu_id_weight_provided->Get(*event);

                auto top_pt_weight_provided = eventWeights_HH->GetProviderT<mc_corrections::TopPtWeight>(
                        mc_corrections::WeightType::TopPt);

                auto btag_weight_provider = eventWeights_HH->GetProviderT<mc_corrections::BTagWeight>(
                    mc_corrections::WeightType::BTag);

                for(const auto wp : btag_wps){
                    btag_weights[wp][UncertaintyScale::Central] = static_cast<float>(btag_weight_provider->Get(*event,
                                                                                     wp, unc_source, unc_scale));
                    if(unc_source == UncertaintySource::None) {
                        btag_weights[wp][UncertaintyScale::Up] = static_cast<float>(btag_weight_provider->Get(*event, wp,
                                UncertaintySource::Eff_b, UncertaintyScale::Up));
                        btag_weights[wp][UncertaintyScale::Down] = static_cast<float>(btag_weight_provider->Get(*event,
                                wp, UncertaintySource::Eff_b, UncertaintyScale::Down));
                    }
                }
                if(unc_source == UncertaintySource::None) {
                    static const std::vector<UncertaintySource> uncs_trigger_weight = { UncertaintySource::EleTriggerUnc,
                        UncertaintySource::MuonTriggerUnc, UncertaintySource::TauTriggerUnc_DM0,
                        UncertaintySource::TauTriggerUnc_DM1, UncertaintySource::TauTriggerUnc_DM10,
                        UncertaintySource::TauTriggerUnc_DM11 };

                    static const std::vector<UncertaintySource> uncs_tau_weight = { UncertaintySource::TauVSjetSF_DM0,
                        UncertaintySource::TauVSjetSF_DM1, UncertaintySource::TauVSjetSF_3prong,
                        UncertaintySource::TauVSjetSF_pt20to25, UncertaintySource::TauVSjetSF_pt25to30,
                        UncertaintySource::TauVSjetSF_pt30to35, UncertaintySource::TauVSjetSF_pt35to40,
                        UncertaintySource::TauVSjetSF_ptgt40, UncertaintySource::TauVSeSF_barrel,
                        UncertaintySource::TauVSeSF_endcap, UncertaintySource::TauVSmuSF_etaLt0p4,
                        UncertaintySource::TauVSmuSF_eta0p4to0p8, UncertaintySource::TauVSmuSF_eta0p8to1p2,
                        UncertaintySource::TauVSmuSF_eta1p2to1p7, UncertaintySource::TauVSmuSF_etaGt1p7,
                        UncertaintySource::EleIdIsoUnc, UncertaintySource::MuonIdIsoUnc };

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
                    if(sample.apply_top_pt_unc){
                        uncs_weight_map[UncertaintySource::TopPt][UncertaintyScale::Up] =
                            static_cast<float>(top_pt_weight_provided->Get(*event));
                        uncs_weight_map[UncertaintySource::TopPt][UncertaintyScale::Central] = 1;
                    }
                }
            }

            const auto [eventCategories, categories_flags] = DetermineEventCategories(*event, pass_vbf_trigger);
            for(auto eventCategory : eventCategories) {
                const EventRegion eventRegion = DetermineEventRegion(*event, eventCategory);
                for(const auto& region : ana_setup.regions){
                    if(!eventRegion.Implies(region)) continue;
                    std::map<SelectionCut, double> mva_scores;
                    const auto eventSubCategory = DetermineEventSubCategory(*event, eventCategory, mva_scores);
                    for(const auto& subCategory : sub_categories_to_process) {
                        if(!eventSubCategory.Implies(subCategory)) continue;
                        SelectionCut mva_cut;
                        double mva_score = 0;
                        if(subCategory.TryGetLastMvaCut(mva_cut))
                            mva_score = mva_scores.at(mva_cut);
                        event->SetMvaScore(mva_score);
                        const EventAnalyzerDataId anaDataId(eventCategory, subCategory, region,
                                                            unc_source, unc_scale, sample_wp.full_name);
                        double shape_weight = 1.;
                        if(sample.sampleType == SampleType::Data) {
                            dataIds[anaDataId] = std::make_tuple(1., mva_score);
                        } else {

                            if(eventCategory.HasBtagConstraint())
                                btag_weight = btag_weights.at(eventCategory.BtagWP()).at(UncertaintyScale::Central);

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
                                dataIds[anaDataId] = std::make_tuple(weight, mva_score);
                            } else
                                ProcessSpecialEvent(sample, sample_wp, anaDataId, *event, weight,
                                                    (*summary)->totalShapeWeight, dataIds, cross_section);
                        }

                        for(size_t n = 0; n < sync_descriptors.size(); ++n) {
                            if(sync_event_selected[n]) continue;
                            const auto& regex_pattern = sync_descriptors.at(n).regex_pattern;
                            for(auto& dataId : dataIds) {
                                if(!boost::regex_match(dataId.first.GetName(), *regex_pattern)) continue;

                                htt_sync::FillSyncTuple(*event, *sync_descriptors.at(n).sync_tree, ana_setup.period,
                                                        ana_setup.use_svFit, std::get<0>(dataId.second),
                                                        lepton_id_iso_weight, trigger_weight, btag_weight,
                                                        shape_weight, jet_pu_id_weight, prescale_weight);
                                sync_event_selected[n] = true;
                                break;
                            }
                        }
                    }
                }
            }
            //dataId
            anaTupleWriter.AddEvent(*event, dataIds, pass_vbf_trigger, categories_flags, btag_weights, uncs_weight_map);
        }
    }
}

void BaseEventAnalyzer::ProcessSpecialEvent(const SampleDescriptor& sample,
                                            const SampleDescriptor::Point& /*sample_wp*/,
                                            const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight,
                                            double shape_weight, bbtautau::AnaTupleWriter::DataIdMap& dataIds,
                                            double cross_section)
{
    if(sample.sampleType == SampleType::DY){
        dymod.at(sample.name)->ProcessEvent(anaDataId,event,weight,dataIds);
    } else if(sample.sampleType == SampleType::ggHH_NonRes) {
        nonResModel->ProcessEvent(anaDataId, event, weight, shape_weight, dataIds, cross_section);
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
