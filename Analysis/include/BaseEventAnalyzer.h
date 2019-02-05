/*! Definition of BaseEventAnalyzer class, the base class for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTools/Run/include/MultiThread.h"
#include "SampleDescriptorConfigEntryReader.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2017.h"
#include "h-tautau/Analysis/include/EventLoader.h"
#include "h-tautau/Analysis/include/SyncTupleHTT.h"
#include "h-tautau/Analysis/include/BTagger.h"
#include "MvaReader.h"
#include "EventAnalyzerData.h"
#include "AnaTuple.h"
#include "EventAnalyzerCore.h"
#include "Analysis/include/DYModel.h"
#include "NonResModel.h"

namespace analysis {

struct AnalyzerArguments : CoreAnalyzerArguments {
    REQ_ARG(std::string, input);
    REQ_ARG(std::string, output);
    OPT_ARG(std::string, output_sync,"sync.root");
};

template<typename _FirstLeg, typename _SecondLeg>
class BaseEventAnalyzer : public EventAnalyzerCore {
public:
    using FirstLeg = _FirstLeg;
    using SecondLeg = _SecondLeg;
    using Event = ntuple::Event;
    using EventPtr = std::shared_ptr<Event>;
    using EventInfo = ::analysis::EventInfo<FirstLeg, SecondLeg>;

    static constexpr Channel ChannelId() { return ChannelInfo::IdentifyChannel<FirstLeg, SecondLeg>(); }

    EventCategorySet DetermineEventCategories(EventInfo& event)
    {
        size_t b_jets = 0;
        static const std::map<DiscriminatorWP, size_t> btag_working_points = {{DiscriminatorWP::Loose, b_jets},
                                                                              {DiscriminatorWP::Medium, b_jets},
                                                                              {DiscriminatorWP::Tight, b_jets }};
        EventCategorySet categories;

        const bool is_boosted = event.SelectFatJet(cuts::hh_bbtautau_2016::fatJetID::mass,
                                                   cuts::hh_bbtautau_2016::fatJetID::deltaR_subjet) != nullptr;
        bool is_VBF = false;
        if(event.HasVBFjetPair()){
            const auto vbf_jet_1 = event.GetVBFJet(1).GetMomentum();
            const auto vbf_jet_2 = event.GetVBFJet(2).GetMomentum();

            const auto deta_jj =  std::abs(vbf_jet_1.Eta() - vbf_jet_2.Eta());
            const auto m_jj = (vbf_jet_1 + vbf_jet_2).M();

            is_VBF = (m_jj > cuts::hh_bbtautau_2017::VBF::mass_jj
                        && deta_jj > cuts::hh_bbtautau_2017::VBF::deltaeta_jj);
        }

        const auto& jets = event.SelectJets(cuts::btag_2017::pt,cuts::btag_2017::eta, ana_setup.jet_ordering);
        std::map<DiscriminatorWP, size_t> bjet_counts = btag_working_points;

        for(const auto& jet : jets) {
            for(const auto& btag_wp : btag_working_points){
                if(bTagger->Pass(*jet, btag_wp.first)) ++bjet_counts[btag_wp.first];
            }
        }

        for(const auto& category : ana_setup.categories) {
            if(category.Contains(jets.size(), bjet_counts, is_VBF ,is_boosted))
                categories.insert(category);
        }
        return categories;
    }

    BaseEventAnalyzer(const AnalyzerArguments& _args) :
        EventAnalyzerCore(_args, ChannelId()), args(_args), anaTupleWriter(args.output(), ChannelId()),
        trigger_patterns(ana_setup.trigger.at(ChannelId()))
    {
        InitializeMvaReader();
        if(ana_setup.syncDataIds.size()){
            outputFile_sync = root_ext::CreateRootFile(args.output_sync());
            for(unsigned n = 0; n < ana_setup.syncDataIds.size(); ++n){
                const EventAnalyzerDataId dataId = ana_setup.syncDataIds.at(n);
                syncTuple_map[dataId] = std::make_shared<htt_sync::SyncTuple>(dataId.GetName("_"),outputFile_sync.get(),false);
            }

        }
        for(auto& x: sample_descriptors) {
            SampleDescriptor& sample = x.second;
             if(sample.sampleType == SampleType::DY)
                 dymod = std::make_shared<DYModel>(sample,args.working_path());
        }
    }

    void Run()
    {
        run::ThreadPull threads(args.n_threads());
        ProcessSamples(ana_setup.signals, "signal");
        ProcessSamples(ana_setup.data, "data");
        ProcessSamples(ana_setup.backgrounds, "background");
        std::cout << "Saving output file..." << std::endl;
        for (auto& sync_iter : syncTuple_map){
            sync_iter.second->Write();
        }
    }

protected:
    virtual EventRegion DetermineEventRegion(EventInfo& event, EventCategory eventCategory) = 0;

    void InitializeMvaReader()
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

    virtual EventSubCategory DetermineEventSubCategory(EventInfo& event, const EventCategory& category,
                                                       std::map<SelectionCut, double>& mva_scores)
    {
        using namespace cuts::hh_bbtautau_2016::hh_tag;
        using MvaKey = mva_study::MvaReader::MvaKey;

        EventSubCategory sub_category;

        double mbb = 0;
        if(event.HasBjetPair()){
            mbb = event.GetHiggsBB().GetMomentum().mass();
            if(category.HasBoostConstraint() && category.IsBoosted()){
                bool isInsideBoostedCut = IsInsideBoostedMassWindow(event.GetHiggsTT(true).GetMomentum().mass(),mbb);
                sub_category.SetCutResult(SelectionCut::mh,isInsideBoostedCut);
                sub_category.SetCutResult(SelectionCut::mhVis,isInsideBoostedCut);
                sub_category.SetCutResult(SelectionCut::mhMET,isInsideBoostedCut);
            }
            else{
                if(ana_setup.massWindowParams.count(SelectionCut::mh))
                    sub_category.SetCutResult(SelectionCut::mh,ana_setup.massWindowParams.at(SelectionCut::mh)
                            .IsInside(event.GetHiggsTT(true).GetMomentum().mass(),mbb));
                if(ana_setup.massWindowParams.count(SelectionCut::mhVis))
                    sub_category.SetCutResult(SelectionCut::mhVis,ana_setup.massWindowParams.at(SelectionCut::mhVis)
                            .IsInside(event.GetHiggsTT(false).GetMomentum().mass(),mbb));
                if(ana_setup.massWindowParams.count(SelectionCut::mhMET))
                    sub_category.SetCutResult(SelectionCut::mhMET,ana_setup.massWindowParams.at(SelectionCut::mhMET)
                            .IsInside((event.GetHiggsTT(false).GetMomentum() + event.GetMET().GetMomentum()).mass(),mbb));
            }
            sub_category.SetCutResult(SelectionCut::KinematicFitConverged,
                                      event.GetKinFitResults().HasValidMass());
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

    void ProcessSamples(const std::vector<std::string>& sample_names, const std::string& sample_set_name)
    {
        std::cout << "Processing " << sample_set_name << " samples... " << std::endl;
        for(size_t sample_index = 0; sample_index < sample_names.size(); ++sample_index) {
            const std::string& sample_name = sample_names.at(sample_index);
            if(!sample_descriptors.count(sample_name))
                throw exception("Sample '%1%' not found.") % sample_name;
            SampleDescriptor& sample = sample_descriptors.at(sample_name);
            if(sample.sampleType == SampleType::QCD || (sample.channels.size() && !sample.channels.count(ChannelId())))
                continue;
            std::cout << '\t' << sample.name << std::endl;

            std::set<std::string> processed_files;
            for(const auto& sample_wp : sample.working_points) {
                if(!sample_wp.file_path.size() || processed_files.count(sample_wp.file_path)) continue;
                auto file = root_ext::OpenRootFile(tools::FullPath({args.input(), sample_wp.file_path}));
                auto tuple = ntuple::CreateEventTuple(ToString(ChannelId()), file.get(), true,
                                                      ntuple::TreeState::Skimmed);
                auto summary_tuple = ntuple::CreateSummaryTuple("summary", file.get(), true,
                                                                ntuple::TreeState::Skimmed);
                const auto prod_summary = ntuple::MergeSummaryTuple(*summary_tuple);
                if(sample.sampleType == SampleType::NonResHH) {
                    std::cout << "\t\tpreparing NonResModel... ";
                    std::cout.flush();
                    nonResModel = std::make_shared<NonResModel>(ana_setup.period, sample, file);
                    std::cout << "done." << std::endl;
                }
                ProcessDataSource(sample, sample_wp, tuple, prod_summary);

                processed_files.insert(sample_wp.file_path);
            }
        }
    }

    void ProcessDataSource(const SampleDescriptor& sample, const SampleDescriptor::Point& sample_wp,
                           std::shared_ptr<ntuple::EventTuple> tuple, const ntuple::ProdSummary& prod_summary)
    {
        const SummaryInfo summary(prod_summary);
        Event prevFullEvent, *prevFullEventPtr = nullptr;
        for(auto tupleEvent : *tuple) {
            if(ntuple::EventLoader::Load(tupleEvent, prevFullEventPtr).IsFull()) {
                prevFullEvent = tupleEvent;
                prevFullEventPtr = &prevFullEvent;
            }
            EventInfo event(tupleEvent, ana_setup.period, ana_setup.jet_ordering, &summary);
            if(!ana_setup.energy_scales.count(event.GetEnergyScale())) continue;
            if(!event.GetTriggerResults().AnyAcceptAndMatch(trigger_patterns)) continue;
            bbtautau::AnaTupleWriter::DataIdMap dataIds;
            const auto eventCategories = DetermineEventCategories(event);
            for(auto eventCategory : eventCategories) {
                //if (!ana_setup.categories.count(eventCategory)) continue;
                const EventRegion eventRegion = DetermineEventRegion(event, eventCategory);
                for(const auto& region : ana_setup.regions){
                    if(!eventRegion.Implies(region)) continue;
                    std::map<SelectionCut, double> mva_scores;
                    const auto eventSubCategory = DetermineEventSubCategory(event, eventCategory, mva_scores);
                    for(const auto& subCategory : sub_categories_to_process) {
                        if(!eventSubCategory.Implies(subCategory)) continue;
                        SelectionCut mva_cut;
                        double mva_score = 0, mva_weight_scale = 1.;
                        if(subCategory.TryGetLastMvaCut(mva_cut)) {
                            mva_score = mva_scores.at(mva_cut);
                            const auto& mva_params = mva_setup->selections.at(mva_cut);
                            if(mva_params.training_range.is_initialized() && mva_params.samples.count(sample.name)) {
                                if(mva_params.training_range->Contains(event->split_id)) continue;
                                mva_weight_scale = double(summary->n_splits)
                                        / (summary->n_splits - mva_params.training_range->size());
                            }
                        }
                        event.SetMvaScore(mva_score);
                        const EventAnalyzerDataId anaDataId(eventCategory, subCategory, region,
                                                            event.GetEnergyScale(), sample_wp.full_name);
                        if(sample.sampleType == SampleType::Data) {
                            dataIds[anaDataId] = std::make_tuple(1., mva_score);
                        } else {
                            const double weight = event->weight_total * sample.cross_section * ana_setup.int_lumi
                                                / summary->totalShapeWeight * mva_weight_scale;
                            if(sample.sampleType == SampleType::MC) {
                                dataIds[anaDataId] = std::make_tuple(weight, mva_score);
                            } else
                                ProcessSpecialEvent(sample, sample_wp, anaDataId, event, weight,
                                                    summary->totalShapeWeight, dataIds);
                        }
                    }
                }
            }
            anaTupleWriter.AddEvent(event, dataIds);
            for (auto& sync_iter : syncTuple_map) {
                if(!dataIds.count(sync_iter.first)) continue;
                htt_sync::FillSyncTuple(event, *sync_iter.second, ana_setup.period);
            }
        }
    }

    virtual void ProcessSpecialEvent(const SampleDescriptor& sample, const SampleDescriptor::Point& /*sample_wp*/,
                                     const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight,
                                     double shape_weight, bbtautau::AnaTupleWriter::DataIdMap& dataIds)
    {
        if(sample.sampleType == SampleType::DY) {
            dymod->ProcessEvent(anaDataId,event,weight,dataIds);
        } else if(sample.sampleType == SampleType::TT) {
            dataIds[anaDataId] = std::make_tuple(weight, event.GetMvaScore());
            if(anaDataId.Get<EventEnergyScale>() == EventEnergyScale::Central) {
//                const double weight_topPt = event->weight_total * sample.cross_section * ana_setup.int_lumi
//                        / event.GetSummaryInfo()->totalShapeWeight_withTopPt;
                // FIXME
                if(ana_setup.energy_scales.count(EventEnergyScale::TopPtUp))
                    dataIds[anaDataId.Set(EventEnergyScale::TopPtUp)] = std::make_tuple(weight * event->weight_top_pt,
                                                                                        event.GetMvaScore());
                if(ana_setup.energy_scales.count(EventEnergyScale::TopPtDown))
                    dataIds[anaDataId.Set(EventEnergyScale::TopPtDown)] = std::make_tuple(weight * event->weight_top_pt,
                                                                                          event.GetMvaScore());
            }
        } else if(sample.sampleType == SampleType::NonResHH) {
            nonResModel->ProcessEvent(anaDataId, event, weight, shape_weight, dataIds);
        } else
            throw exception("Unsupported special event type '%1%'.") % sample.sampleType;
    }

protected:
    AnalyzerArguments args;
    bbtautau::AnaTupleWriter anaTupleWriter;
    mva_study::MvaReader mva_reader;
    std::shared_ptr<TFile> outputFile_sync;
    std::map<EventAnalyzerDataId, std::shared_ptr<htt_sync::SyncTuple>> syncTuple_map;
    std::shared_ptr<DYModel> dymod;
    std::shared_ptr<NonResModel> nonResModel;
    const std::vector<std::string> trigger_patterns;
};

} // namespace analysis
