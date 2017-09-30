/*! Definition of BaseEventAnalyzer class, the base class for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Run/include/program_main.h"
#include "EventAnalyzerDataCollection.h"
#include "SampleDescriptorConfigEntryReader.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"
#include "StackedPlotsProducer.h"

namespace analysis {

struct AnalyzerArguments {
    REQ_ARG(std::string, setup);
    REQ_ARG(std::string, sources);
    REQ_ARG(std::string, input);
    REQ_ARG(std::string, output);
    OPT_ARG(bool, draw, true);
    OPT_ARG(bool, saveFullOutput, false);
    OPT_ARG(unsigned, n_threads, 1);
};

template<typename _FirstLeg, typename _SecondLeg>
class BaseEventAnalyzer {
public:
    using FirstLeg = _FirstLeg;
    using SecondLeg = _SecondLeg;
    using EventInfo = ::analysis::EventInfo<FirstLeg, SecondLeg>;
    using AnaData = ::analysis::EventAnalyzerData<FirstLeg, SecondLeg>;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection<AnaData>;
    using PlotsProducer = ::analysis::StackedPlotsProducer<FirstLeg, SecondLeg>;

    static constexpr Channel ChannelId() { return ChannelInfo::IdentifyChannel<FirstLeg, SecondLeg>(); }
    static const std::string& ChannelNameLatex() { return __Channel_names_latex.EnumToString(ChannelId()); }

    virtual const EventCategorySet& EventCategoriesToProcess() const
    {
        static const EventCategorySet categories = {
            EventCategory::TwoJets_Inclusive(), EventCategory::TwoJets_ZeroBtag(),
            EventCategory::TwoJets_OneBtag(), EventCategory::TwoJets_OneLooseBtag(),
            EventCategory::TwoJets_TwoBtag(), EventCategory::TwoJets_TwoLooseBtag()
        };
        return categories;
    }

    virtual const EventSubCategorySet& EventSubCategoriesToProcess() const
    {
        static const EventSubCategorySet sub_categories = {
            EventSubCategory().SetCutResult(SelectionCut::InsideMassWindow, true)
//                              .SetCutResult(SelectionCut::MVA, true)
                              .SetCutResult(SelectionCut::KinematicFitConverged, true)
        };
        return sub_categories;
    }

    virtual const EventRegionSet& EventRegionsToProcess() const
    {
        static const EventRegionSet regions = {
            EventRegion::OS_Isolated(), EventRegion::OS_AntiIsolated(),
            EventRegion::SS_Isolated(), EventRegion::SS_AntiIsolated()
        };
        return regions;
    }

    virtual const EventEnergyScaleSet& EventEnergyScaleToProcess() const
    {
        static const EventEnergyScaleSet scales = {
            EventEnergyScale::Central, EventEnergyScale::TauUp, EventEnergyScale::TauDown,
            EventEnergyScale::JetUp, EventEnergyScale::JetDown
        };
        return scales;
    }

    static EventCategorySet DetermineEventCategories(EventInfo& event)
    {
        static const std::map<DiscriminatorWP, double> btag_working_points = {
            { DiscriminatorWP::Loose, cuts::btag_2016::CSVv2L },
            { DiscriminatorWP::Medium, cuts::btag_2016::CSVv2M }
        };

        EventCategorySet categories;
        categories.insert(EventCategory::Inclusive());

        if(event.HasBjetPair()) {
            categories.insert(EventCategory::TwoJets_Inclusive());
            const std::vector<const JetCandidate*> jets = {
                &event.GetHiggsBB().GetFirstDaughter(), &event.GetHiggsBB().GetSecondDaughter(),
            };

            std::map<DiscriminatorWP, size_t> bjet_counts;
            for(const auto& jet : jets) {
                for(const auto& btag_wp : btag_working_points) {
                    if((*jet)->csv() > btag_wp.second)
                        ++bjet_counts[btag_wp.first];
                }
            }
            for(const auto& wp_entry : btag_working_points)
                categories.emplace(2, bjet_counts[wp_entry.first], wp_entry.first);
        }
        return categories;
    }

    static EventSubCategory DetermineEventSubCategory(EventInfo& event)
    {
        using namespace cuts::hh_bbtautau_2016::hh_tag;
        EventSubCategory sub_category;
        sub_category.SetCutResult(SelectionCut::InsideMassWindow,
                                  IsInsideEllipse(event.GetHiggsTT(true).GetMomentum().mass(),
                                                  event.GetHiggsBB().GetMomentum().mass()));
        sub_category.SetCutResult(SelectionCut::KinematicFitConverged,
                                  event.GetKinFitResults().HasValidMass());
        return sub_category;
    }

    BaseEventAnalyzer(const AnalyzerArguments& _args)
        : args(_args), anaDataCollection(args.output() + "_full.root", args.saveFullOutput())
    {
        ROOT::EnableThreadSafety();
        if(args.n_threads() > 1)
            ROOT::EnableImplicitMT(args.n_threads());

        ConfigReader config_reader;

        AnalyzerSetupCollection ana_setup_collection;
        AnalyzerConfigEntryReader ana_entry_reader(ana_setup_collection);
        config_reader.AddEntryReader("ANA_DESC", ana_entry_reader, false);

        SampleDescriptorConfigEntryReader sample_entry_reader(sample_descriptors);
        config_reader.AddEntryReader("SAMPLE", sample_entry_reader, true);

        CombineSampleDescriptorConfigEntryReader combine_entry_reader(cmb_sample_descriptors,sample_descriptors);
        config_reader.AddEntryReader("SAMPLE_CMB", combine_entry_reader, false);

        config_reader.ReadConfig(args.sources());

        if(!ana_setup_collection.count(args.setup()))
            throw exception("Setup '%1%' not found in the configuration file '%2%'.") % args.setup() % args.sources();
        ana_setup = ana_setup_collection.at(args.setup());

        RemoveUnusedSamples();
    }

    virtual ~BaseEventAnalyzer() {}

    void Run()
    {
        ProcessSamples(ana_setup.signals, "signals");
        ProcessSamples(ana_setup.data, "data");
        ProcessSamples(ana_setup.backgrounds, "backgrounds");

        if(args.draw()) {
            std::cout << "Creating plots..." << std::endl;
            const auto samplesToDraw = PlotsProducer::CreateOrderedSampleCollection(
                        ana_setup.draw_sequence, sample_descriptors, cmb_sample_descriptors, ana_setup.signals,
                        ana_setup.data);
            PlotsProducer plotsProducer(anaDataCollection, samplesToDraw, false, true);
            std::set<std::string> signal_names(ana_setup.signals.begin(), ana_setup.signals.end());
            plotsProducer.PrintStackedPlots(args.output(), EventRegion::SignalRegion(), EventCategoriesToProcess(),
                                            EventSubCategoriesToProcess(), signal_names);
        }

        std::cout << "Saving output file..." << std::endl;
    }

protected:
    virtual EventRegion DetermineEventRegion(EventInfo& event, EventCategory eventCategory) = 0;

    void ProcessSamples(const std::vector<std::string>& sample_names, const std::string& sample_set_name)
    {
        std::cout << "Processing " << sample_set_name << " samples... " << std::endl;
        for(size_t sample_index = 0; sample_index < sample_names.size(); ++sample_index) {
            const std::string& sample_name = sample_names.at(sample_index);
            if(!sample_descriptors.count(sample_name))
                throw exception("Sample '%1%' not found.") % sample_name;
            SampleDescriptor& sample = sample_descriptors.at(sample_name);
            if(sample.channels.size() && !sample.channels.count(ChannelId())) continue;
            std::cout << sample.name << std::endl;
            sample.CreateWorkingPoints();
            if(sample.sampleType == SampleType::QCD) {
                if(sample_index != sample_names.size() - 1)
                    throw exception("QCD sample should be the last in the background list.");
                EstimateQCD(sample);
                continue;
            }
            std::set<std::string> processed_files;
            for(const auto& sample_wp : sample.working_points) {
                if(!sample_wp.file_path.size() || processed_files.count(sample_wp.file_path)) continue;
                auto file = root_ext::OpenRootFile(tools::FullPath({args.input(), sample_wp.file_path}));
                auto tuple = ntuple::CreateEventTuple(ToString(ChannelId()), file.get(), true,
                                                      ntuple::TreeState::Skimmed);
                auto summary_tuple = ntuple::CreateSummaryTuple("summary", file.get(), true,
                                                                ntuple::TreeState::Skimmed);
                const auto prod_summary = ntuple::MergeSummaryTuple(*summary_tuple);
                ProcessDataSource(sample, sample_wp, tuple, prod_summary);
            }
        }
    }

    void ProcessDataSource(const SampleDescriptor& sample, const SampleDescriptor::Point& sample_wp,
                           std::shared_ptr<ntuple::EventTuple> tuple, const ntuple::ProdSummary& prod_summary)
    {
        const SummaryInfo summary(prod_summary);
        for(const auto& tupleEvent : *tuple) {
            EventInfo event(tupleEvent, ntuple::JetPair{0, 1}, summary);
            if(!EventEnergyScaleToProcess().count(event.GetEnergyScale())) continue;
            const auto eventCategories = DetermineEventCategories(event);
            for(auto eventCategory : eventCategories) {
                if (!EventCategoriesToProcess().count(eventCategory)) continue;
                const EventRegion eventRegion = DetermineEventRegion(event, eventCategory);
                if(!EventRegionsToProcess().count(eventRegion)) continue;

                const auto eventSubCategory = DetermineEventSubCategory(event);
                for(const auto& subCategory : EventSubCategoriesToProcess()) {
                    if(!eventSubCategory.Implies(subCategory)) continue;
                    const EventAnalyzerDataId anaDataId(eventCategory, subCategory, eventRegion,
                                                        event.GetEnergyScale(), sample_wp.full_name);
                    if(sample.sampleType == SampleType::Data) {
                        anaDataCollection.Fill(anaDataId, event, 1);
                    } else {
                        const double weight = event->weight_total * sample.cross_section * ana_setup.int_lumi
                                            / summary->totalShapeWeight;
                        if(sample.sampleType == SampleType::MC)
                            anaDataCollection.Fill(anaDataId, event, weight);
                        else
                            ProcessSpecialEvent(sample, sample_wp, anaDataId, event, weight);
                    }
                }
            }
        }
    }

    virtual void EstimateQCD(const SampleDescriptor& /*qcd_sample*/)
    {
        // TODO: implement QCD estimation.
    }

    virtual void ProcessSpecialEvent(const SampleDescriptor& /*sample*/, const SampleDescriptor::Point& /*sample_wp*/,
                                     const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight)
    {
        // TODO: implement DY estimation.
        anaDataCollection.Fill(anaDataId, event, weight);
    }

    template<typename SampleCollection>
    static std::vector<std::string> FilterInactiveSamples(const SampleCollection& samples,
                                                          const std::vector<std::string>& sample_names)
    {
        std::vector<std::string> filtered;
        for(const auto& name : sample_names) {
            if(!samples.count(name))
                throw exception("Sample '%1%' not found.") % name;
            const auto& sample = samples.at(name);
            if(sample.channels.size() && !sample.channels.count(ChannelId())) continue;
            filtered.push_back(name);
        }
        return filtered;
    }

    void RemoveUnusedSamples()
    {
        RemoveUnusedSamples(sample_descriptors, { &ana_setup.signals, &ana_setup.backgrounds, &ana_setup.data });
        RemoveUnusedSamples(cmb_sample_descriptors, { &ana_setup.cmb_samples });
    }

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
    AnalyzerArguments args;
    AnaDataCollection anaDataCollection;
    AnalyzerSetup ana_setup;
    SampleDescriptorCollection sample_descriptors;
    CombineSampleDescriptorCollection cmb_sample_descriptors;
};

} // namespace analysis
