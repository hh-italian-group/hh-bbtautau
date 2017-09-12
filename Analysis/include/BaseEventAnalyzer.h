/*! Definition of BaseEventAnalyzer class, the base class for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "EventAnalyzerDataCollection.h"
#include "SampleDescriptorConfigEntryReader.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"

namespace analysis {

struct AnalyzerArguments {
    REQ_ARG(std::string, setup);
    REQ_ARG(std::string, sources);
    REQ_ARG(std::string, input);
    REQ_ARG(std::string, output);
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
    using PhysicalValueMap = std::map<EventRegion, PhysicalValue>;

    static constexpr Channel ChannelId() { return ChannelInfo::IdentifyChannel<FirstLeg, SecondLeg>(); }

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
            for(const auto& n_bjets : bjet_counts)
                categories.emplace(2, n_bjets.second, n_bjets.first);
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

        CombineSampleDescriptorConfigEntryReader combine_entry_reader(combine_descriptors,sample_descriptors);
        config_reader.AddEntryReader("SAMPLE_CMB", combine_entry_reader, false);

        config_reader.ReadConfig(args.sources());

        if(!ana_setup_collection.count(args.setup()))
            throw exception("Setup '%1%' not found in the configuration file '%2%'.") % args.setup() % args.sources();
        ana_setup = ana_setup_collection.at(args.setup());
    }

    virtual ~BaseEventAnalyzer() {}

    void Run()
    {
        ProcessSamples(ana_setup.signals, "signals");
        ProcessSamples(ana_setup.backgrounds, "backgrounds");
        ProcessSamples(ana_setup.data, "data");
        std::cout << "Saving output file..." << std::endl;
    }

protected:
    virtual EventRegion DetermineEventRegion(EventInfo& event, EventCategory eventCategory) = 0;

    const std::string& ChannelName() const { return __Channel_names<>::names.EnumToString(ChannelId()); }
    const std::string& ChannelNameLatex() const { return __Channel_names_latex.EnumToString(ChannelId()); }

    void ProcessSamples(const std::vector<std::string>& sample_names, const std::string& sample_set_name)
    {
        std::cout << "Processing " << sample_set_name << " samples... " << std::endl;
        for(const auto& sample_name : sample_names) {
            if(!sample_descriptors.count(sample_name))
                throw exception("Sample '%1%' not found.") % sample_name;
            const SampleDescriptor& sample = sample_descriptors.at(sample_name);
            if(sample.channels.size() && !sample.channels.count(ChannelId())) continue;
            std::cout << sample.name << std::endl;
            std::vector<std::pair<std::string, std::string>> fileNames;
            for(const auto& file : sample.file_paths)
                fileNames.emplace_back(sample.name, file);
            for(size_t n = 0; n < sample.GetNSignalPoints(); ++n)
                fileNames.emplace_back(sample.GetFullName(n), sample.GetFileName(n));
            for(const auto& fileEntry : fileNames) {
                auto file = root_ext::OpenRootFile(tools::FullPath({args.input(), fileEntry.second}));
                auto tuple = ntuple::CreateEventTuple(ChannelName(), file.get(), true,
                                                      ntuple::TreeState::Skimmed);
                auto summary_tuple = ntuple::CreateSummaryTuple("summary", file.get(), true,
                                                                ntuple::TreeState::Skimmed);
                const auto prod_summary = ntuple::MergeSummaryTuple(*summary_tuple);
                ProcessDataSource(sample, tuple, prod_summary, fileEntry.first);
            }
        }
    }

    void ProcessDataSource(const SampleDescriptor& sample, std::shared_ptr<ntuple::EventTuple> tuple,
                           const ntuple::ProdSummary& prod_summary, const std::string& full_name)
    {
        auto summary = std::make_shared<SummaryInfo>(prod_summary);
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
                                                     event.GetEnergyScale(), full_name);
                    anaDataCollection.Fill(anaDataId, event, event->weight_total * sample.cross_section);
                }
            }
        }
    }

protected:
    AnalyzerArguments args;
    AnaDataCollection anaDataCollection;
    AnalyzerSetup ana_setup;
    SampleDescriptorCollection sample_descriptors;
    CombineSampleDescriptorCollection combine_descriptors;
};

} // namespace analysis
