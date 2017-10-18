/*! Definition of BaseEventAnalyzer class, the base class for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "EventAnalyzerDataCollection.h"
#include "SampleDescriptorConfigEntryReader.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"
#include "h-tautau/Analysis/include/EventLoader.h"
#include "StackedPlotsProducer.h"
#include "LimitsInputProducer.h"
#include "MvaReader.h"

namespace analysis {

struct AnalyzerArguments {
    REQ_ARG(std::string, setup);
    REQ_ARG(std::string, sources);
    REQ_ARG(std::string, input);
    REQ_ARG(std::string, output);
    OPT_ARG(std::string, mva_setup, "");
    OPT_ARG(bool, draw, true);
    OPT_ARG(bool, saveFullOutput, false);
    OPT_ARG(unsigned, event_set, 0);
    OPT_ARG(unsigned, n_threads, 1);
};

template<typename _FirstLeg, typename _SecondLeg>
class BaseEventAnalyzer {
public:
    using FirstLeg = _FirstLeg;
    using SecondLeg = _SecondLeg;
    using Event = ntuple::Event;
    using EventPtr = std::shared_ptr<Event>;
    using EventInfo = ::analysis::EventInfo<FirstLeg, SecondLeg>;
    using AnaData = ::analysis::EventAnalyzerData<FirstLeg, SecondLeg>;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection<AnaData>;
    using PlotsProducer = ::analysis::StackedPlotsProducer<FirstLeg, SecondLeg>;

    static constexpr Channel ChannelId() { return ChannelInfo::IdentifyChannel<FirstLeg, SecondLeg>(); }
    static const std::string& ChannelNameLatex() { return __Channel_names_latex.EnumToString(ChannelId()); }

    static EventCategorySet DetermineEventCategories(EventInfo& event)
    {
        static const std::map<DiscriminatorWP, double> btag_working_points = {
            { DiscriminatorWP::Loose, cuts::btag_2016::CSVv2L },
            { DiscriminatorWP::Medium, cuts::btag_2016::CSVv2M }
        };

        EventCategorySet categories;
        categories.insert(EventCategory::Inclusive());

        const bool is_boosted = event.SelectFatJet(cuts::hh_bbtautau_2016::fatJetID::mass,
                                                   cuts::hh_bbtautau_2016::fatJetID::deltaR_subjet) != nullptr;

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
            for(const auto& wp_entry : btag_working_points) {
                categories.emplace(2, bjet_counts[wp_entry.first], wp_entry.first);
                categories.emplace(2, bjet_counts[wp_entry.first], wp_entry.first, is_boosted);
            }
        }
        return categories;
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

        MvaReaderSetupCollection mva_setup_collection;
        MvaReaderSetupEntryReader mva_entry_reader(mva_setup_collection);
        config_reader.AddEntryReader("MVA", mva_entry_reader, false);

        SampleDescriptorConfigEntryReader sample_entry_reader(sample_descriptors);
        config_reader.AddEntryReader("SAMPLE", sample_entry_reader, true);

        CombinedSampleDescriptorConfigEntryReader combined_entry_reader(cmb_sample_descriptors, sample_descriptors);
        config_reader.AddEntryReader("SAMPLE_CMB", combined_entry_reader, false);

        config_reader.ReadConfig(args.sources());

        if(!ana_setup_collection.count(args.setup()))
            throw exception("Setup '%1%' not found in the configuration file '%2%'.") % args.setup() % args.sources();
        ana_setup = ana_setup_collection.at(args.setup());

        if(args.mva_setup().size()) {
            if(!mva_setup_collection.count(args.mva_setup()))
                throw exception("MVA setup '%1%' not found.") % args.mva_setup();
            mva_setup = mva_setup_collection.at(args.mva_setup());
        }
        RemoveUnusedSamples();

        CreateEventSubCategoriesToProcess();
        InitializeMvaReader();
    }

    virtual ~BaseEventAnalyzer() {}

    void Run()
    {
        ProcessSamples(ana_setup.signals, "signal");
        ProcessSamples(ana_setup.data, "data");
        ProcessSamples(ana_setup.backgrounds, "background");
        ProcessCombinedSamples(ana_setup.cmb_samples);

        std::cout << "Producing inputs for limits..." << std::endl;
        LimitsInputProducer<FirstLeg, SecondLeg> limitsInputProducer(anaDataCollection, sample_descriptors,
                                                                     cmb_sample_descriptors);
        for(const std::string& hist_name : ana_setup.final_variables) {
            for(const auto& subCategory : EventSubCategoriesToProcess()) {
                std::cout << '\t' << hist_name << "/" << subCategory << std::endl;
                limitsInputProducer.Produce(args.output(), hist_name, ana_setup.limit_categories, subCategory,
                                            ana_setup.energy_scales, EventRegionsToProcess());
            }
        }

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
    virtual const EventCategorySet& EventCategoriesToProcess() const
    {
        static const EventCategorySet categories = {
            EventCategory::TwoJets_Inclusive(), EventCategory::TwoJets_ZeroBtag_Resolved(),
            EventCategory::TwoJets_OneBtag_Resolved(), /*EventCategory::TwoJets_OneLooseBtag(),*/
            EventCategory::TwoJets_TwoBtag_Resolved(), /*EventCategory::TwoJets_TwoLooseBtag()*/
            EventCategory::TwoJets_TwoLooseBtag_Boosted()
        };
        return categories;
    }

    virtual const EventSubCategorySet& EventSubCategoriesToProcess() const { return sub_categories_to_process; }

    virtual const EventRegionSet& EventRegionsToProcess() const
    {
        static const EventRegionSet regions = {
            EventRegion::OS_Isolated(), EventRegion::OS_AntiIsolated(),
            EventRegion::SS_Isolated(), EventRegion::SS_AntiIsolated(),
            EventRegion::SS_LooseIsolated()
        };
        return regions;
    }

    void CreateEventSubCategoriesToProcess()
    {
        const auto base = EventSubCategory().SetCutResult(SelectionCut::mh, true);
        sub_categories_to_process.insert(base);
        if(mva_setup.is_initialized()) {
            for(const auto& mva_sel : mva_setup->selections) {
                auto param = mva_sel.second;
                std::cout<<"Selection cut: "<<mva_sel.first<<" - name: "<<param.name
                        <<" spin: "<<param.spin<<" mass: "<<param.mass<<" cut: "<<param.cut<<std::endl;
                auto sub_category = EventSubCategory(base).SetCutResult(mva_sel.first, true);
                sub_categories_to_process.insert(sub_category);

            }
        }
    }

    void InitializeMvaReader()
    {
        if(!mva_setup.is_initialized()) return;
        for(const auto& method : mva_setup->trainings) {
            const auto& name = method.first;
            const auto& file = method.second;
            const auto& vars = mva_setup->variables.at(name);
            const auto& masses = mva_setup->masses.at(name);
            const auto& mass_range_pair= std::minmax_element(masses.begin(), masses.end());
            const Range<int> mass_range(static_cast<int>(*mass_range_pair.first),
                                        static_cast<int>(*mass_range_pair.second));
            const bool legacy = mva_setup->legacy.count(name);
            const bool legacy_lm = legacy && mva_setup->legacy.at(name) == "lm";
            mva_reader.AddRange(mass_range, name, file, vars, legacy, legacy_lm);
        }
    }

    virtual EventSubCategory DetermineEventSubCategory(EventInfo& event, const EventCategory& category,
                                                       std::map<SelectionCut, double>& mva_scores)
    {
        using namespace cuts::hh_bbtautau_2016::hh_tag;
        using MvaKey = std::tuple<std::string, int, int>;

        EventSubCategory sub_category;
        sub_category.SetCutResult(SelectionCut::mh,
                                  IsInsideMassWindow(event.GetHiggsTT(true).GetMomentum().mass(),
                                                     event.GetHiggsBB().GetMomentum().mass(),
                                                     category.HasBoostConstraint() && category.IsBoosted()));
        sub_category.SetCutResult(SelectionCut::KinematicFitConverged,
                                  event.GetKinFitResults().HasValidMass());

        if(mva_setup.is_initialized()) {

            std::map<MvaKey, double> scores;
            for(const auto& mva_sel : mva_setup->selections) {
                const auto& params = mva_sel.second;
                const MvaKey key{params.name, static_cast<int>(params.mass), params.spin};
                if(!scores.count(key))
                    scores[key] = mva_reader.Evaluate(*event, static_cast<int>(params.mass), params.name, params.spin);
                const double score = scores.at(key);
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
            if(sample.channels.size() && !sample.channels.count(ChannelId())) continue;
            std::cout << '\t' << sample.name << std::endl;
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
                processed_files.insert(sample_wp.file_path);
            }
        }
    }

    void ProcessDataSource(const SampleDescriptor& sample, const SampleDescriptor::Point& sample_wp,
                           std::shared_ptr<ntuple::EventTuple> tuple, const ntuple::ProdSummary& prod_summary)
    {
        const bool is_signal = ana_setup.IsSignal(sample.name);
        const bool need_to_blind = args.event_set() && (sample.sampleType == SampleType::TT || is_signal);
        const unsigned event_set = args.event_set(), half_split = prod_summary.n_splits / 2;
        const SummaryInfo summary(prod_summary);
        Event prevFullEvent, *prevFullEventPtr = nullptr;
        for(auto tupleEvent : *tuple) {
            if(ntuple::EventLoader::Load(tupleEvent, prevFullEventPtr).IsFull()) {
                prevFullEvent = tupleEvent;
                prevFullEventPtr = &prevFullEvent;
            }
            if(need_to_blind){
                if((event_set == 1 && tupleEvent.split_id >= half_split) || tupleEvent.split_id < half_split)
                    continue;
                tupleEvent.weight_total *= 2;
            }
            EventInfo event(tupleEvent, ntuple::JetPair{0, 1}, &summary);
            if(!ana_setup.energy_scales.count(event.GetEnergyScale())) continue;

            const auto eventCategories = DetermineEventCategories(event);
            for(auto eventCategory : eventCategories) {
                if (!EventCategoriesToProcess().count(eventCategory)) continue;
                const EventRegion eventRegion = DetermineEventRegion(event, eventCategory);
                for(const auto& region : EventRegionsToProcess()){
                    if(!eventRegion.Implies(region)) continue;

                    std::map<SelectionCut, double> mva_scores;
                    const auto eventSubCategory = DetermineEventSubCategory(event, eventCategory, mva_scores);
                    for(const auto& subCategory : EventSubCategoriesToProcess()) {
                        if(!eventSubCategory.Implies(subCategory)) continue;
                        SelectionCut mva_cut;
                        double mva_score = 0;
                        if(subCategory.TryGetLastMvaCut(mva_cut))
                            mva_score = mva_scores.at(mva_cut);
                        event.SetMvaScore(mva_score);
                        const EventAnalyzerDataId anaDataId(eventCategory, subCategory, region,
                                                            event.GetEnergyScale(), sample_wp.full_name);
                        if(sample.sampleType == SampleType::Data) {
                            ProcessDataEvent(anaDataId, event);
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
    }

    virtual void EstimateQCD(const SampleDescriptor& qcd_sample)
    {
        static const EventRegionSet sidebandRegions = {
            EventRegion::OS_AntiIsolated(), EventRegion::SS_Isolated(), EventRegion::SS_AntiIsolated(),
            EventRegion::SS_LooseIsolated()
        };
        static const EventEnergyScaleSet qcdEnergyScales = { EventEnergyScale::Central };

        for(const EventAnalyzerDataId& metaDataId : EventAnalyzerDataId::MetaLoop(EventCategoriesToProcess(),
                EventSubCategoriesToProcess(), sidebandRegions, qcdEnergyScales)) {
            const auto qcdAnaDataId = metaDataId.Set(qcd_sample.name);
            auto& qcdAnaData = anaDataCollection.Get(qcdAnaDataId);
            for(const auto& sample_name : sample_descriptors) {
                const SampleDescriptor& sample =  sample_name.second;
                if(sample.sampleType == SampleType::QCD) continue;
                if(ana_setup.IsSignal(sample.name)) continue;
                double factor = sample.sampleType == SampleType::Data ? +1 : -1;
                for(const auto& sample_wp : sample.working_points) {
                    const auto anaDataId = metaDataId.Set(sample_wp.full_name);
                    auto& anaData = anaDataCollection.Get(anaDataId);
                    for(const auto& sub_entry : anaData.template GetEntriesEx<TH1D>()) {
                        auto& entry = qcdAnaData.template GetEntryEx<TH1D>(sub_entry.first);
                        for(const auto& hist : sub_entry.second->GetHistograms()) {
                            entry(hist.first).Add(hist.second.get(), factor);
                        }
                    }
                }
            }
        }

        for(const EventAnalyzerDataId& metaDataId : EventAnalyzerDataId::MetaLoop(EventCategoriesToProcess(),
                EventSubCategoriesToProcess(), qcdEnergyScales)) {
            const auto anaDataId = metaDataId.Set(qcd_sample.name);
            auto& osIsoData = anaDataCollection.Get(anaDataId.Set(EventRegion::OS_Isolated()));
            auto& ssIsoData = anaDataCollection.Get(anaDataId.Set(EventRegion::SS_Isolated()));
            auto& ssLooseIsoData = anaDataCollection.Get(anaDataId.Set(EventRegion::SS_LooseIsolated()));
            auto& osAntiIsoData = anaDataCollection.Get(anaDataId.Set(EventRegion::OS_AntiIsolated()));
            auto& ssAntiIsoData = anaDataCollection.Get(anaDataId.Set(EventRegion::OS_AntiIsolated()));
            for(const auto& sub_entry : ssIsoData.template GetEntriesEx<TH1D>()) {
                auto& entry_osIso = osIsoData.template GetEntryEx<TH1D>(sub_entry.first);
                auto& entry_ss_looseIso = ssLooseIsoData.template GetEntryEx<TH1D>(sub_entry.first);
                auto& entry_osAntiIso = osAntiIsoData.template GetEntryEx<TH1D>(sub_entry.first);
                auto& entry_ssAntiIso = ssAntiIsoData.template GetEntryEx<TH1D>(sub_entry.first);
                for(const auto& hist : sub_entry.second->GetHistograms()) {
                    std::string debug_info, negative_bins_info;
                    if(!HasNegativeContribution(*hist.second,debug_info, negative_bins_info)) continue;
                    const auto osAntiIso_integral = analysis::Integral(entry_osAntiIso(hist.first), true);
                    const auto ssAntiIso_integral = analysis::Integral(entry_ssAntiIso(hist.first), true);
                    if (osAntiIso_integral.GetValue() <= 0){
                        std::cout << "Warning: OS Anti Iso integral less or equal 0 for " << hist.first << std::endl;
                        continue;
                    }

                    if (ssAntiIso_integral.GetValue() <= 0){
                        std::cout << "Warning: SS Anti Iso integral less or equal 0 for " << hist.first << std::endl;
                        continue;
                    }
                    const double k_factor = osAntiIso_integral.GetValue()/ssAntiIso_integral.GetValue();
                    const auto ssIso_integral = analysis::Integral(*hist.second, true);
                    if (ssIso_integral.GetValue() <= 0){
                        std::cout << "Warning: SS Iso integral less or equal 0 for " << hist.first << std::endl;
                        continue;
                    }
                    const double total_yield = ssIso_integral.GetValue() * k_factor;
                    entry_osIso(hist.first).CopyContent(entry_ss_looseIso(hist.first));
                    analysis::RenormalizeHistogram(entry_osIso(hist.first),total_yield,true);


                }
            }
        }
    }

    bool HasNegativeContribution(/*const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                      EventRegion eventRegion,*/ TH1D& histogram, /*const std::string& current_category,*/
                                      std::string& debug_info, std::string& negative_bins_info)
    {
        static const double correction_factor = 0.0000001;

        std::ostringstream ss_debug;

//        ss_debug << "\nSubtracting background for '" << histogram.GetName() << "' in region " << eventRegion
//                 << " for Event category '" << anaDataMetaId.eventCategory
//                 << "' for data category '" << current_category
//                 << "'.\nInitial integral: " << Integral(histogram, true) << ".\n";


        const PhysicalValue original_Integral = Integral(histogram, true);
        ss_debug << "\nSubtracted hist for '" << histogram.GetName() << ".\n";
        ss_debug << "Integral after bkg subtraction: " << original_Integral << ".\n";
        debug_info = ss_debug.str();
        if (original_Integral.GetValue() < 0) {
            std::cout << debug_info << std::endl;
            std::cout << "Integral after bkg subtraction is negative for histogram '"
                << histogram.GetName() << /*"' in event category " << anaDataMetaId.eventCategory
                << " for event region " << eventRegion << "." <<*/ std::endl;
            return false;
        }

        std::ostringstream ss_negative;

        for (Int_t n = 1; n <= histogram.GetNbinsX(); ++n) {
            if (histogram.GetBinContent(n) >= 0) continue;
            const std::string prefix = histogram.GetBinContent(n) + histogram.GetBinError(n) >= 0 ? "WARNING" : "ERROR";

            ss_negative << prefix << ": " << histogram.GetName() << " Bin " << n << ", content = "
                        << histogram.GetBinContent(n) << ", error = " << histogram.GetBinError(n)
                        << ", bin limits=[" << histogram.GetBinLowEdge(n) << "," << histogram.GetBinLowEdge(n+1)
                        << "].\n";
            const double error = correction_factor - histogram.GetBinContent(n);
            const double new_error = std::sqrt(std::pow(error,2) + std::pow(histogram.GetBinError(n),2));
            histogram.SetBinContent(n, correction_factor);
            histogram.SetBinError(n, new_error);
        }
        analysis::RenormalizeHistogram(histogram, original_Integral.GetValue(), true);
        negative_bins_info = ss_negative.str();
        return true;
    }

    void ProcessDataEvent(const EventAnalyzerDataId& anaDataId, EventInfo& event)
    {
        for(const auto& es : ana_setup.energy_scales) {
            anaDataCollection.Fill(anaDataId.Set<EventEnergyScale>(es), event, 1);
        }
    }

    virtual void ProcessSpecialEvent(const SampleDescriptor& sample, const SampleDescriptor::Point& /*sample_wp*/,
                                     const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight)
    {
        if(sample.sampleType == SampleType::DY) {
            static auto const find_b_index = [&]() {
                const auto& param_names = sample.GetModelParameterNames();
                const auto b_param_iter = param_names.find("b");
                if(b_param_iter == param_names.end())
                    throw exception("Unable to find b_parton WP for DY sample");
                return b_param_iter->second;
            };
            static const size_t b_index = find_b_index();

            bool wp_found = false;
            static constexpr double pt_cut =18, b_Flavour = 5;

            for(const auto& sample_wp : sample.working_points) {
                const size_t n_b_partons = static_cast<size_t>(sample_wp.param_values.at(b_index));
                size_t n_genJets = 0;
                for(const auto& b_Candidates : event.GetHiggsBB().GetDaughterMomentums()) {
                    for(size_t i=0; i<event->genJets_p4.size(); i++){
                        const auto& jet_p4 = event->genJets_p4.at(i);
                        const auto& jet_hadronFlavour = event->genJets_hadronFlavour.at(i);
                        double deltaR = ROOT::Math::VectorUtil::DeltaR(b_Candidates, jet_p4);
                        if (jet_p4.Pt() <= pt_cut || jet_hadronFlavour != b_Flavour || deltaR >= 0.3) continue;
                        n_genJets++;
                    }
                }
                if(n_genJets == n_b_partons ||
                        (n_b_partons == sample.GetNWorkingPoints() - 1
                         && n_genJets > n_b_partons)) {
                    const auto finalId = anaDataId.Set(sample_wp.full_name);
                    anaDataCollection.Fill(finalId, event, weight * sample_wp.norm_sf);
                    wp_found = true;
                    break;
                }
            }
            if(!wp_found)
                throw exception("Unable to find WP for DY event with lhe_n_b_partons = %1%") % event->lhe_n_b_partons;

        } else if(sample.sampleType == SampleType::TT) {
            anaDataCollection.Fill(anaDataId, event, weight / event->weight_top_pt);
            if(anaDataId.Get<EventEnergyScale>() == EventEnergyScale::Central) {
//                const double weight_topPt = event->weight_total * sample.cross_section * ana_setup.int_lumi
//                        / event.GetSummaryInfo()->totalShapeWeight_withTopPt;
                anaDataCollection.Fill(anaDataId.Set(EventEnergyScale::TopPtUp), event, weight); // FIXME
                anaDataCollection.Fill(anaDataId.Set(EventEnergyScale::TopPtDown), event, weight); // FIXME
            }
        } else
            throw exception("Unsupported special event type '%1%'.") % sample.sampleType;
    }

    void ProcessCombinedSamples(const std::vector<std::string>& sample_names)
    {
        std::cout << "Processing combined samples... " << std::endl;
        for(const std::string& sample_name : sample_names) {
            if(!cmb_sample_descriptors.count(sample_name))
                throw exception("Combined sample '%1%' not found.") % sample_name;
            CombinedSampleDescriptor& sample = cmb_sample_descriptors.at(sample_name);
            if(sample.channels.size() && !sample.channels.count(ChannelId())) continue;
            std::cout << sample.name << std::endl;
            for(const std::string& sub_sample_name : sample.sample_descriptors) {
                if(!sample_descriptors.count(sub_sample_name))
                    throw exception("Unable to create '%1%': sub-sample '%2%' not found.")
                        % sample_name % sub_sample_name;
                SampleDescriptor& sub_sample =  sample_descriptors.at(sub_sample_name);
                AddSampleToCombined(sample, sub_sample);
            }
        }
    }

    void AddSampleToCombined(CombinedSampleDescriptor& sample, SampleDescriptor& sub_sample)
    {
        for(const EventAnalyzerDataId& metaDataId : EventAnalyzerDataId::MetaLoop(EventCategoriesToProcess(),
                EventSubCategoriesToProcess(), EventRegionsToProcess(), ana_setup.energy_scales)) {
            const auto anaDataId = metaDataId.Set(sample.name);
            auto& anaData = anaDataCollection.Get(anaDataId);
            for(const auto& sub_sample_wp : sub_sample.working_points) {
                const auto subDataId = metaDataId.Set(sub_sample_wp.full_name);
                auto& subAnaData = anaDataCollection.Get(subDataId);
                for(const auto& sub_entry : subAnaData.template GetEntriesEx<TH1D>()) {
                    auto& entry = anaData.template GetEntryEx<TH1D>(sub_entry.first);
                    for(const auto& hist : sub_entry.second->GetHistograms()) {
                        entry(hist.first).Add(hist.second.get(), 1);
                    }
                }
            }
        }
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
    boost::optional<MvaReaderSetup> mva_setup;
    SampleDescriptorCollection sample_descriptors;
    CombinedSampleDescriptorCollection cmb_sample_descriptors;
    EventSubCategorySet sub_categories_to_process;
    mva_study::MvaReader mva_reader;
};

} // namespace analysis
