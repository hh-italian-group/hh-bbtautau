/*! Post-processing of the analysis results.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/EventAnalyzerCore.h"
#include "hh-bbtautau/Analysis/include/AnaTuple.h"
#include "hh-bbtautau/Analysis/include/EventAnalyzerDataCollection.h"
#include "hh-bbtautau/Analysis/include/StackedPlotsProducer.h"
#include "hh-bbtautau/Analysis/include/LimitsInputProducer.h"

namespace analysis {

struct AnalyzerArguments : CoreAnalyzerArguments {
    REQ_ARG(Channel, channel);
    REQ_ARG(std::string, input);
    REQ_ARG(std::string, output);
    OPT_ARG(bool, shapes, true);
    OPT_ARG(bool, draw, true);
    OPT_ARG(std::string, vars, "");
};

class ProcessAnaTuple : public EventAnalyzerCore {
public:
    using AnaData = ::analysis::EventAnalyzerData;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection;
    using PlotsProducer = ::analysis::StackedPlotsProducer;

    ProcessAnaTuple(const AnalyzerArguments& _args) :
        EventAnalyzerCore(_args, _args.channel()), args(_args), activeVariables(ParseVarSet(args.vars())),
        tupleReader(args.input(), args.channel(), activeVariables),
        outputFile(root_ext::CreateRootFile(args.output() + "_full.root"))
    {

        histConfig.Parse(ana_setup.hist_cfg);
    }

    void Run()
    {
        const std::set<std::string> signal_names(ana_setup.signals.begin(), ana_setup.signals.end());
        const auto samplesToDraw = PlotsProducer::CreateOrderedSampleCollection(
                    ana_setup.draw_sequence, sample_descriptors, cmb_sample_descriptors, ana_setup.signals,
                    ana_setup.data, args.channel());

        for(const auto& subCategory : EventSubCategoriesToProcess()) {
            std::cout << subCategory << std::endl;
            AnaDataCollection anaDataCollection(outputFile, channelId, &tupleReader.GetAnaTuple(), activeVariables,
                                                histConfig.GetItems());
            std::cout << "\tCreating histograms..." << std::endl;
            ProduceHistograms(anaDataCollection, subCategory);
            std::cout << "\tProcessing combined samples... " << std::endl;
            ProcessCombinedSamples(anaDataCollection, subCategory, ana_setup.cmb_samples);
            for(const auto& sample : sample_descriptors) {
                if(sample.second.sampleType == SampleType::QCD) {
                    std::cout << "\tEstimating QCD..." << std::endl;
                    EstimateQCD(anaDataCollection, subCategory, sample.second);
                    break;
                }
//                std::cout << "\t\t\t" << hist.first << " -> " << entry.Name() << ", "
//                          << subDataId << ", " << anaDataId << std::endl;

            }

            if(args.shapes()) {
                std::cout << "\tProducing inputs for limits..." << std::endl;
                LimitsInputProducer limitsInputProducer(anaDataCollection, sample_descriptors, cmb_sample_descriptors);
                for(const std::string& hist_name : ana_setup.final_variables) {
                    if(!activeVariables.count(hist_name)) continue;
                    limitsInputProducer.Produce(args.output(), hist_name, ana_setup.limit_categories, subCategory,
                                                ana_setup.energy_scales, EventRegionsToProcess());
                }
            }

            if(args.draw()) {
                std::cout << "\tCreating plots..." << std::endl;
                PlotsProducer plotsProducer(anaDataCollection, samplesToDraw, false, true);
                const std::string pdf_prefix = args.output() + "_" + ToString(subCategory);
                plotsProducer.PrintStackedPlots(pdf_prefix, EventRegion::SignalRegion(), EventCategoriesToProcess(),
                                                {subCategory}, signal_names);
            }
        }

        std::cout << "Saving output file..." << std::endl;
    }

private:

    void ProduceHistograms(AnaDataCollection& anaDataCollection, const EventSubCategory& subCategory)
    {
        auto& tuple = tupleReader.GetAnaTuple();
        for(Long64_t entryIndex = 0; entryIndex < tuple.GetEntries(); ++entryIndex) {
            tuple.GetEntry(entryIndex);
            for(size_t n = 0; n < tuple().dataIds.size(); ++n) {
                const auto& dataId = tupleReader.GetDataIdByIndex(n);
                if(dataId.Get<EventSubCategory>() != subCategory) continue;
                tupleReader.UpdateSecondaryBranches(n);
                anaDataCollection.Fill(dataId, tuple().weight);
            }
        }
    }

    void EstimateQCD(AnaDataCollection& anaDataCollection, const EventSubCategory& subCategory,
                     const SampleDescriptor& qcd_sample)
    {
        static const EventRegionSet sidebandRegions = {
            EventRegion::OS_AntiIsolated(), EventRegion::SS_Isolated(), EventRegion::SS_AntiIsolated(),
            EventRegion::SS_LooseIsolated()
        };
        static const EventEnergyScaleSet qcdEnergyScales = { EventEnergyScale::Central };
        const EventSubCategorySet subCategories = { subCategory };

        for(const EventAnalyzerDataId& metaDataId : EventAnalyzerDataId::MetaLoop(EventCategoriesToProcess(),
                subCategories, sidebandRegions, qcdEnergyScales)) {
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
                subCategories, qcdEnergyScales)) {
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
                    if(!FixNegativeContributions(*hist.second,debug_info, negative_bins_info)) continue;
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

    void ProcessCombinedSamples(AnaDataCollection& anaDataCollection, const EventSubCategory& subCategory,
                                const std::vector<std::string>& sample_names)
    {

        for(const std::string& sample_name : sample_names) {
            if(!cmb_sample_descriptors.count(sample_name))
                throw exception("Combined sample '%1%' not found.") % sample_name;
            CombinedSampleDescriptor& sample = cmb_sample_descriptors.at(sample_name);
            if(sample.channels.size() && !sample.channels.count(channelId)) continue;
            std::cout << "\t\t" << sample.name << std::endl;
            for(const std::string& sub_sample_name : sample.sample_descriptors) {
                if(!sample_descriptors.count(sub_sample_name))
                    throw exception("Unable to create '%1%': sub-sample '%2%' not found.")
                        % sample_name % sub_sample_name;
                SampleDescriptor& sub_sample =  sample_descriptors.at(sub_sample_name);
                AddSampleToCombined(anaDataCollection, subCategory, sample, sub_sample);
            }
        }
    }

    void AddSampleToCombined(AnaDataCollection& anaDataCollection, const EventSubCategory& subCategory,
                             CombinedSampleDescriptor& sample, SampleDescriptor& sub_sample)
    {
        for(const EventAnalyzerDataId& metaDataId : EventAnalyzerDataId::MetaLoop(EventCategoriesToProcess(),
                EventRegionsToProcess(), ana_setup.energy_scales)) {
            const auto anaDataId = metaDataId.Set(sample.name).Set(subCategory);
            auto& anaData = anaDataCollection.Get(anaDataId);
            for(const auto& sub_sample_wp : sub_sample.working_points) {
                const auto subDataId = metaDataId.Set(sub_sample_wp.full_name).Set(subCategory);
                auto& subAnaData = anaDataCollection.Get(subDataId);
                for(const auto& sub_entry : subAnaData.GetEntriesEx<TH1D>()) {
                    auto& entry = anaData.GetEntryEx<TH1D>(sub_entry.first);
                    for(const auto& hist : sub_entry.second->GetHistograms()) {
                        entry(hist.first).Add(hist.second.get(), 1);
                    }
                }
            }
        }
    }

    static std::set<std::string> ParseVarSet(const std::string& active_vars_str)
    {
        const auto list = SplitValueList(active_vars_str, false, ", \t", true);
        return std::set<std::string>(list.begin(), list.end());
    }

private:
    AnalyzerArguments args;
    std::set<std::string> activeVariables;
    bbtautau::AnaTupleReader tupleReader;
    std::shared_ptr<TFile> outputFile;
    PropertyConfigReader histConfig;
};

} // namespace analysis

PROGRAM_MAIN(analysis::ProcessAnaTuple, analysis::AnalyzerArguments)
