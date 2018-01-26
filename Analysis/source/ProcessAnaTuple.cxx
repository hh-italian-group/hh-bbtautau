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
    OPT_ARG(size_t, n_parallel, 10);
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
        histConfig.Parse(FullPath(ana_setup.hist_cfg));
        if(!ana_setup.unc_cfg.empty()) {
            ConfigReader config_reader;
            unc_collection = std::make_shared<ModellingUncertaintyCollection>();
            ModellingUncertaintyEntryReader unc_reader(*unc_collection);
            config_reader.AddEntryReader("UNC", unc_reader, true);
            config_reader.ReadConfig(FullPath(ana_setup.unc_cfg));
        }
    }

    void Run()
    {
        const std::set<std::string> signal_names(ana_setup.signals.begin(), ana_setup.signals.end());
        const std::set<std::string> bkg_names(ana_setup.backgrounds.begin(), ana_setup.backgrounds.end());
        const auto samplesToDraw = PlotsProducer::CreateOrderedSampleCollection(
                    ana_setup.draw_sequence, sample_descriptors, cmb_sample_descriptors, ana_setup.signals,
                    ana_setup.data, args.channel());

        std::ofstream qcd_out(args.output() +"_QCD.txt");

        const std::vector<EventSubCategory> all_subCategories(sub_categories_to_process.begin(),
                                                              sub_categories_to_process.end());

        for(size_t n = 0; n * args.n_parallel() < all_subCategories.size(); ++n) {
            AnaDataCollection anaDataCollection(outputFile, channelId, &tupleReader.GetAnaTuple(), activeVariables,
                                                histConfig.GetItems(), bkg_names, unc_collection);
            EventSubCategorySet subCategories;
            for(size_t k = 0; k < args.n_parallel() && n * args.n_parallel() + k < all_subCategories.size(); ++k) {
                const auto& subCategory = all_subCategories.at(n * args.n_parallel() + k);
                subCategories.insert(subCategory);
                std::cout << subCategory << " ";
            }
            std::cout << std::endl;

            std::cout << "\tCreating histograms..." << std::endl;
            ProduceHistograms(anaDataCollection, subCategories);

            std::cout << "\tProcessing combined samples and QCD... " << std::endl;
            for(const auto& subCategory : subCategories) {
                ProcessCombinedSamples(anaDataCollection, subCategory, ana_setup.cmb_samples);
                for(const auto& sample : sample_descriptors) {
                    if(sample.second.sampleType == SampleType::QCD) {
                        EstimateQCD(anaDataCollection, subCategory, sample.second, qcd_out);
                        break;
                    }
                }
            }

            if(args.shapes()) {
                std::cout << "\t\tProducing inputs for limits..." << std::endl;
                LimitsInputProducer limitsInputProducer(anaDataCollection, sample_descriptors,
                                                        cmb_sample_descriptors);
                for(const std::string& hist_name : ana_setup.final_variables) {
                    if(!activeVariables.count(hist_name)) continue;
                    for(const auto& subCategory : subCategories)
                        limitsInputProducer.Produce(args.output(), hist_name, ana_setup.limit_categories, subCategory,
                                                    ana_setup.energy_scales, ana_setup.regions, mva_sel_aliases);
                }
            }

            if(args.draw()) {
                std::cout << "\t\tCreating plots..." << std::endl;
                PlotsProducer plotsProducer(anaDataCollection, samplesToDraw, FullPath(ana_setup.plot_cfg),
                                            ana_setup.plot_page_opt);
                std::string pdf_prefix = args.output();
                if(n != 0)
                    pdf_prefix += "_part" + ToString(n + 1);
                plotsProducer.PrintStackedPlots(pdf_prefix, EventRegion::SignalRegion(), ana_setup.categories,
                                                subCategories, signal_names);
            }
        }

        std::cout << "Saving output file..." << std::endl;
    }

private:

    void ProduceHistograms(AnaDataCollection& anaDataCollection, const EventSubCategorySet& subCategories)
    {
        auto& tuple = tupleReader.GetAnaTuple();
        for(Long64_t entryIndex = 0; entryIndex < tuple.GetEntries(); ++entryIndex) {
            tuple.GetEntry(entryIndex);
            for(size_t n = 0; n < tuple().dataIds.size(); ++n) {
                const auto& dataId = tupleReader.GetDataIdByIndex(n);
                if(!subCategories.count(dataId.Get<EventSubCategory>())) continue;
                tupleReader.UpdateSecondaryBranches(dataId, n);
                anaDataCollection.Fill(dataId, tuple().weight);
            }
        }
    }

    void EstimateQCD(AnaDataCollection& anaDataCollection, const EventSubCategory& subCategory,
                     const SampleDescriptor& qcd_sample, std::ostream& log)
    {
        static const EventRegionSet sidebandRegions = {
            EventRegion::OS_AntiIsolated(), EventRegion::SS_Isolated(), EventRegion::SS_AntiIsolated(),
            EventRegion::SS_LooseIsolated()
        };
        static const EventEnergyScaleSet qcdEnergyScales = { EventEnergyScale::Central };
        const EventSubCategorySet subCategories = { subCategory };

        for(const EventAnalyzerDataId& metaDataId : EventAnalyzerDataId::MetaLoop(ana_setup.categories,
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

        for(const EventAnalyzerDataId& metaDataId : EventAnalyzerDataId::MetaLoop(ana_setup.categories,
                subCategories, qcdEnergyScales)) {
            const auto anaDataId = metaDataId.Set(qcd_sample.name);
            auto& osIsoData = anaDataCollection.Get(anaDataId.Set(EventRegion::OS_Isolated()));
            auto& ssIsoData = anaDataCollection.Get(anaDataId.Set(EventRegion::SS_Isolated()));
            auto& ssLooseIsoData = anaDataCollection.Get(anaDataId.Set(EventRegion::SS_LooseIsolated()));
            auto& osAntiIsoData = anaDataCollection.Get(anaDataId.Set(EventRegion::OS_AntiIsolated()));
            auto& ssAntiIsoData = anaDataCollection.Get(anaDataId.Set(EventRegion::SS_AntiIsolated()));
            for(const auto& sub_entry : ssIsoData.template GetEntriesEx<TH1D>()) {
                auto& entry_osIso = osIsoData.template GetEntryEx<TH1D>(sub_entry.first);
                auto& entry_ss_looseIso = ssLooseIsoData.template GetEntryEx<TH1D>(sub_entry.first);
                auto& entry_osAntiIso = osAntiIsoData.template GetEntryEx<TH1D>(sub_entry.first);
                auto& entry_ssAntiIso = ssAntiIsoData.template GetEntryEx<TH1D>(sub_entry.first);
                for(const auto& hist : sub_entry.second->GetHistograms()) {
                    log << anaDataId << ": " << sub_entry.first << " " << hist.first << "\n";
                    std::string debug_info, negative_bins_info;
                    if(!FixNegativeContributions(*hist.second,debug_info, negative_bins_info)) {
                        log << debug_info << "\n" << negative_bins_info << "\n";
                        continue;
                    }
                    const auto osAntiIso_integral = Integral(entry_osAntiIso(hist.first), true);
                    const auto ssAntiIso_integral = Integral(entry_ssAntiIso(hist.first), true);
                    if (osAntiIso_integral.GetValue() <= 0 || osAntiIso_integral.IsCompatible(PhysicalValue::Zero)){
                        log << "Warning: OS Anti Iso integral is too small " << hist.first << std::endl;
                        continue;
                    }

                    if (ssAntiIso_integral.GetValue() <= 0 || ssAntiIso_integral.IsCompatible(PhysicalValue::Zero)){
                        log << "Warning: SS Anti Iso integral is too small " << hist.first << std::endl;
                        continue;
                    }
                    const auto k_factor = osAntiIso_integral / ssAntiIso_integral;
                    const auto ssIso_integral = analysis::Integral(*hist.second, true);
                    if (ssIso_integral.GetValue() <= 0){
                        log << "Warning: SS Iso integral less or equal 0 for " << hist.first << std::endl;
                        continue;
                    }
                    const auto total_yield = ssIso_integral * k_factor;
                    log << anaDataId << ": osAntiIso integral = " << osAntiIso_integral
                        << ", ssAntiIso integral = " << ssAntiIso_integral << ", os/ss sf = " << k_factor
                        << ", ssIso integral = " << ssIso_integral << ", total yield = " << total_yield << std::endl;
                    entry_osIso(hist.first).CopyContent(entry_ss_looseIso(hist.first));
                    analysis::RenormalizeHistogram(entry_osIso(hist.first), total_yield.GetValue(), true);
                }
            }
        }
        log << std::endl;
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
        for(const EventAnalyzerDataId& metaDataId : EventAnalyzerDataId::MetaLoop(ana_setup.categories,
                ana_setup.regions, ana_setup.energy_scales)) {
            const auto anaDataId = metaDataId.Set(sample.name).Set(subCategory);
            auto& anaData = anaDataCollection.Get(anaDataId);
            for(const auto& sub_sample_wp : sub_sample.working_points) {
                const auto subDataId = metaDataId.Set(sub_sample_wp.full_name).Set(subCategory);
                auto& subAnaData = anaDataCollection.Get(subDataId);
                for(const auto& sub_entry : subAnaData.GetEntriesEx<TH1D>()) {
                    auto& entry = anaData.GetEntryEx<TH1D>(sub_entry.first);
                    for(const auto& hist : sub_entry.second->GetHistograms()) {
                        entry(hist.first).AddHistogram(*hist.second);
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
    std::shared_ptr<ModellingUncertaintyCollection> unc_collection;
};

} // namespace analysis

PROGRAM_MAIN(analysis::ProcessAnaTuple, analysis::AnalyzerArguments)
