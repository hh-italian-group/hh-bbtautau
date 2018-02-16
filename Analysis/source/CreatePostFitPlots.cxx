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
    REQ_ARG(std::string, signal_name);
    REQ_ARG(std::string, output);
    OPT_ARG(std::string, var, "");
};

class CreatePostFitPlots : public EventAnalyzerCore {
public:
    using AnaData = ::analysis::EventAnalyzerData;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection;
    using PlotsProducer = ::analysis::StackedPlotsProducer;

    CreatePostFitPlots(const AnalyzerArguments& _args) :
        EventAnalyzerCore(_args, _args.channel()), args(_args), activeVariables({args.var()}),
        outputFile(root_ext::CreateRootFile(args.output() + "_postfit.root"))
    {
        histConfig.Parse(FullPath(ana_setup.hist_cfg));
    }

    void Run()
    {
        const std::set<std::string> signal_names(ana_setup.signals.begin(), ana_setup.signals.end());
        const std::set<std::string> bkg_names(ana_setup.backgrounds.begin(), ana_setup.backgrounds.end());
        const std::set<std::string> data_names(ana_setup.data.begin(), ana_setup.data.end());
        std::set<std::string> all_samples;
        all_samples.insert(signal_names.begin(),signal_names.end());
        all_samples.insert(bkg_names.begin(),bkg_names.end());
        all_samples.insert(data_names.begin(),data_names.end());

        const auto samplesToDraw = PlotsProducer::CreateOrderedSampleCollection(
                    ana_setup.draw_sequence, sample_descriptors, cmb_sample_descriptors, ana_setup.signals,
                    ana_setup.data, args.channel());

        std::map<analysis::EventCategory,size_t> categories_map = {{ analysis::EventCategory::TwoJets_OneBtag_Resolved(), 0 },
                                                                   { analysis::EventCategory::TwoJets_TwoBtag_Resolved(), 1 },
                                                                   { analysis::EventCategory::TwoJets_TwoLooseBtag_Boosted(), 2}};

        auto inputFile = root_ext::OpenRootFile(args.input());
        AnaDataCollection anaDataCollection(outputFile, channelId, nullptr, activeVariables,
                                            histConfig.GetItems(), all_samples, nullptr);

        for(const auto& category : ana_setup.categories){
            if(!categories_map.count(category))
                continue;
            size_t n = categories_map.at(category);
            for(const auto& sample_name : all_samples){
                SampleDescriptor& sample = sample_descriptors.at(sample_name);
                for(const SampleDescriptorBase::Point& item : sample.working_points){
                    for (const auto& subcategory : sub_categories_to_process){

                        const EventAnalyzerDataId anaDataId(category, subcategory, analysis::EventRegion::OS_Isolated(),
                                                            analysis::EventEnergyScale::Central, item.full_name);
                        if(item.datacard_name.empty()) continue;

                        if(ana_setup.IsSignal(sample.name) && item.full_name != args.signal_name()) {
//                            std::cout << "Item full name: " << item.full_name << ", args: " << args.signal_name() << std::endl;
                            continue;
                        }

                        const std::string hist_dir_name = ana_setup.IsSignal(sample.name) ?
                                "hh_ttbb_"+ToString(args.channel())+"_"+ToString(n)+"_13TeV_postfit/"+sample.postfit_name :
                                "hh_ttbb_"+ToString(args.channel())+"_"+ToString(n)+"_13TeV_postfit/"+item.datacard_name ;

                        auto hist = std::shared_ptr<TH1>
                                (root_ext::TryReadObject<TH1>(*inputFile,hist_dir_name));
                        if (!hist){
                            std::cout << "WARNING! Datacard Name doesn't correspond to sample name: " << anaDataId <<
                                         ", hist dir name: "<< hist_dir_name <<  std::endl;
                            continue;
                        }
                        std::cout << "Loading " << anaDataId << " ..." << std::endl;
                        auto& histAnaData = anaDataCollection.Get(anaDataId);
                        auto& histAnaData_entry = histAnaData.GetEntryEx<TH1D>(args.var());
                        histAnaData_entry().CopyContent(*hist);

                    }
                }
            }
        }

        std::cout << std::endl;
        std::cout << "\tProcessing combined samples " << std::endl;

        for(const auto& subcategory : sub_categories_to_process) {
            ProcessCombinedSamples(anaDataCollection, subcategory, ana_setup.cmb_samples);
        }

        std::cout << "\t\tCreating plots..." << std::endl;
        PlotsProducer plotsProducer(anaDataCollection, samplesToDraw, FullPath(ana_setup.plot_cfg),
                                    ana_setup.plot_page_opt);
        std::string pdf_prefix = args.output();
        plotsProducer.PrintStackedPlots(pdf_prefix, EventRegion::SignalRegion(), ana_setup.categories,
                                        sub_categories_to_process, signal_names, &sample_descriptors.at("TotalBkg"));



        std::cout << "Saving output file..." << std::endl;
    }

private:

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
    std::shared_ptr<TFile> outputFile;
    PropertyConfigReader histConfig;
};

} // namespace analysis

PROGRAM_MAIN(analysis::CreatePostFitPlots, analysis::AnalyzerArguments)
