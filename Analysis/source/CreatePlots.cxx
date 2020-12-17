/*! Post-processing of the analysis results.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Core/include/ProgressReporter.h"
#include "hh-bbtautau/Analysis/include/AnaTuple.h"
#include "hh-bbtautau/Analysis/include/EventAnalyzerCore.h"
#include "hh-bbtautau/Analysis/include/EventAnalyzerDataCollection.h"
#include "hh-bbtautau/Analysis/include/LimitsInputProducer.h"
#include "hh-bbtautau/Analysis/include/SampleDescriptorConfigEntryReader.h"
#include "hh-bbtautau/Analysis/include/StackedPlotsProducer.h"

namespace analysis {

struct AnalyzerArguments : CoreAnalyzerArguments {
    REQ_ARG(Channel, channel);
    REQ_ARG(std::string, input);
    REQ_ARG(std::string, output);
    REQ_ARG(std::string, vars);
    OPT_ARG(size_t, n_parallel, 10);
};

class CreatePlots : public EventAnalyzerCore {
public:
    using AnaData = ::analysis::EventAnalyzerData;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection;
    using PlotsProducer = ::analysis::StackedPlotsProducer;

    CreatePlots(const AnalyzerArguments& _args) :
        EventAnalyzerCore(_args, _args.channel(), false), args(_args), activeVariables(ParseVarSet(args.vars())),
        inputFile(root_ext::OpenRootFile(args.input() + "_full.root"))
    {
        for (auto& hist_config : ana_setup.hist_cfg)
            histConfig.Parse(FullPath(hist_config));
        //histConfig.Parse(FullPath(ana_setup.hist_cfg));
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

        const std::vector<EventSubCategory> all_subCategories(sub_categories_to_process.begin(),
                                                              sub_categories_to_process.end());

        for(size_t n = 0; n * args.n_parallel() < all_subCategories.size(); ++n) {
            AnaDataCollection anaDataCollection(inputFile, channelId, activeVariables, histConfig.GetItems(),
                                                true, bkg_names, unc_collection);
            EventSubCategorySet subCategories;
            for(size_t k = 0; k < args.n_parallel() && n * args.n_parallel() + k < all_subCategories.size(); ++k) {
                const auto& subCategory = all_subCategories.at(n * args.n_parallel() + k);
                subCategories.insert(subCategory);
                std::cout << subCategory << " ";
            }
            std::cout << std::endl;

            std::cout << "\t\tCreating plots..." << std::endl;
            PlotsProducer plotsProducer(anaDataCollection, samplesToDraw, FullPath(ana_setup.plot_cfg),
                                        ana_setup.plot_page_opt, activeVariables);
            std::string pdf_prefix = args.output();
            if(n != 0)
                pdf_prefix += "_part" + ToString(n + 1);
            plotsProducer.PrintStackedPlots(pdf_prefix, EventRegion::SignalRegion(), ana_setup.categories,
                                            subCategories, signal_names);
        }
    }

private:

    static std::set<std::string> ParseVarSet(const std::string& active_vars_str)
    {
        const auto list = SplitValueList(active_vars_str, false, ", \t", true);
        return std::set<std::string>(list.begin(), list.end());
    }

private:
    AnalyzerArguments args;
    std::set<std::string> activeVariables;
    std::shared_ptr<TFile> inputFile;
    PropertyConfigReader histConfig;
    std::shared_ptr<ModellingUncertaintyCollection> unc_collection;
};

} // namespace analysis

PROGRAM_MAIN(analysis::CreatePlots, analysis::AnalyzerArguments)
