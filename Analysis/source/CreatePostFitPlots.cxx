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

        std::map<analysis::EventCategory,size_t> categories_map = {{ analysis::EventCategory::Parse("2j1bR_noVBF"), 0 },
                                                                   { analysis::EventCategory::Parse("2j2b+R_noVBF")  , 1 },
                                                                   { analysis::EventCategory::Parse("2j1b+B_noVBF"), 2}};

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
                                                            UncertaintySource::None, UncertaintyScale::Central,
                                                            item.full_name);
                        if(item.datacard_name.empty()) continue;

                        if(ana_setup.IsSignal(sample.name) && item.full_name != args.signal_name())
                            continue;

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

    }

//private:



private:
    AnalyzerArguments args;
    std::set<std::string> activeVariables;
    std::shared_ptr<TFile> outputFile;
    PropertyConfigReader histConfig;
};

} // namespace analysis

PROGRAM_MAIN(analysis::CreatePostFitPlots, analysis::AnalyzerArguments)
