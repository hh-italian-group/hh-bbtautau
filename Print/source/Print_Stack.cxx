/*! Print stack with specified name superimposing several files.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "hh-bbtautau/Analysis/include/EventAnalyzerCore.h"
#include "hh-bbtautau/Analysis/include/StackedPlotsProducer.h"

namespace analysis {

struct Arguments : CoreAnalyzerArguments {
    REQ_ARG(std::string, input);
    REQ_ARG(std::string, output);
    REQ_ARG(Channel, channel);
    REQ_ARG(std::string, vars);
};

class PrintStack : public EventAnalyzerCore {
public:
    using AnaData = ::analysis::EventAnalyzerData;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection;
    using PlotsProducer = ::analysis::StackedPlotsProducer;

    PrintStack(const Arguments& _args) :
        EventAnalyzerCore(_args, _args.channel()), args(_args), inputFile(root_ext::OpenRootFile(args.input())),
        vars(ParseVarSet(args.vars()))
    {
        histConfig.Parse(ana_setup.hist_cfg);
        anaDataCollection = std::make_shared<AnaDataCollection>(inputFile, args.channel(), nullptr, vars,
                                                                histConfig.GetItems());
    }

    void Run()
    {
        std::cout << "Creating plots..." << std::endl;
        const std::set<std::string> signal_names(ana_setup.signals.begin(), ana_setup.signals.end());
        const auto samplesToDraw = PlotsProducer::CreateOrderedSampleCollection(
                    ana_setup.draw_sequence, sample_descriptors, cmb_sample_descriptors, ana_setup.signals,
                    ana_setup.data, channelId);

        PlotsProducer plotsProducer(*anaDataCollection, samplesToDraw, ana_setup.plot_cfg, ana_setup.plot_page_opt,
                                    vars);
        plotsProducer.PrintStackedPlots(args.output(), EventRegion::SignalRegion(), ana_setup.categories,
                                        sub_categories_to_process, signal_names);
    }

private:
    static std::set<std::string> ParseVarSet(const std::string& active_vars_str)
    {
        const auto list = SplitValueList(active_vars_str, false, ", \t", true);
        return std::set<std::string>(list.begin(), list.end());
    }

private:
    Arguments args;
    std::shared_ptr<TFile> inputFile;
    std::set<std::string> vars;
    PropertyConfigReader histConfig;
    std::shared_ptr<AnaDataCollection> anaDataCollection;
};

} // namespace analysis

PROGRAM_MAIN(analysis::PrintStack, analysis::Arguments)
