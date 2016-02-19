/*! Analyzer which recalls three ananlysis to produce flatTree for Htautau using looser selection.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "h-tautau/Analysis/include/TreeExtractor.h"
#include "Analysis/include/Config.h"

#include "FlatTreeProducer_etau.cxx"
#include "FlatTreeProducer_mutau.cxx"
#include "FlatTreeProducer_tautau.cxx"

class FlatTreeProducer {
public:
    FlatTreeProducer(const std::string& inputFileName, const std::string& outputMuTauFile,
                     const std::string& outputETauFile, const std::string& outputTauTauFile,
                     const std::string& configFileName, const std::string& _prefix = "none",
                     size_t _maxNumberOfEvents = 0)
        : config(configFileName), timer(config.ReportInterval(), std::cout), maxNumberOfEvents(_maxNumberOfEvents),
          treeExtractor(_prefix == "none" ? "" : _prefix, inputFileName, config.extractMCtruth(), config.MaxTreeVersion()),
          HmutauAnalyzer(inputFileName, outputMuTauFile, configFileName, "external", _maxNumberOfEvents),
          HetauAnalyzer(inputFileName, outputETauFile, configFileName, "external", _maxNumberOfEvents),
          HtautauAnalyzer(inputFileName, outputTauTauFile, configFileName, "external", _maxNumberOfEvents)
    { }

    virtual void Run()
    {
        size_t n = 0;
        auto _event = std::shared_ptr<analysis::EventDescriptor>(new analysis::EventDescriptor());
        for(; ( !maxNumberOfEvents || n < maxNumberOfEvents ) && treeExtractor.ExtractNext(*_event); ++n) {
            timer.Report(n);
//            std::cout << _event->eventId().eventId << std::endl;
            if(config.RunSingleEvent() && _event->eventId().eventId != config.SingleEventId()) continue;
            HmutauAnalyzer.ProcessEventWithEnergyUncertainties(_event);
            HetauAnalyzer.ProcessEventWithEnergyUncertainties(_event);
            HtautauAnalyzer.ProcessEventWithEnergyUncertainties(_event);
            if(config.RunSingleEvent()) break;
        }
        timer.Report(n, true);
    }


private:
    analysis::Config config;
    analysis::tools::ProgressReporter timer;
    size_t maxNumberOfEvents;
    std::shared_ptr<const analysis::EventDescriptor> event;
    analysis::TreeExtractor treeExtractor;
    FlatTreeProducer_mutau HmutauAnalyzer;
    FlatTreeProducer_etau HetauAnalyzer;
    FlatTreeProducer_tautau HtautauAnalyzer;
    double eventWeight;
};

namespace make_tools {
template<>
struct Factory<FlatTreeProducer> {
    static FlatTreeProducer* Make(int argc, char *argv[])
    {
        std::cout << "Command line: ";
        for(int n = 0; n < argc; ++n)
            std::cout << argv[n] << " ";
        std::cout << std::endl;
        if(argc < 6 || argc > 8)
            throw std::runtime_error("Invalid number of command line arguments.");

        int n = 0;
        const std::string inputFileName = argv[++n];
        const std::string outputMuTauFile = argv[++n];
        const std::string outputETauFile = argv[++n];
        const std::string outputTauTauFile = argv[++n];
        const std::string configFileName = argv[++n];
        if(argc <= n)
            return new FlatTreeProducer(inputFileName, outputMuTauFile, outputETauFile, outputTauTauFile,
                                        configFileName);

        const std::string prefix = argv[++n];
        if(argc <= n)
            return new FlatTreeProducer(inputFileName, outputMuTauFile, outputETauFile, outputTauTauFile,
                                        configFileName, prefix);

        char c;
        size_t maxNumberOfEvents;
        std::istringstream ss_nEvents(argv[++n]);
        ss_nEvents >> c >> maxNumberOfEvents;
        if(c != '@')
            throw std::runtime_error("Bad command line format.");

        return new FlatTreeProducer(inputFileName, outputMuTauFile, outputETauFile, outputTauTauFile,
                                    configFileName, prefix, maxNumberOfEvents);
    }
};
} // make_tools
