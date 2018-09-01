/*! Analyzer
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */


#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include <string>
#include <iostream>
#include <unordered_map>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
};

namespace analysis {


namespace fs = boost::filesystem;
class ReadingTimeAnalyzer {
public:

    ReadingTimeAnalyzer(const Arguments& _args): args(_args) {}

    void Run()
    {
            for (const auto & entry : fs::directory_iterator(args.input_path())){
                auto file = entry.path().string();

                std::size_t pos = file.find(".");
                std::string partial = file.substr(pos,pos+5);
                std::cout<<partial<<std::endl;
                if (partial != ".root") continue;
                auto inputFile = root_ext::OpenRootFile(file);
                auto tuple = ntuple::CreateEventTuple("muMu", inputFile.get(), true, ntuple::TreeState::Full);
                std::cout << entry << " number of events: " << tuple->size() << std::endl;
            }

    }

private:
    Arguments args;
};

}

PROGRAM_MAIN(analysis::ReadingTimeAnalyzer, Arguments) // definition of the main program function
