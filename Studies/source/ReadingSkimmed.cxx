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
#include <TSystemDirectory.h>
#include <TSystemFile.h>


struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, output_file);
};

namespace analysis {


class ReadingTimeAnalyzer{
public:

    ReadingTimeAnalyzer(const Arguments& _args): args(_args),
        outfile(root_ext::CreateRootFile(args.output_file())) {}

    void Run()
    {
        using Event = ntuple::Event;
        using EventTuple = ntuple::EventTuple;

        namespace fs = boost::filesystem;
        for (const auto & entry : fs::directory_iterator(args.input_path())){
            auto file = entry.path().string();
            std::cout<<file<<std::endl;
            std::size_t pos = file.find(".");
            std::string partial = file.substr(pos,pos+5);
            std::cout<<partial<<std::endl;
            if (partial != ".root") continue;

            auto inputFile = root_ext::OpenRootFile(file);

            std::string tree;
            TList *dirlist = inputFile->GetListOfKeys();
            TIterator *iter = dirlist->MakeIterator();
            TObject *key = iter->Next();
            while (key) {
               std::string avoid=key->GetName();
               if (avoid == "summary") {
                   key = iter->Next();
                continue;
               }
               tree = key->GetName();
               std::cout<<tree<<std::endl;
               key = iter->Next();
            }

            std::cout<<tree<<std::endl;

            auto tuple = ntuple::CreateEventTuple(tree, inputFile.get(), true, ntuple::TreeState::Skimmed);

            auto histo = std::make_shared<TH1D>((file+tree).c_str(), (file+tree).c_str(), 25,0,200);
            histo->SetCanExtend(TH1::kAllAxes);

            for(const Event& event : *tuple) {
                histo->Fill(event.p4_1.px());
            }
            auto directory = root_ext::GetDirectory(*outfile, file+tree);
            root_ext::WriteObject(*histo, directory);
            std::cout << entry << " number of events: " << tuple->size() << std::endl;

        }
    }

private:
    Arguments args;
    std::shared_ptr<TFile> outfile;
};

}

PROGRAM_MAIN(analysis::ReadingTimeAnalyzer, Arguments) // definition of the main program function
