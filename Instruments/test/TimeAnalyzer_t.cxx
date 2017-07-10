/*! Merge DYJets files.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include <iostream>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "Instruments/include/TimeFileConfigEntryReader.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "Instruments/include/TimeFileDescriptor.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"



struct Arguments {
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> input_path{"input_path", "Input path of the samples"};
    run::Argument<std::string> cfg_name{"cfg_name", "file path cfg"};
    run::Argument<std::string> output_file{"output_file", "Output root file"};
    //run::Argument<std::string> output_cfg_name{"output_cfg_name", "output file cfg"};
};

namespace analysis {

class TimeAnalyzerData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    TH1D_ENTRY(m_sv, 200, 0, 400)

};


class TimeAnalyzer_t {
public:

    TimeAnalyzer_t(const Arguments& _args) : args(_args),
        anaData(output) , output(root_ext::CreateRootFile(args.output_file()))
    {
        std::cout << "Starting..." << std::endl;
        LoadInputs();
    }

    void Run() {}


private:
    Arguments args;
    TimeAnalyzerData anaData;
    std::shared_ptr<TFile> output;

    void LoadInputs()
    {
        analysis::ConfigReader config_reader;

        TimeFileDescriptorCollection file_descriptors;
        TimeFileConfigEntryReader file_entry_reader(file_descriptors);
        config_reader.AddEntryReader("FILE", file_entry_reader, true);

        config_reader.ReadConfig(args.cfg_name());

        for (auto file_descriptor : file_descriptors){ //loop on files
            const TimeFileDescriptor file_descriptor_element = file_descriptor.second;

            for (auto single_file_path : file_descriptor_element.file_paths){ //loop on files


                std::cout << "File descriptor characteristics: " << file_descriptor.first << ", " <<
                             single_file_path << std::endl;
                std::string filename = args.input_path()  + "/" + single_file_path;
                auto inputFile = root_ext::OpenRootFile(filename);
                ntuple::SummaryTuple summaryTuple(args.tree_name(), inputFile.get(), true);
                const Long64_t n_entries = summaryTuple.GetEntries();

                for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on entries
                    summaryTuple.GetEntry(current_entry);

                    std::cout << "Nprocessed events" << summaryTuple.data().numberOfProcessedEvents << std::endl;

                } //end loop on entries

            } // end loop on files

        } //end loop n file_descriptors

    }


};

} //namespace analysis

PROGRAM_MAIN(analysis::TimeAnalyzer_t, Arguments)
