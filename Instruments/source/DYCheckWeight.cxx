/*! Merge DYJets files.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "Instruments/include/DYFileConfigEntryReader.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "Instruments/include/SampleDescriptor.h"
#include "McCorrections/include/DY_weight.h"
#include "McCorrections/include/EventWeights_HH.h"


struct Arguments {
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> input_path{"input_path", "Input path of the samples"};
    run::Argument<std::string> cfg_name{"cfg_name", "cfg bin splitting"};
    run::Argument<std::string> file_cfg_name{"file_cfg_name", "DY file cfg"};
    run::Argument<std::string> output_file{"output_file", "Output file"};
};

namespace analysis {

namespace sample_merging{

class DYFileMerger {
public:
    using GenMap = ntuple::GenEventCountMap;
    using VectorSampleDescriptor = std::vector<SampleDescriptor<DYBinDescriptor, ntuple::GenEventCountMap>>;
    using VectorDYBinDescriptor = std::vector<DYBinDescriptor>;
    DYFileMerger(const Arguments& _args) : args(_args)
    {
        LoadInputs();
    }

public:

    void Run()
    {

    }


private:

    Arguments args;


    void LoadInputs()
    {
        analysis::ConfigReader config_reader;

        DYBinDescriptorCollection file_descriptors;
        DYFileConfigEntryReader file_entry_reader(file_descriptors);
        config_reader.AddEntryReader("FILE", file_entry_reader, true);

        config_reader.ReadConfig(args.file_cfg_name());

        Int_t totalWeight = 0;
        for (auto file_descriptor : file_descriptors){ //loop on DYJets files
            const DYBinDescriptor file_descriptor_element = file_descriptor.second;

            for (auto single_file_path : file_descriptor_element.file_paths){ //loop on files


                std::cout << "File descriptor characteristics: " << file_descriptor.first << ", " <<
                             single_file_path << ", " << file_descriptor_element.fileType
                          << std::endl;
                std::string filename = args.input_path()  + "/" + single_file_path;
                auto inputFile = root_ext::OpenRootFile(filename);
                ntuple::SummaryTuple summaryTuple(args.tree_name(), inputFile.get(), true);
                const Long64_t n_entries = summaryTuple.GetEntries();

                for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on entries
                    summaryTuple.GetEntry(current_entry);
                    for (size_t i = 0; i < summaryTuple.data().lhe_n_partons.size(); ++i){ // loop on gen event info

                        size_t nevents = summaryTuple.data().lhe_n_events.at(i);
                        totalWeight +=

                    } // end loop on gen event info
                } //end loop on entries

            } // end loop on files

        } //end loop n file_descriptors
    }


};

} //namespace sample_merging

} //namespace analysis

PROGRAM_MAIN(analysis::sample_merging::DYFileMerger, Arguments)
