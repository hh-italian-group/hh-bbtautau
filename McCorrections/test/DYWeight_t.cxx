/*! Merge DYJets files.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "hh-bbtautau/McCorrections/include/NJets_HT_BinFileConfigEntryReader.h"
#include "hh-bbtautau/McCorrections/include/NJets_HT_weight.h"

struct Arguments {
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> input_path{"input_path", "Input path of the samples"};
    run::Argument<std::string> cfg_name{"cfg_name", "cfg bin splitting"};
    run::Argument<std::string> file_cfg_name{"file_cfg_name", "DY file cfg"};
};

namespace analysis {

namespace sample_merging{

class DYWeight_t {
public:
    using GenMap = ntuple::GenEventCountMap;
    using VectorSampleDescriptor = std::vector<SampleDescriptor<NJets_HT_BinFileDescriptor, ntuple::GenEventCountMap>>;
    using VectorDYBinDescriptor = std::vector<NJets_HT_BinFileDescriptor>;
    using Event = ntuple::SummaryTuple;
    using DY_weight = mc_corrections::NJets_HT_weight;

    DYWeight_t(const Arguments& _args) : args(_args), dy_weight("DY", args.cfg_name())
    {
        std::cout << "Ciao" << std::endl;
    }

public:

    void Run()
    {
        std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        analysis::ConfigReader config_reader;

        NJets_HT_BinDescriptorCollection file_descriptors;
        NJets_HT_BinFileConfigEntryReader file_entry_reader(file_descriptors);
        config_reader.AddEntryReader("FILE", file_entry_reader, true);

        config_reader.ReadConfig(args.file_cfg_name());

        double totalWeight = 0;
        for (auto file_descriptor : file_descriptors){ //loop on DYJets files
            const NJets_HT_BinFileDescriptor file_descriptor_element = file_descriptor.second;

            for (auto single_file_path : file_descriptor_element.file_paths){ //loop on files


                std::cout << "File descriptor characteristics: " << file_descriptor.first << ", " <<
                             single_file_path << ", " << file_descriptor_element.fileType
                          << std::endl;
                std::string filename = args.input_path()  + "/" + single_file_path;
                auto inputFile = root_ext::OpenRootFile(filename);
                ntuple::ExpressTuple summaryTuple(args.tree_name(), inputFile.get(), true);
                const Long64_t n_entries = summaryTuple.GetEntries();

                size_t nevents = 0;
                double weight = 0;
                for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on entries
                    summaryTuple.GetEntry(current_entry);

                    ++nevents;
                    UInt_t n_partons = summaryTuple.data().lhe_n_partons;
                    UInt_t n_b_partons =  summaryTuple.data().lhe_n_b_partons;
                    UInt_t ht =  summaryTuple.data().lhe_ht10_bin;
                    weight += dy_weight.GetWeight(n_partons,n_b_partons,ht);


                } //end loop on entries
                //double weight_prime = nevents * weight;
                totalWeight += weight;
            } // end loop on files

        } //end loop n file_descriptors
        std::cout << totalWeight << ": Total weight" << std::endl;
        if (std::abs(totalWeight - 1) > 1e-3)
            throw analysis::exception("Total Weight not compatible with 1!");
    }


private:
    Arguments args;
    DY_weight dy_weight;

};

} //namespace sample_merging

} //namespace analysis

PROGRAM_MAIN(analysis::sample_merging::DYWeight_t, Arguments)
