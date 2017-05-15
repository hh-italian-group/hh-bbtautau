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


struct Arguments {
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> input_path{"input_path", "Input path of the samples"};
    run::Argument<std::string> cfg_name{"cfg_name", "cfg bin splitting"};
    run::Argument<std::string> file_cfg_name{"file_cfg_name", "DY file cfg"};
};

namespace analysis {

namespace sample_merging{

class DYCheckWeight {
public:
    using GenMap = ntuple::GenEventCountMap;
    using VectorSampleDescriptor = std::vector<SampleDescriptor<DYBinDescriptor, ntuple::GenEventCountMap>>;
    using VectorDYBinDescriptor = std::vector<DYBinDescriptor>;
    using Range_weight_map = std::map<Range<int>, double>;
    using DoubleRange_map = std::map<Range<int>, Range_weight_map>;
    using Event = ntuple::SummaryTuple;

    DYCheckWeight(const Arguments& _args) : args(_args)
    {
        std::vector<analysis::sample_merging::DYBinDescriptor> dy_descriptors =
                analysis::sample_merging::DYBinDescriptor::LoadConfig(args.cfg_name());
        for (unsigned n = 0; n < dy_descriptors.size(); ++n){
            const analysis::sample_merging::DYBinDescriptor dybin_descriptor = dy_descriptors.at(n);
            dy_weight_map[dybin_descriptor.n_jet][dybin_descriptor.n_bjet][dybin_descriptor.n_ht]
                    = dybin_descriptor.weight.GetValue()/dybin_descriptor.inclusive_integral;
        }
        LoadInputs();
    }

public:

    void Run()
    {

    }


private:

    Arguments args;
    std::map<Range<int>, DoubleRange_map> dy_weight_map;


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
                        UInt_t n_partons = summaryTuple.data().lhe_n_partons.at(i);
                        UInt_t n_b_partons =  summaryTuple.data().lhe_n_b_partons.at(i);
                        UInt_t ht =  summaryTuple.data().lhe_ht10_bin.at(i);
                        double weight = GetWeight(n_partons,n_b_partons,ht);
                        double weight_prime = nevents * weight;
                        totalWeight += weight_prime;

                    } // end loop on gen event info
                } //end loop on entries

            } // end loop on files

        } //end loop n file_descriptors
        std::cout << "Total weight: " << totalWeight << std::endl;

    }

    double GetWeight(Int_t n_partons, Int_t n_b_partons, Int_t ht)
    {
        for (auto iter : dy_weight_map){
            Range<int> n_jet = iter.first;
            if (!(n_jet.Contains(n_partons))) continue;
            DoubleRange_map doubleRange_map = iter.second;
            for (auto iter_1 : doubleRange_map){
                Range<int> n_bjet = iter_1.first;
                if (!(n_bjet.Contains(n_b_partons))) continue;
                Range_weight_map range_weight_map = iter_1.second;
                for (auto iter_2 : range_weight_map){
                    Range<int> n_ht = iter_2.first;
                    if (n_ht.Contains(ht))
                        return iter_2.second;
                }
            }
        }
        throw analysis::exception("drell-yan merge weight not found.");
    }


};

} //namespace sample_merging

} //namespace analysis

PROGRAM_MAIN(analysis::sample_merging::DYCheckWeight, Arguments)
