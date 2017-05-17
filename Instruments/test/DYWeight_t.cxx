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
#include <iostream>


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
    using VectorSampleDescriptor = std::vector<SampleDescriptor<DYBinDescriptor, ntuple::GenEventCountMap>>;
    using VectorDYBinDescriptor = std::vector<DYBinDescriptor>;
    using Range_weight_map = std::map<int, double>;
    using DoubleRange_map = std::map<int, Range_weight_map>;
    using Event = ntuple::SummaryTuple;
    using DY_weightPtr = std::shared_ptr<mc_corrections::DY_weight>;

    DYWeight_t(const Arguments& _args) : args(_args)
    {
        dy_weight = DY_weightPtr(new class mc_corrections::DY_weight("hh-bbtautau/Analysis/config/dyjets_weights.cfg"));
        std::vector<analysis::sample_merging::DYBinDescriptor> dy_descriptors =
                analysis::sample_merging::DYBinDescriptor::LoadConfig(args.cfg_name());
        std::cout << "Loaded descriptor" << std::endl;

        for (unsigned n = 0; n < dy_descriptors.size(); ++n){
            const analysis::sample_merging::DYBinDescriptor dybin_descriptor = dy_descriptors.at(n);
            const double weight = dybin_descriptor.weight.GetValue()/dybin_descriptor.inclusive_integral;
            std::cout << "weight before map: " << weight << std::endl;
            for(int n_jet = dybin_descriptor.n_jet.min(); n_jet <= dybin_descriptor.n_jet.max(); ++n_jet) {
//                if(!dy_weight_map.count(n_jet)) {
                    dy_weight_map[n_jet] = DoubleRange_map();
//                }
                DoubleRange_map& njet_map = dy_weight_map.at(n_jet);
                for(int n_bjet = dybin_descriptor.n_bjet.min(); n_bjet <= dybin_descriptor.n_bjet.max(); ++n_bjet) {
//                    if(!njet_map.count(n_bjet))
                        njet_map[n_bjet] = Range_weight_map();
                    Range_weight_map& nbjet_map = njet_map.at(n_bjet);
                    for(int ht = dybin_descriptor.n_ht.min(); ht <= dybin_descriptor.n_ht.max(); ++ht) {
                        if(nbjet_map.count(ht))
                            throw exception("Repeated bin");
                        nbjet_map[ht] = weight;
                    }
                }
            }
        }
        LoadInputs();
    }

public:

    void Run()
    {

    }


private:

    Arguments args;
    std::map<int, DoubleRange_map> dy_weight_map;
    DY_weightPtr dy_weight;


    void LoadInputs()
    {
        std::setprecision(std::numeric_limits<double>::digits10 + 1);
        analysis::ConfigReader config_reader;

        DYBinDescriptorCollection file_descriptors;
        DYFileConfigEntryReader file_entry_reader(file_descriptors);
        config_reader.AddEntryReader("FILE", file_entry_reader, true);

        config_reader.ReadConfig(args.file_cfg_name());

        double totalWeight = 0;
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
                        double weight = dy_weight->GetWeight(n_partons,n_b_partons,ht);
                        std::cout.unsetf ( std::ios::floatfield );
                        std::cout << weight << ": weight" <<  std::endl;
                        double weight_prime = nevents * weight;
                        std::cout << weight_prime << ": weight_prime" <<std::endl;
                        totalWeight += weight_prime;

                    } // end loop on gen event info
                } //end loop on entries

            } // end loop on files

        } //end loop n file_descriptors
        std::cout << totalWeight << ": Total weight" << std::endl;

    }


};

} //namespace sample_merging

} //namespace analysis

PROGRAM_MAIN(analysis::sample_merging::DYWeight_t, Arguments)
