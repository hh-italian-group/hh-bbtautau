/*! Merge DYJets files.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "Instruments/include/NJets_HT_BinFileConfigEntryReader.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "Instruments/include/SampleDescriptor.h"


struct Arguments {
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> input_path{"input_path", "Input path of the samples"};
    run::Argument<std::string> cfg_name{"cfg_name", "cfg bin for DY or Wjets splitting"};
    run::Argument<std::string> file_cfg_name{"file_cfg_name", "DY or Wjets file cfg"};
    run::Argument<std::string> output_file{"output_file", "Output file"};
};

namespace analysis {

namespace sample_merging{

class NJets_HT_BinFileMerger {
public:
    using GenMap = ntuple::GenEventCountMap;
    using VectorSampleDescriptor = std::vector<SampleDescriptor<NJets_HT_BinFileDescriptor, ntuple::GenEventCountMap>>;
    using VectorDYBinDescriptor = std::vector<NJets_HT_BinFileDescriptor>;
    NJets_HT_BinFileMerger(const Arguments& _args) : args(_args)
    {
        LoadInputs();
        std::cout << "Done LoadInputs" << std::endl;
        output_bins = NJets_HT_BinFileDescriptor::LoadConfig(args.cfg_name());
        std::cout << "Done LoadCfg" << std::endl;
    }

public:

    void Run()
    {
        for(auto& output_bin : output_bins)
        {
            CalculateWeight(output_bin);
        }
        std::cout << "Done CalculateWeight" << std::endl;
        NJets_HT_BinFileDescriptor::SaveCfg(args.output_file(), output_bins);
        std::cout << "Done SaveCfg" << std::endl;
    }


private:
    SampleDescriptor<NJets_HT_BinFileDescriptor, ntuple::GenEventCountMap> global_map;
    SampleDescriptor<NJets_HT_BinFileDescriptor, ntuple::GenEventCountMap> inclusive;
    VectorSampleDescriptor all_samples;
    VectorDYBinDescriptor output_bins;
    Arguments args;


    void LoadInputs()
    {
        analysis::ConfigReader config_reader;

        NJets_HT_BinDescriptorCollection file_descriptors;
        NJets_HT_BinFileConfigEntryReader file_entry_reader(file_descriptors);
        config_reader.AddEntryReader("FILE", file_entry_reader, true);

        config_reader.ReadConfig(args.file_cfg_name());

        for (auto file_descriptor : file_descriptors){ //loop on DYJets files
            const NJets_HT_BinFileDescriptor file_descriptor_element = file_descriptor.second;

            SampleDescriptor<NJets_HT_BinFileDescriptor, ntuple::GenEventCountMap> sample_desc;
            sample_desc.bin = file_descriptor_element;
            //global_map.bin = file_descriptor_element;
            if (file_descriptor_element.fileType == FileType::inclusive)
                inclusive.bin = file_descriptor_element;

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
                        const ntuple::GenId genId(summaryTuple.data().lhe_n_partons.at(i),
                                                  summaryTuple.data().lhe_n_b_partons.at(i),
                                                  summaryTuple.data().lhe_ht10_bin.at(i));
                        size_t nevents = summaryTuple.data().lhe_n_events.at(i);
                        sample_desc.gen_counts[genId] += nevents;
                        global_map.gen_counts[genId] += nevents;
                        if (file_descriptor_element.fileType == FileType::inclusive)
                            inclusive.gen_counts[genId] += nevents;

                    } // end loop on gen event info
                } //end loop on entries

            } // end loop on files
            all_samples.push_back(sample_desc);
        } //end loop n file_descriptors
    }


    void CalculateWeight(NJets_HT_BinFileDescriptor& output_bin) const
    {
        double all_events = global_map.Integral(output_bin);
        for(auto& sample : all_samples) {
            double contribution = sample.Integral(output_bin);
            if(!contribution) continue;
            //formula 2
            PhysicalValue nu ( contribution , sqrt(contribution));
            PhysicalValue weight (nu.GetValue()/all_events, (all_events - contribution)/std::pow(all_events,2)*sqrt(contribution));

            if(!(sample.bin.fileType == FileType::inclusive)) {
                double sample_contribution = inclusive.Integral(sample.bin);
                // formula 3
                PhysicalValue nu_incl(sample_contribution/ sample.Integral(),
                                      sqrt(sample_contribution)/ sample.Integral());

                nu *= nu_incl;
                weight *= nu_incl;
            }

            if(output_bin.nu.GetStatisticalError() > nu.GetStatisticalError()) {
                output_bin.nu = nu;
                output_bin.ref_sample = sample.bin.name;
                output_bin.weight = weight;
                output_bin.fileType = sample.bin.fileType;
                output_bin.inclusive_integral = inclusive.Integral();
            }

        }

        if(output_bin.nu.GetStatisticalError() == std::numeric_limits<double>::infinity())
            throw exception("ref not found");
    }



};

} //namespace sample_merging

} //namespace analysis

PROGRAM_MAIN(analysis::sample_merging::NJets_HT_BinFileMerger, Arguments)
