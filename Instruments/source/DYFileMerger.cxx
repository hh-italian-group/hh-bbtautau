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


//std::string tree_name{"tree_name", "Tree on which we work"};
struct Arguments {
    REQ_ARG(std::string, tree_name);
    REQ_ARG(std::string, cfg_name);
    REQ_ARG(std::string, file_cfg_name);
    REQ_ARG(std::string, output_file);
};

namespace analysis {

class DYFileMerger { // simple analyzer definition
public:
    using GenMap = ntuple::GenEventCountMap;
    using VSD = std::vector<SampleDescriptor<DYBinDescriptor, ntuple::GenEventCountMap>>;
    using VBD = std::vector<DYBinDescriptor>;
    DYFileMerger(const Arguments& _args) : args(_args)
    {
        // Analyzer initialization (e.g. open input/output files, parse configs...)
        LoadInputs();
        LoadBins();
    }


private:
    SampleDescriptor<DYBinDescriptor, ntuple::GenEventCountMap> global_map;
    SampleDescriptor<DYBinDescriptor, ntuple::GenEventCountMap> inclusive;
    VSD all_samples;
    VBD output_bins;


    void LoadInputs()
    {
        // Fill inclusive and all_samples
        // analyzer code
        analysis::ConfigReader config_reader;

        analysis::DYBinDescriptorCollection file_descriptors;
        analysis::DYFileConfigEntryReader file_entry_reader(file_descriptors);
        config_reader.AddEntryReader("FILE", file_entry_reader, true);

        config_reader.ReadConfig(args.file_cfg_name());

        for (auto file_descriptor : file_descriptors){ //loop on DYJets files
            const DYBinDescriptor file_descriptor_element = file_descriptor.second;
            std::cout << "File descriptor characteristics: " << file_descriptor.first << ", " <<
                         file_descriptor_element.file_path << ", " << file_descriptor_element.fileType
                      << std::endl;
            SampleDescriptor<DYBinDescriptor, ntuple::GenEventCountMap> sample_desc;
            sample_desc.bin = file_descriptor_element;
            //global_map.bin = file_descriptor_element;
            if (file_descriptor_element.fileType == analysis::FileType::inclusive)
                inclusive.bin = file_descriptor_element;

            auto inputFile = root_ext::OpenRootFile(file_descriptor_element.file_path);
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
                    if (file_descriptor_element.fileType == analysis::FileType::inclusive)
                        inclusive.gen_counts[genId] += nevents;

                } // end loop on gen event info
            } //end loop on entries
            all_samples.push_back(sample_desc);
        } //end loop n file_descriptors
    }

    void LoadBins()
    {
       output_bins = DYBinDescriptor::LoadConfig(args.cfg_name());
    }

    void CalculateWeight(DYBinDescriptor& output_bin) const
    {
        size_t all_events = global_map.Integral(output_bin);
        for(unsigned i = 0; i < all_samples.size(); ++i) {
            const SampleDescriptor<DYBinDescriptor, ntuple::GenEventCountMap> sample = all_samples.at(i);
            size_t contribution = sample.Integral(output_bin);
            if(!contribution) continue;
            //formula 2
            PhysicalValue nu ( contribution / double(sample.Integral()), sqrt(double(contribution))/double(sample.Integral()));
            PhysicalValue weight ((double)nu.GetValue()/(double)all_events, ((double)(all_events - contribution)/(double)std::pow(all_events,2))*sqrt(double(contribution))/double(sample.Integral()));

            if(!(sample.bin.fileType == analysis::FileType::inclusive)) {
                size_t sample_contribution = inclusive.Integral(sample.bin);
                // formula 3
                PhysicalValue nu_incl(sample_contribution/double(inclusive.Integral()),
                                      sqrt(double(sample_contribution))/double(inclusive.Integral()));
                nu *= nu_incl;
                weight *= nu_incl;
            }

            if(output_bin.nu.GetStatisticalError() > nu.GetStatisticalError()) {
                output_bin.nu = nu;
                output_bin.ref_sample = sample.bin.name;
                output_bin.weight = weight;
                output_bin.fileType = sample.bin.fileType;
            }

        }

        if(output_bin.nu.GetStatisticalError() == std::numeric_limits<double>::infinity())
            throw exception("ref not found");
    }

    //old definition
    //output_bin.weight = output_bin.nu.GetValue() / global_map.Integral(output_bin);

    public:


    void Run()
    {
        for(unsigned i = 0; i < output_bins.size(); ++i)
        {
            DYBinDescriptor& output_bin = output_bins.at(i);
            CalculateWeight(output_bin);
        }
        DYBinDescriptor::SaveCfg(args.output_file(), output_bins);

    }
private:
    Arguments args;
};

}

PROGRAM_MAIN(analysis::DYFileMerger, Arguments) // definition of the main program function
