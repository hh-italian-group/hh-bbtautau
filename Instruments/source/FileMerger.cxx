/*! Merge DYJets files.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h" // definition of wrappers for the program main and program arguments.
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "Instruments/include/FileConfigEntryReader.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"

struct Arguments { // list of all program arguments
    //REQ_ARG(std::string, input_file); // required argument "input_file"
    REQ_ARG(std::string, tree_name); // tree name
    REQ_ARG(std::string, cfg_name); // cfg name
    REQ_ARG(std::string, file_cfg_name); // DY cfg name
    REQ_ARG(std::string, output_file); // output file
};

namespace analysis {

struct SampleDescriptor {
    DYBinDescriptor bin;
    using GenMap = ntuple::GenEventCountMap;
    GenMap gen_counts;

    size_t Integral() const
    {
        size_t total_n_events = 0;
        for (auto genEventCount : gen_counts){
            const ntuple::GenId genId = genEventCount.first;
            const size_t nevents = genEventCount.second;
            total_n_events += nevents;
            std::cout << "GenId: " << genId.n_partons << " " << genId.n_b_partons << " " << genId.ht10_bin <<
                         ", nevents: " << nevents << std::endl;
        }
        std::cout << "Total events in the sumple: " << total_n_events << std::endl;
        return total_n_events;
    }
    size_t Integral(const DYBinDescriptor& binDescriptor_range) const
    {     
        size_t totalEvents = 0;
        for (auto genEventCount : gen_counts){
            const ntuple::GenId genId = genEventCount.first;
            const size_t nevents = genEventCount.second;
            if(binDescriptor_range.n_jet.Contains(genId.n_partons) &&
                  binDescriptor_range.n_bjet.Contains(genId.n_b_partons) &&
                  binDescriptor_range.n_ht.Contains(genId.ht10_bin))
                totalEvents += nevents;
            continue;

        }
        return totalEvents;
    }
};

class FileMerger { // simple analyzer definition
public:
    using GenMap = SampleDescriptor::GenMap;
    using VSD = std::vector<SampleDescriptor>;
    using VBD = std::vector<DYBinDescriptor>;
    FileMerger(const Arguments& _args) : args(_args)
    {
        // Analyzer initialization (e.g. open input/output files, parse configs...)
        LoadInputs();
        LoadBins();
    }


private:
    SampleDescriptor global_map;
    SampleDescriptor inclusive;
    VSD all_samples;
    VBD output_bins;


    void LoadInputs()
    {
        // Fill inclusive and all_samples
        // analyzer code
        analysis::ConfigReader config_reader;

        analysis::DYBinDescriptorCollection file_descriptors;
        analysis::FileConfigEntryReader file_entry_reader(file_descriptors);
        config_reader.AddEntryReader("FILE", file_entry_reader, true);

        config_reader.ReadConfig(args.file_cfg_name());

        for (auto file_descriptor : file_descriptors){ //loop on DYJets files
            const DYBinDescriptor file_descriptor_element = file_descriptor.second;
            std::cout << "File descriptor characteristics: " << file_descriptor.first << ", " <<
                         file_descriptor_element.file_path << ", " << file_descriptor_element.fileType
                      << std::endl;
            global_map.bin = file_descriptor_element;
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
                    global_map.gen_counts[genId] += nevents;
                    if (file_descriptor_element.fileType == analysis::FileType::inclusive)
                        inclusive.gen_counts[genId] += nevents;

                } // end loop on gen event info
            } //end loop on entries
            all_samples.push_back(global_map);
        } //end loop n file_descriptors
    }

    void LoadBins()
    {
       output_bins = DYBinDescriptor::LoadConfig(args.cfg_name());
    }

    void FindBestBin(DYBinDescriptor& output_bin) const
    {
        for(unsigned i = 0; i < all_samples.size(); ++i) {
            const SampleDescriptor sample = all_samples.at(i);
            const size_t contribution = sample.Integral(output_bin);
            if(!contribution) continue;
            //formula 2
            PhysicalValue nu ( contribution / double(sample.Integral()), sqrt(double(contribution)/double(sample.Integral())));

            if(!(sample.bin.fileType == analysis::FileType::inclusive)) {
                const size_t sample_contribution = inclusive.Integral(sample.bin);
                // formula 3
                PhysicalValue nu_incl(sample_contribution/double(inclusive.Integral()),
                                      sqrt(double(sample_contribution)/double(inclusive.Integral())));
                nu *= nu_incl;
            }

            if(output_bin.nu.GetStatisticalError() > nu.GetStatisticalError()) {
                output_bin.nu = nu;
                output_bin.ref_sample = sample.bin.name;
            }
        }

//        if(output_bin.nu.GetStatisticalError() == std::numeric_limits<double>::infinity())
//            throw exception("ref not found");
    }

    void CalcWeight(DYBinDescriptor& output_bin) const
    {
        output_bin.weight = output_bin.nu.GetValue() / global_map.Integral(output_bin);
    }

    public:


    void Run()
    {

        std::ofstream cfg(args.output_file());
        if(cfg.fail())
            throw analysis::exception("Unable to create outputfile'");
        cfg.exceptions(std::ofstream::failbit | std::ofstream::badbit);

        cfg << "#n_jet_min n_jet_max n_bjet_min n_bjet_max ht_bin_min ht_bin_max weight nu err_nu ref_sample\n";

        for(unsigned i = 0; i < output_bins.size(); ++i)
        {
            DYBinDescriptor output_bin = output_bins.at(i);
            FindBestBin(output_bin);
            CalcWeight(output_bin);
            cfg << output_bin.n_jet.min() << " " <<output_bin.n_jet.max()  << " " <<
                   output_bin.n_bjet.min() << " " << output_bin.n_bjet.max() << " " <<
                   output_bin.n_ht.min() << " " << output_bin.n_ht.max() << " " <<
                   output_bin.weight << " " << output_bin.nu << " " << output_bin.nu.GetStatisticalError() <<
                " " << output_bin.ref_sample <<  "\n";
        }

    }
private:
    Arguments args;
};

}

PROGRAM_MAIN(analysis::FileMerger, Arguments) // definition of the main program function
