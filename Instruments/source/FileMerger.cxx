/*! Merge DYJets files.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h" // definition of wrappers for the program main and program arguments.
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "Instruments/include/FileConfigEntryReader.h"

struct Arguments { // list of all program arguments
    //REQ_ARG(std::string, input_file); // required argument "input_file"
    REQ_ARG(std::string, tree_name); // tree name
    REQ_ARG(std::string, cfg_name); // cfg name
    REQ_ARG(std::string, file_cfg_name); // cfg name
};

namespace analysis {

struct SampleDescriptor {
    using GenMap = std::map<ntuple::GenId, size_t>;
    DYBinDescriptor bin;
    GenMap gen_counts;

    size_t Integral() const;
    size_t Integral(const DYBinDescriptor& range) const
    {

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
    }


private:
    SampleDescriptor global_map;
    SampleDescriptor inclusive;
    VSD all_samples;
    VBD output_bins;

    void LoadInputs()
    {
        // Fill inclusive and all_samples
    }

    void LoadBins()
    {
       output_bins = DYBinDescriptor::LoadConfig(args.cfg_name);
    }

    void FindBestBin(DYBinDescriptor& output_bin) const
    {
        for(sample : all_samples) {
            const size_t contribution = sample.Integrate(output_bin);
            if(!contribution) continue;
            PhysicalValue nu ( contribution / double(sample.Inegral()), sqrt(double(con)/integral));
            if(sample == inclusive) {
                // formula 2
                PhysicalValue nu = nu_excl;
            } else {
                const size_t sample_contirbution = inclusive.Integral(sample.bin);
                // formula 3
                PV nu_inc(sample_contirbution/incl.Integral);
                nu *= nu_inc;
            }

            if(oputput_bin.nu.Errror() > nu.error) {
                output_bin.nu = nu;
                output_bin.ref_sample = sample.name;
            }
        }

        if(output_bin.nu.Error == inf)
            throw exception("ref not found");
    }

    void CalcWeight(out_bin) const
    {
        out_bin.w = out_bin.nu / global_map.Integral(out_bin);
    }

    public:


    void Run()
    {

        for(out_bin : output_bins)
        {
            FindBestBin(out_bin);
            CalcWeight(out_bin);
            DBV out_bin_descs;
            for(sample : all_samples) {

            }

        }

        DYDescriptor::SaveCfg(args.ouput);


        // analyzer code
        analysis::ConfigReader config_reader;

        analysis::FileDescriptorCollection file_descriptors;
        analysis::FileConfigEntryReader file_entry_reader(file_descriptors);
        config_reader.AddEntryReader("FILE", file_entry_reader, true);

        config_reader.ReadConfig(args.file_cfg_name());

        JetSplitting::JetParametersMap jetParametersMap = jetSplitting.LoadConfig(args.cfg_name());


        using JetParameters_nuvalues_Map = std::map<JetSplitting::JetParameters, exclusive_samples_info::exclusive_samples_vector>;

        JetParameters_nuvalues_Map jetParameters_nuvalues_Map;

        for (auto file_descriptor : file_descriptors){ //loop on DYJets files
            const FileDescriptor file_descriptor_element = file_descriptor.second;
            std::cout << "File descriptor characteristics: " << file_descriptor.first << ", " <<
                         file_descriptor_element.file_path << ", " << file_descriptor_element.fileType
                      << std::endl;
            //if(!(file_descriptor_element.fileType == analysis::FileType::inclusive)) continue;

            const std::string outputFileName = "splitting" + file_descriptor_element.name + ".txt";
            std::ofstream cfg(outputFileName);
            if(cfg.fail())
                throw analysis::exception("Unable to create outputfile'");
            cfg.exceptions(std::ofstream::failbit | std::ofstream::badbit);

            cfg << "#n_jet_min n_jet_max n_bjet_min n_bjet_max ht_bin_min ht_bin_max n_events nu eps_nu\n";

            //auto inputFile = root_ext::OpenRootFile(args.input_file());
            auto inputFile = root_ext::OpenRootFile(file_descriptor_element.file_path);
            ntuple::SummaryTuple summaryTuple(args.tree_name(), inputFile.get(), true);
            const Long64_t n_entries = summaryTuple.GetEntries();
//            std::cout << "N entries in Summary Tuple " << n_entries << std::endl;

            ntuple::GenEventCountMap genEventCountMap;
            size_t total_n_processed_events = 0;
            for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on entries
                summaryTuple.GetEntry(current_entry);
//                std::cout << "N events processed in Summary Tuple " << summaryTuple.data().numberOfProcessedEvents << std::endl;
                //std::shared_ptr<analysis::SummaryInfo> summaryInfo(new analysis::SummaryInfo(summaryTuple.data()));

// da rimuovere tale struttura
                LHE_event_info::LHE_event_infos lhe_event_infos;
                total_n_processed_events += summaryTuple.data().numberOfProcessedEvents;
                for (size_t i = 0; i < summaryTuple.data().lhe_n_partons.size(); ++i){ // loop on gen event info
                    const ntuple::GenId genId(summaryTuple.data().lhe_n_partons.at(i),
                                              summaryTuple.data().lhe_n_b_partons.at(i),
                                              summaryTuple.data().lhe_ht10_bin.at(i));
                    size_t nevents = summaryTuple.data().lhe_n_events.at(i);
                    genEventCountMap[genId] += nevents;

                    LHE_event_info lhe_event_info(summaryTuple.data().lhe_n_partons.at(i),
                                 summaryTuple.data().lhe_n_b_partons.at(i),
                                 summaryTuple.data().lhe_ht10_bin.at(i),
                                 summaryTuple.data().lhe_n_events.at(i));
//                    std::cout << lhe_event_info << std::endl;
                    lhe_event_infos.push_back(lhe_event_info);
                } // end loop on gen event info

            } //end loop on entries
            std::cout << "Total processed events in the sumple: " << total_n_processed_events << std::endl;
            std::cout << "Print out map" << std::endl;
            size_t total_n_events = 0;
            for (auto genEventCount : genEventCountMap){
                const ntuple::GenId genId = genEventCount.first;
                const size_t nevents = genEventCount.second;
                total_n_events += nevents;
                std::cout << "GenId: " << genId.n_partons << " " << genId.n_b_partons << " " << genId.ht10_bin <<
                             ", nevents: " << nevents << std::endl;
            }
            std::cout << "Total events in the sumple: " << total_n_events << std::endl;


            for (auto jetParam : jetParametersMap){ // loop on splitting
                JetSplitting::JetParameters jetParameters = jetParam.second;
//                std::cout << "param ht min " << jetParameters.n_ht.min() << ",  max " << jetParameters.n_ht.max() << std::endl;
                size_t totalEvents = 0;
                for (auto genEventCount : genEventCountMap){
                    const ntuple::GenId genId = genEventCount.first;
                    const size_t nevents = genEventCount.second;
                    if(jetParameters.n_jet.Contains(genId.n_partons) &&
                          jetParameters.n_bjet.Contains(genId.n_b_partons) &&
                          jetParameters.n_ht.Contains(genId.ht10_bin))
                        totalEvents += nevents;
                    continue;

                }
                const PhysicalValue nu(((double)totalEvents)/total_n_events, 1/(std::sqrt(totalEvents)));
                double eps_nu = 1/(std::sqrt(totalEvents));

                if(file_descriptor_element.fileType == analysis::FileType::inclusive){
                    const PhysicalValue nu_incl(((double)totalEvents)/total_n_events,std::sqrt((double)totalEvents)/total_n_events);
                }

                if(file_descriptor_element.fileType == analysis::FileType::exclusive){
                    const PhysicalValue nu_excl(((double)totalEvents)/total_n_events,std::sqrt((double)totalEvents)/total_n_events);
                }

                const PhysicalValue nu_prod = nu_incl * nu_excl;
                const double rel_nu_err = nu_prod.GetRelativeStatisticalError();

                exclusive_samples_info excl_info(nu_prod,rel_nu_err,totalEvents,file_descriptor_element.name);
                exclusive_samples_info excl_info(nu,eps_nu,totalEvents,file_descriptor_element.name);
                jetParameters_nuvalues_Map[jetParameters].push_back(excl_info);

                cfg << jetParameters.n_jet.min() << " " <<jetParameters.n_jet.max()  << " " <<
                       jetParameters.n_bjet.min() << " " << jetParameters.n_bjet.max() << " " <<
                       jetParameters.n_ht.min() << " " << jetParameters.n_ht.max() << " " <<
                       totalEvents << " " << nu << " " << eps_nu <<  "\n";
            }
        }

        for (auto splitting_param : jetParameters_nuvalues_Map){
            JetSplitting::JetParameters jetParameters = splitting_param.first;
            exclusive_samples_info::exclusive_samples_vector excl_vector = splitting_param.second;
            std::cout << jetParameters.n_jet.min() << " " <<jetParameters.n_jet.max()  << " " <<
                         jetParameters.n_bjet.min() << " " << jetParameters.n_bjet.max() << " " <<
                         jetParameters.n_ht.min() << " " << jetParameters.n_ht.max() << ", size excl vector: " <<
                         excl_vector.size() << std::endl;
        }

    }
private:
    Arguments args;
    JetSplitting jetSplitting;
};

}

PROGRAM_MAIN(analysis::FileMerger, Arguments) // definition of the main program function
