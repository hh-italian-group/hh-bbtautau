/*! Merge DYJets files.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h" // definition of wrappers for the program main and program arguments.
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "Instruments/include/FileConfigEntryReader.h"

struct Arguments { // list of all program arguments
    //REQ_ARG(std::string, input_file); // required argument "input_file"
    REQ_ARG(std::string, tree_name); // tree name
    REQ_ARG(std::string, cfg_name); // cfg name
    REQ_ARG(std::string, file_cfg_name); // cfg name
};

struct JetSplitting{
public:
    struct JetParameters{
      int n_jet_min;
      double n_jet_max;
      double n_bjet_min;
      double n_bjet_max;
      double n_ht_min;
      double n_ht_max;
    };

    using JetParametersMap = std::map<size_t, JetParameters>;

    static JetParametersMap LoadConfig(const std::string& config_file)
    {
        static std::mutex m;
        std::lock_guard<std::mutex> lock(m);
        JetParametersMap jetParameters;
        jetParameters.clear();
        size_t line_number = 0;
        std::ifstream cfg(config_file);
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            if (!cfgLine.size() || cfgLine.at(0) == '#') continue;
            ++line_number;
            std::istringstream ss(cfgLine);
            JetParameters jetParam;
            ss >> jetParam.n_jet_min >> jetParam.n_jet_max >> jetParam.n_bjet_min >> jetParam.n_bjet_max >>
                    jetParam.n_ht_min >> jetParam.n_ht_max;
            jetParameters[line_number] = jetParam;
        }
        return jetParameters;
    }

};

struct LHE_event_info{
public:

    LHE_event_info(const double lhe_n_partons, const double lhe_n_b_partons, const double lhe_ht_bin,
                   const double lhe_n_events) : _lhe_n_partons(lhe_n_partons), _lhe_n_b_partons(lhe_n_b_partons),
        _lhe_ht_bin(lhe_ht_bin), _lhe_n_events(lhe_n_events) {}

    using LHE_event_infos = std::vector<LHE_event_info>;

    const double lhe_n_partons() const {return _lhe_n_partons; }
    const double lhe_n_b_partons() const {return _lhe_n_b_partons; }
    const double lhe_ht_bin() const {return _lhe_ht_bin; }
    const double lhe_n_events() const {return _lhe_n_events; }


private:
    double _lhe_n_partons; //stesso_tipo_delle ntuple
    double _lhe_n_b_partons;
    double _lhe_ht_bin;
    double _lhe_n_events;

};

std::ostream& operator<<(std::ostream& s, const LHE_event_info& lhe_event_info)
{
    s << lhe_event_info.lhe_n_partons() << " " << lhe_event_info.lhe_n_b_partons() << " " <<
         lhe_event_info.lhe_ht_bin() << " " << lhe_event_info.lhe_n_events();
    return s;
}

namespace analysis {

class FileMerger { // simple analyzer definition
public:
    FileMerger(const Arguments& _args) : args(_args)
    {
        // Analyzer initialization (e.g. open input/output files, parse configs...)
    }
    void Run()
    {
        // analyzer code
        analysis::ConfigReader config_reader;

        analysis::FileDescriptorCollection file_descriptors;
        analysis::FileConfigEntryReader file_entry_reader(file_descriptors);
        config_reader.AddEntryReader("FILE", file_entry_reader, true);

        config_reader.ReadConfig(args.file_cfg_name());

        JetSplitting::JetParametersMap jetParametersMap = jetSplitting.LoadConfig(args.cfg_name());
        const std::string outputFileName = "splitting.txt";
        std::ofstream cfg(outputFileName);
        if(cfg.fail())
            throw analysis::exception("Unable to create outputfile'");
        cfg.exceptions(std::ofstream::failbit | std::ofstream::badbit);

        cfg << "#n_jet_min n_jet_max	n_bjet_min n_bjet_max ht_bin_min ht_bin_max inclusive_n_events\n";

        for (auto file_descriptor : file_descriptors){ //loop on DYJets files
            const FileDescriptor file_descriptor_element = file_descriptor.second;
            std::cout << "File descriptor characteristics: " << file_descriptor.first << ", " <<
                         file_descriptor_element.file_path << ", " << file_descriptor_element.fileType
                      << std::endl;
            if(!(file_descriptor_element.fileType == analysis::FileType::inclusive)) continue;
            //auto inputFile = root_ext::OpenRootFile(args.input_file());
            auto inputFile = root_ext::OpenRootFile(file_descriptor_element.file_path);
            ntuple::SummaryTuple summaryTuple(args.tree_name(), inputFile.get(), true);
            const Long64_t n_entries = summaryTuple.GetEntries();
            std::cout << "N entries in Summary Tuple " << n_entries << std::endl;

            for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on events
                summaryTuple.GetEntry(current_entry);
//                std::cout << "N events processed in Summary Tuple " << summaryTuple.data().numberOfProcessedEvents << std::endl;
                //std::shared_ptr<analysis::SummaryInfo> summaryInfo(new analysis::SummaryInfo(summaryTuple.data()));
    //            std::cout << "N partons in Summary Tuple " << summaryTuple.data().lhe_n_partons.size() //or summaryTuple.operator()()
    //                      << std::endl;
    //            if (summaryTuple.data().lhe_n_partons.size() == summaryTuple.data().lhe_n_partons.size())
    //                std::cout << "Same size!" << std::endl;        

                LHE_event_info::LHE_event_infos lhe_event_infos;
                    for (size_t i = 0; i < summaryTuple.data().lhe_n_partons.size(); ++i){ // loop on gen event info
                        LHE_event_info lhe_event_info(summaryTuple.data().lhe_n_partons.at(i),
                                     summaryTuple.data().lhe_n_b_partons.at(i),
                                     summaryTuple.data().lhe_ht10_bin.at(i),
                                     summaryTuple.data().lhe_n_events.at(i));
//                        std::cout << lhe_event_info << std::endl;
                        lhe_event_infos.push_back(lhe_event_info);
                    }
//                std::cout << "LHE infos size: " << lhe_event_infos.size() << std::endl;

                for (auto jetParam : jetParametersMap){ // loop on splitting
                    JetSplitting::JetParameters jetParameters = jetParam.second;

                    for (size_t m = 0; m < lhe_event_infos.size(); ++m){
                        LHE_event_info lhe_event_element = lhe_event_infos.at(m);
                        if (lhe_event_element.lhe_n_partons() >= jetParameters.n_jet_min
                                && lhe_event_element.lhe_n_partons() <= jetParameters.n_jet_max
                                && lhe_event_element.lhe_n_b_partons() >= jetParameters.n_bjet_min &&
                                lhe_event_element.lhe_n_b_partons() <= jetParameters.n_bjet_max
                                && lhe_event_element.lhe_ht_bin() >= jetParameters.n_ht_min &&
                                lhe_event_element.lhe_ht_bin() <= jetParameters.n_ht_max)
//                            std::cout << "FOUND!!" << lhe_event_element.lhe_ht_bin()  << " " <<
//                                      jetParameters.n_ht_min << " " << jetParameters.n_ht_max << ", b partons: " <<
//                                         lhe_event_element.lhe_n_b_partons() << " " <<  jetParameters.n_bjet_min << " " <<
//                                         jetParameters.n_bjet_max << ", partons: " <<
//                                         lhe_event_element.lhe_n_partons() << " " <<  jetParameters.n_jet_min << " " <<
//                                         jetParameters.n_jet_max << std::endl;
                            cfg << jetParameters.n_jet_min << " " <<jetParameters.n_jet_max  << " " <<
                                   jetParameters.n_bjet_min << " " << jetParameters.n_bjet_max << " " <<
                                   jetParameters.n_ht_min << " " << jetParameters.n_ht_max << " " <<
                                   lhe_event_element.lhe_n_events() << "\n";
                    }
                }
            }
        }

    }
private:
    Arguments args;
    JetSplitting jetSplitting;
};

}

PROGRAM_MAIN(analysis::FileMerger, Arguments) // definition of the main program function
