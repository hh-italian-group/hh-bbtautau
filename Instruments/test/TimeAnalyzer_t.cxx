/*! Merge DYJets files.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include <iostream>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "Instruments/include/TimeFileDescriptor.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"



struct Arguments {
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> input_path{"input_path", "Input path of the samples"};
    run::Argument<std::string> cfg_name{"cfg_name", "file path cfg"};
    run::Argument<std::string> output_file{"output_file", "Output file"};
};

namespace analysis {

class TimeAnalyzer_t {
public:
    using VectorTimeFileDescriptor = std::vector<TimeFileDescriptor>;

    TimeAnalyzer_t(const Arguments& _args) : args(_args)
    {
        std::cout << "Starting..." << std::endl;
        outputs = TimeFileDescriptor::LoadConfig(args.cfg_name());
    }

    void Run()
    {
        for (TimeFileDescriptor& output : outputs){

            std::cout << "File name: " << output.file << std::endl;
            std::string filename = args.input_path()  + "/" + output.file;
            auto inputFile = root_ext::OpenRootFile(filename);
            ntuple::SummaryTuple summaryTuple(args.tree_name(), inputFile.get(), true);

            const Long64_t n_entries = summaryTuple.GetEntries();
            TH1F *n_process_event   = new TH1F("n_process_event","number processed events",1000,0,50000);
            TH1F *exeTime   = new TH1F("exeTime","exeTime",100,0,5000);

            for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on entries
                summaryTuple.GetEntry(current_entry);
                n_process_event->Fill(summaryTuple.data().numberOfProcessedEvents);
                exeTime->Fill(summaryTuple.data().exeTime);

            } //end loop on entries

            int crab_dead_time = 39450;
            double integral_calc = 0.999 * exeTime->Integral();
            output.scale_factor = crab_dead_time/integral_calc;
            output.n_evt_per_job_prod_v2 = static_cast<size_t>(n_process_event->GetMean()/exeTime->GetEntries());
            output.n_evt_per_job_prod_v3 = static_cast<size_t>(output.scale_factor * output.n_evt_per_job_prod_v2);
            std::cout << "scale factor: " << output.scale_factor << std::endl;
            std::cout << "n_evt_per_job_prod_v2: " << output.n_evt_per_job_prod_v2 << std::endl;
            std::cout << "n_evt_per_job_prod_v3: " << output.n_evt_per_job_prod_v3 << std::endl;

        }
        TimeFileDescriptor::SaveCfg(args.output_file(), outputs);
        std::cout << "Done SaveCfg" << std::endl;
    }


private:
    Arguments args;
    VectorTimeFileDescriptor outputs;
};

} //namespace analysis

PROGRAM_MAIN(analysis::TimeAnalyzer_t, Arguments)
