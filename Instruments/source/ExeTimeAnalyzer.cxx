/*! Merge DYJets files.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include <iostream>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"


struct Arguments {
    run::Argument<std::string> input{"input", "Input file."};
};

namespace analysis {

class ExeTimeAnalyzer {
public:
    ExeTimeAnalyzer(const Arguments& _args) : args(_args)
    {
    }

    void Run()
    {
        auto inputFile = root_ext::OpenRootFile(args.input());
        auto summaryTuple = ntuple::CreateSummaryTuple("summary", inputFile.get(), true, ntuple::TreeState::Full);

        auto n_process_event = std::make_shared<TH1F>("n_process_event","number processed events",50000,-0.5,50000-0.5);
        auto exeTime = std::make_shared<TH1F>("exeTime","exeTime",10000,-0.5,100000-0.5);
        auto time_per_event = std::make_shared<TH1F>("time_per_event","time_per_event",10000,-0.5,100000-0.5);

        for(const auto& summary : *summaryTuple) {
            n_process_event->Fill(summary.numberOfProcessedEvents);
            exeTime->Fill(summary.exeTime);
            if(summary.numberOfProcessedEvents == 0) continue;
            double time_event = double(summary.exeTime)/summary.numberOfProcessedEvents;
            time_per_event->Fill(time_event);
        }

        double common_exe_time = 1.0;
        double full_integral = exeTime->Integral(0,exeTime->GetXaxis()->GetNbins()+1);
        for(Int_t n = 0; n < exeTime->GetXaxis()->GetNbins(); ++n){
            if (exeTime->Integral(0,n)/full_integral > 0.995){
                common_exe_time = exeTime->GetBinCenter(n);
                break;
            }
        }

        static constexpr int crab_dead_time = 12 * 60 * 60; // s
        static constexpr char sep = ','; // ,

        double scale_factor = crab_dead_time/common_exe_time;
        size_t n_evt_per_job_prod_v2 = static_cast<size_t>(std::ceil(n_process_event->GetMean()));
        size_t n_evt_per_job_prod_v3 = static_cast<size_t>(std::ceil(scale_factor * n_evt_per_job_prod_v2));
        double time_event_calc = time_per_event->GetMean();

        std::cout << args.input() << sep << scale_factor << sep << n_evt_per_job_prod_v2 << sep
                  << n_evt_per_job_prod_v3 << sep
                  << time_event_calc << std::endl;
    }

private:
    Arguments args;
};

} //namespace analysis

PROGRAM_MAIN(analysis::ExeTimeAnalyzer, Arguments)
