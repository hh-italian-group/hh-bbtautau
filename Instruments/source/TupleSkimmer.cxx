/*! Skim EventTuple.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <thread>
#include <functional>

#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Run/include/EntryQueue.h"
#include "AnalysisTools/Core/include/ProgressReporter.h"

struct Arguments {
    REQ_ARG(std::string, treeName);
    REQ_ARG(std::string, originalFileName);
    REQ_ARG(std::string, outputFileName);
};

namespace analysis {

class TupleSkimmer {
public:
    using Event = ntuple::Event;
    using EventPtr = std::shared_ptr<Event>;
    using EventTuple = ntuple::EventTuple;
    using EventQueue = run::EntryQueue<EventPtr>;

    TupleSkimmer(const Arguments& _args) : args(_args), processQueue(100000), writeQueue(100000) {}

    void Run()
    {
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,6,0)
        ROOT::EnableThreadSafety();
#endif

        std::thread process_thread(std::bind(&TupleSkimmer::ProcessThread, this));
        std::thread writer_thread(std::bind(&TupleSkimmer::WriteThread, this, args.treeName(),
                                            args.outputFileName()));


        ReadThread(args.treeName(), args.originalFileName());

        std::cout << "Waiting for process and write threads to finish..." << std::endl;
        process_thread.join();
        writer_thread.join();
    }

private:
    void ReadThread(const std::string& treeName, const std::string& originalFileName)
    {
        auto originalFile = root_ext::OpenRootFile(originalFileName);
        std::shared_ptr<EventTuple> originalTuple(new EventTuple(treeName, originalFile.get(), true,
                { "lhe_n_partons", "lhe_HT" }));

        tools::ProgressReporter reporter(10, std::cout, "Starting skimming...");
        const Long64_t n_entries = originalTuple->GetEntries();
        reporter.SetTotalNumberOfEvents(n_entries);
        for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) {
            originalTuple->GetEntry(current_entry);
            reporter.Report(current_entry);
            EventPtr event(new Event(originalTuple->data()));
            processQueue.Push(event);
        }
        processQueue.SetAllDone();
        reporter.Report(n_entries, true);
    }

    void ProcessThread()
    {
        EventPtr event;
        while(processQueue.Pop(event)) {
            if(ProcessEvent(*event))
                writeQueue.Push(event);
        }
        writeQueue.SetAllDone();
    }

    void WriteThread(const std::string& treeName, const std::string& outputFileName)
    {
        auto outputFile = root_ext::CreateRootFile(outputFileName);
        std::shared_ptr<EventTuple> outputTuple(new EventTuple(treeName, outputFile.get(), false,
                { "lhe_particle_pdg", "lhe_particle_p4" } ));

        EventPtr event;
        while(writeQueue.Pop(event)) {
            (*outputTuple)() = *event;
            outputTuple->Fill();
        }

        outputTuple->Write();
    }

    static bool ProcessEvent(Event& event)
    {
        const EventEnergyScale es = static_cast<EventEnergyScale>(event.eventEnergyScale);
        if(es != EventEnergyScale::Central) return false;

        static const std::set<std::string> tauID_Names = {
            "againstMuonTight3", "againstElectronVLooseMVA6", "againstElectronTightMVA6", "againstMuonLoose3",
            "byTightIsolationMVArun2v1DBoldDMwLT", "byVTightIsolationMVArun2v1DBoldDMwLT"
        };
        decltype(event.tauIDs_1) tauIDs_1, tauIDs_2;
        for(const auto& name : tauID_Names) {
            if(event.tauIDs_1.count(name))
                tauIDs_1[name] = event.tauIDs_1.at(name);
            if(event.tauIDs_2.count(name))
                tauIDs_2[name] = event.tauIDs_2.at(name);

        }
        event.tauIDs_1 = tauIDs_1;
        event.tauIDs_2 = tauIDs_2;

        static const std::set<int> quarks_and_gluons = { 1, 2, 3, 4, 5, 6, 21 };
        if(event.lhe_particle_p4.size()) {
            double lhe_HT2 = 0;
            event.lhe_n_partons = 0;
            for(size_t n = 0; n < event.lhe_particle_p4.size(); ++n) {
                const int abs_pdg = std::abs(event.lhe_particle_pdg.at(n));
                if(!quarks_and_gluons.count(abs_pdg)) continue;
                lhe_HT2 += event.lhe_particle_p4.at(n).perp2();
                ++event.lhe_n_partons;
            }
            event.lhe_HT = std::sqrt(lhe_HT2);
            event.lhe_particle_p4.clear();
            event.lhe_particle_pdg.clear();
        }

        return true;
    }

private:
    Arguments args;
    EventQueue processQueue, writeQueue;
};

} // namespace analysis

PROGRAM_MAIN(analysis::TupleSkimmer, Arguments)
