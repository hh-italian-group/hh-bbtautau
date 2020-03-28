/*! Study of a kinematic fit performance.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "h-tautau/Analysis/include/Htautau_2015.h"

struct Arguments {
    REQ_ARG(std::string, input_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(std::string, event_id);
};

namespace analysis {

class KinFitStudy {
public:
    KinFitStudy(const Arguments& _args) : args(_args), eventId(args.event_id()), kinfitProducer(100) {}

    void Run()
    {
        auto inputFile = root_ext::OpenRootFile(args.input_file());
        ntuple::EventTuple eventTuple(args.tree_name(), inputFile.get(), true, { "lhe_n_partons", "lhe_HT" });
        const Long64_t n_entries = eventTuple.GetEntries();
        for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) {
            eventTuple.GetEntry(current_entry);
            const auto bjet_pair = EventInfo::SelectBjetPair(eventTuple.data(), cuts::Htautau_2015::btag::pt,
                                                                 cuts::Htautau_2015::btag::eta, JetOrdering::CSV);
            EventInfo event(eventTuple.data(), bjet_pair);
            if(eventId != event.GetEventId()) continue;
            if(event.GetEnergyScale() != EventEnergyScale::Central) continue;

            kinfitProducer.Fit(event->p4_1, event->p4_2, event.GetHiggsBB().GetFirstDaughter().GetMomentum(),
                               event.GetHiggsBB().GetSecondDaughter().GetMomentum(), event.GetMET());
            return;
        }
    }

private:
    Arguments args;
    EventIdentifier eventId;
    kin_fit::FitProducer kinfitProducer;
};

} // namespace analysis

PROGRAM_MAIN(analysis::KinFitStudy, Arguments)
