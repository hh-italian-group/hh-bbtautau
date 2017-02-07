/*! My first study
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_file); // required argument "input_file"
    REQ_ARG(std::string, output_file); // required argument "output_file"
};

namespace analysis {

class MvaData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;

    TH1D_ENTRY(PtTau, 100, 0, 200)
    TH1D_ENTRY(PhiTau, 100, -3.5, 3.5)
    TH1D_ENTRY(EtaTau, 100, -2.5, 2.5)
};


class MvaPreparation { // simple analyzer definition
public:
    using Event = ntuple::Event;
    using EventPtr = std::shared_ptr<Event>;
    using EventTuple = ntuple::EventTuple;

    static const std::set<std::string>& GetDisabledBranches()
    {
        static const std::set<std::string> DisabledBranches_read = {
            "dphi_mumet", "dphi_metsv", "dR_taumu", "mT1", "mT2", "dphi_bbmet", "dphi_bbsv", "dR_bb", "m_bb", "n_jets",
            "btag_weight", "ttbar_weight",  "PU_weight", "shape_denominator_weight", "trigger_accepts", "trigger_matches"
        };
        return DisabledBranches_read;
    }

    MvaPreparation(const Arguments& _args) :
        args(_args), input(root_ext::OpenRootFile(args.input_file())),
        output(root_ext::CreateRootFile(args.output_file())),
        anaData(output), tuple("muTau", input.get(), true, GetDisabledBranches())
    {
    }

    void Run()
    {
        std::cout << boost::format("Processing input file '%1%' into output file '%2%'.\n")
                     % args.input_file() % args.output_file();
        const Long64_t n_entries = tuple.GetEntries();
        for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) {
            tuple.GetEntry(current_entry);
            const Event& event = tuple.data();
            anaData.PtTau("signal").Fill(event.p4_1.pt());
            //anaData.PtTau("bkg").Fill(event.p4_2.pt());
            anaData.PhiTau("signal").Fill(event.p4_1.phi());
            anaData.EtaTau("signal").Fill(event.p4_1.eta());
        }
    }
private:
    Arguments args;
    std::shared_ptr<TFile> input, output;
    MvaData anaData;
    EventTuple tuple;
};

}


PROGRAM_MAIN(analysis::MvaPreparation, Arguments) // definition of the main program function
