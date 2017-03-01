/*! My first study
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_signal_file); // required argument "input_signal_file"
    REQ_ARG(std::string, input_bkg_file); // required argument "input_bkg_file"
    REQ_ARG(std::string, output_file); // required argument "output_file"
};

namespace analysis {


class MvaData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    TH1D_ENTRY(PtTau, 100, 0, 500)
    TH1D_ENTRY(PhiTau, 100, -3.5, 3.5)
    TH1D_ENTRY(EtaTau, 100, -2.5, 2.5)
    TH2D_ENTRY(MttMbb, 100, 0, 500, 100, 0, 500)
    TH2D_ENTRY(MttMbb_cut, 100, 0, 500, 100, 0, 500)
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
        args(_args), input_signal(root_ext::OpenRootFile(args.input_signal_file())), input_bkg(root_ext::OpenRootFile(args.input_bkg_file())),
        output(root_ext::CreateRootFile(args.output_file())),
        anaData(output), tuple_signal("muTau", input_signal.get(), true, GetDisabledBranches()), tuple_bkg("muTau", input_bkg.get(), true, GetDisabledBranches())
    {
    }

    void Run()
    {
        std::cout << boost::format("Processing input signal file '%1%' and input signal file '%2%' into output file '%3%'.\n")
                     % args.input_signal_file()   % args.input_bkg_file() % args.output_file();
        const Long64_t n_entries_signal = tuple_signal.GetEntries(), n_entries_bkg = tuple_bkg.GetEntries();
        for(Long64_t current_entry = 0; current_entry < n_entries_signal; ++current_entry) {
            tuple_signal.GetEntry(current_entry);
            const Event& event = tuple_signal.data();

            if (event.eventEnergyScale!=0 || (event.q_1+event.q_2)!=0 || event.jets_p4.size() < 2
                || event.extraelec_veto==true || event.extramuon_veto==true) continue;


            anaData.PtTau("signal").SetMarkerColor(4);
            anaData.PtTau("signal").Fill(event.p4_1.pt());

            anaData.PhiTau("signal").SetMarkerColor(4);
            anaData.PhiTau("signal").Fill(event.p4_1.phi());

            anaData.EtaTau("signal").SetMarkerColor(4);
            anaData.EtaTau("signal").Fill(event.p4_1.eta());

            LorentzVectorE_Float bb= event.jets_p4[0] + event.jets_p4[1];
            anaData.MttMbb("signal").Fill(event.SVfit_p4.mass(),bb.M());

            double circular_cut=std::sqrt(pow(event.SVfit_p4.mass()-116.,2)+pow(bb.M()-111,2));
            if (circular_cut>40) continue;
            anaData.MttMbb_cut("signal").Fill(event.SVfit_p4.mass(),bb.M());

        }

        Double_t scale = 1/anaData.PtTau("signal").Integral();
        anaData.PtTau("signal").Scale(scale);
        scale = 1/anaData.PhiTau("signal").Integral();
        anaData.PhiTau("signal").Scale(scale);
        scale = 1/anaData.EtaTau("signal").Integral();
        anaData.EtaTau("signal").Scale(scale);

        for(Long64_t current_entry = 0; current_entry < n_entries_bkg; ++current_entry) {
            tuple_bkg.GetEntry(current_entry);
            const Event& event = tuple_bkg.data();
            if (event.eventEnergyScale!=0 || (event.q_1+event.q_2)!=0 || event.jets_p4.size() < 2
                || event.extraelec_veto==true || event.extramuon_veto==true) continue;

            anaData.PtTau("bkg").SetMarkerColor(2);
            anaData.PtTau("bkg").Fill(event.p4_1.pt());

            anaData.PhiTau("bkg").SetMarkerColor(2);
            anaData.PhiTau("bkg").Fill(event.p4_1.phi());

            anaData.EtaTau("bkg").SetMarkerColor(2);
            anaData.EtaTau("bkg").Fill(event.p4_1.eta());

            LorentzVectorE_Float bb= event.jets_p4[0] + event.jets_p4[1];
            anaData.MttMbb("bkg").Fill(event.SVfit_p4.mass(),bb.M());
            double circular_cut=std::sqrt(pow(event.SVfit_p4.mass()-116.,2)+pow(bb.M()-111,2));
            if (circular_cut>40) continue;
            anaData.MttMbb_cut("bkg").Fill(event.SVfit_p4.mass(),bb.M());


        }
        scale = 1/anaData.PtTau("bkg").Integral();
        anaData.PtTau("bkg").Scale(scale);
        scale = 1/anaData.PhiTau("bkg").Integral();
        anaData.PhiTau("bkg").Scale(scale);
        scale = 1/anaData.EtaTau("bkg").Integral();
        anaData.EtaTau("bkg").Scale(scale);


    }
private:
    Arguments args;
    std::shared_ptr<TFile> input_signal,input_bkg, output;
    MvaData anaData;
    EventTuple tuple_signal, tuple_bkg;
};


}


PROGRAM_MAIN(analysis::MvaPreparation, Arguments) // definition of the main program function
// ./run.sh MvaPreparation --input_signal_file ~/Desktop/tuples/GluGluToRadionToHHTo2B2Tau_M-250_narrow.root
//   --input_bkg_file ~/Desktop/tuples/TT_ext3_muTau.root  --output_file mr.root
