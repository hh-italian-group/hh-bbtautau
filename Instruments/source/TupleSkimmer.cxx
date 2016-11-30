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
#include "h-tautau/McCorrections/include/EventWeights.h"

struct Arguments {
    REQ_ARG(std::string, treeName);
    REQ_ARG(std::string, originalFileName);
    REQ_ARG(std::string, outputFileName);
	REQ_ARG(std::string, sample_type);
};

namespace analysis {

class TupleSkimmer {
public:
    using Event = ntuple::Event;
    using EventPtr = std::shared_ptr<Event>;
    using EventTuple = ntuple::EventTuple;
    using EventQueue = run::EntryQueue<EventPtr>;

    TupleSkimmer(const Arguments& _args) : args(_args), processQueue(100000), writeQueue(100000), eventWeights(Period::Run2016, DiscriminatorWP::Medium) {}

    void Run()
    {
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,6,0)
        ROOT::EnableThreadSafety();
#endif

        std::thread process_thread(std::bind(&TupleSkimmer::ProcessThread, this, args.sample_type(), eventWeights));
        std::thread writer_thread(std::bind(&TupleSkimmer::WriteThread, this, args.treeName(), args.outputFileName()));

        ReadThread(args.treeName(), args.originalFileName());

        std::cout << "Waiting for process and write threads to finish..." << std::endl;
        process_thread.join();
        writer_thread.join();
		
		std::cout << "Copying the summary tree..." << std::endl;
		SaveSummaryTree(args.originalFileName(), args.outputFileName());
    }

private:
    void ReadThread(const std::string& treeName, const std::string& originalFileName)
    {
        auto originalFile = root_ext::OpenRootFile(originalFileName);
        std::shared_ptr<EventTuple> originalTuple(new EventTuple(treeName, originalFile.get(), true,
                { "lhe_n_partons", "lhe_HT", "dphi_mumet", "dphi_metsv", "dR_taumu", "mT1", "mT2", "dphi_bbmet", "dphi_bbsv", "dR_bb",
				"n_jets", "btag_weight", "ttbar_weight", "trigger_accepts", "trigger_matches"}));

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

	void SaveSummaryTree(const std::string& originalFileName, const std::string& outputFileName)
	{
		const char * in_name = originalFileName.c_str();
		TFile *input_file = new TFile(in_name);
		TTree *input_tree = (TTree*)input_file->Get("summary");
		
		const char * out_name = outputFileName.c_str();
		TFile *output_file = new TFile(out_name,"update");
		//TTree *output_tree =
		input_tree->CloneTree();
		output_file->Write();
		delete input_file;
		delete output_file;

	}


    void ProcessThread(const std::string& sample_type, analysis::mc_corrections::EventWeights eventWeights)
    {
        EventPtr event;
        while(processQueue.Pop(event)) {
            if(ProcessEvent(*event, sample_type, eventWeights))
                writeQueue.Push(event);
        }
        writeQueue.SetAllDone();
    }

    void WriteThread(const std::string& treeName, const std::string& outputFileName)
    {
        auto outputFile = root_ext::CreateRootFile(outputFileName);
        std::shared_ptr<EventTuple> outputTuple(new EventTuple(treeName, outputFile.get(), false,
                { "lhe_particle_pdg", "lhe_particle_p4", "pfMET_cov", "genJets_nTotal", "genJets_partoFlavour", "genJets_hadronFlavour",
				 "genJets_p4", "genParticles_p4", "genParticles_pdg"} ));

        EventPtr event;
        while(writeQueue.Pop(event)) {
            (*outputTuple)() = *event;
            outputTuple->Fill();
        }

        outputTuple->Write();
    }

    static bool ProcessEvent(Event& event, const std::string& sample_type, analysis::mc_corrections::EventWeights eventWeights)
    {
	
		if (event.jets_p4.size() < 2) return false;
	
        //const EventEnergyScale es = static_cast<EventEnergyScale>(event.eventEnergyScale);
        //if(es != EventEnergyScale::Central) return false;

        static const std::set<std::string> tauID_Names = {
            "againstMuonTight3", "againstElectronVLooseMVA6", "againstElectronTightMVA6", "againstMuonLoose3",
            "byTightIsolationMVArun2v1DBoldDMwLT", "byVTightIsolationMVArun2v1DBoldDMwLT"
        };
		
        /*decltype(event.tauIDs_1) tauIDs_1, tauIDs_2;
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
        }*/

		/*std::cout<<"--- Event ---"<<std::endl;
		for (UInt_t i=0; i<event.jets_p4.size(); i++)
		{
			std::cout<<"  jet "<<i<<" - pT: "<<event.jets_p4[i].Pt()<<" - csv: "<<(float)event.jets_csv[i]<<std::endl;
		}*/

		// Event Variables
		event.n_jets = event.jets_p4.size();

		// Event Weights Variables
		event.btag_weight = eventWeights.GetBtagWeight(event);
		
		
		double topWeight = 1.;
		if(sample_type == "ttbar") {
			for(size_t n = 0; n < event.genParticles_pdg.size(); ++n) {
				if(std::abs(event.genParticles_pdg.at(n)) != 6) continue;
				const double pt = event.genParticles_p4.at(n).pt();
				topWeight *= std::sqrt(std::exp(0.0615 - 0.0005 * pt));
			}
		}
		event.ttbar_weight = topWeight;
		
		// BDT Variables
		event.dphi_mumet = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1    , event.pfMET_p4));
		event.dphi_metsv = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4));
		event.dR_taumu = std::abs(ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2));
		event.mT1 = Calculate_MT(event.p4_1, event.pfMET_p4);
		event.mT2 = Calculate_MT(event.p4_2, event.pfMET_p4);
		
		if (event.jets_p4.size() >= 2)
		{
			ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>> bsum = event.jets_p4[0] + event.jets_p4[1];
			event.dphi_bbmet = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bsum, event.pfMET_p4));
			event.dphi_bbsv = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bsum, event.SVfit_p4));
			event.dR_bb = std::abs(ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]));
		}
		else
		{
			event.dphi_bbmet = -1.;
			event.dphi_bbsv = -1.;
			event.dR_bb = -1.;
		}
		
		// Jets Variables Resizing
		event.jets_csv.resize(2);
		event.jets_rawf.resize(2);
		event.jets_mva.resize(2);
		
        return true;
    }

private:
    Arguments args;
    EventQueue processQueue, writeQueue;
	mc_corrections::EventWeights eventWeights;
};

} // namespace analysis

PROGRAM_MAIN(analysis::TupleSkimmer, Arguments)
