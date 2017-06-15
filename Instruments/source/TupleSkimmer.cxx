/*! Skim EventTuple.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <thread>
#include <functional>

#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Run/include/EntryQueue.h"
#include "AnalysisTools/Core/include/ProgressReporter.h"
#include "h-tautau/McCorrections/include/EventWeights.h"
#include "hh-bbtautau/McCorrections/include/EventWeights_HH.h"

struct Arguments {
    REQ_ARG(std::string, treeNames);
    REQ_ARG(std::string, originalFileName);
    REQ_ARG(std::string, outputFileName);
	REQ_ARG(std::string, sample_type);
	OPT_ARG(std::string, additionalInputFile, "");
	OPT_ARG(std::string, additional2InputFile, "");
};

namespace analysis {

class TupleSkimmer {
public:
    using Event = ntuple::Event;
    using EventPtr = std::shared_ptr<Event>;
    using EventTuple = ntuple::EventTuple;
    using EventQueue = run::EntryQueue<EventPtr>;
	
	using ExpressTuple = ntuple::ExpressTuple;
	using ExpressEvent = ntuple::ExpressEvent;
	using ExpressPtr   = std::shared_ptr<ExpressEvent>;

    TupleSkimmer(const Arguments& _args) : args(_args), processQueue(100000), writeQueue(100000),
        eventWeights_HH(Parse<Channel>(args.treeName()),Period::Run2016, DiscriminatorWP::Medium)
    { }

    void Run()
    {
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,6,0)
        ROOT::EnableThreadSafety();
#endif

		outputFile = root_ext::CreateRootFile(args.outputFileName());
		originalFile = root_ext::OpenRootFile(args.originalFileName());
		inputFiles.push_back(originalFile);
		if (args.additionalInputFile() != "")
		{
			additionalFile = root_ext::OpenRootFile(args.additionalInputFile());
			inputFiles.push_back(additionalFile);
			std::cout << "  --> there is one additional file <-- " <<std::endl;
		}
		
		if (args.additional2InputFile() != "")
		{
			additional2File = root_ext::OpenRootFile(args.additional2InputFile());
			inputFiles.push_back(additional2File);
			std::cout << "  --> there is another additional file <-- " <<std::endl;
		}
        //BDT variables: "dphi_mumet", "dphi_metsv", "dR_taumu", "mT1", "mT2", "dphi_bbmet", "dphi_bbsv", "dR_bb", "m_bb",
        DisabledBranches_read = {"n_jets", "weight_btag", "weight_ttbar_pt",  "weight_PU", "weight_dy", "weight_ttbar_merge",
                                 "weight_sm", "shape_denominator_weight", "trigger_accepts", "trigger_matches"};

        DisabledBranches_write = { "lhe_particle_pdg", "lhe_particle_p4", "pfMET_cov", "genJets_partoFlavour",
                                   "genJets_hadronFlavour", "genJets_p4", "genParticles_p4", "genParticles_pdg",
                                   "trigger_accepts", "trigger_matches", "tauId_keys_1", "tauId_values_1",
                                   "tauId_keys_2", "tauId_values_2"};
		
		denominator = GetShapeDenominatorWeight(args.originalFileName());

        std::thread process_thread(std::bind(&TupleSkimmer::ProcessThread, this));
        std::thread writer_thread(std::bind(&TupleSkimmer::WriteThread, this, args.treeName(), args.outputFileName()));

        //ReadThread(args.treeName(), args.originalFileName());
		ReadThread2(args.treeName());

        std::cout << "Waiting for process and write threads to finish..." << std::endl;
        process_thread.join();
        writer_thread.join();
		
		SaveSummaryTree();
    }

private:
    std::pair<double,double> GetShapeDenominatorWeight(const std::string& /*originalFileName*/)
	{
		std::cout << "Calculating denominator for shape changing weights..." << std::endl;
        std::pair<double,double> tot_weight; //first without ttbar_pt, second with ttbar_pt
		std::shared_ptr<ExpressTuple> AllEventTuple(new ExpressTuple("all_events", originalFile.get(), true));
		
		const Long64_t all_entries = AllEventTuple->GetEntries();
		for(Long64_t current_entry = 0; current_entry < all_entries; ++current_entry)
		{
			AllEventTuple->GetEntry(current_entry);
			ExpressPtr event(new ExpressEvent(AllEventTuple->data()));
            double pu = eventWeights_HH.GetPileUpWeight(*event);
            double mc = event->genEventWeight;
			
            double weight_ttbar_pt = 1.;
            if(args.sample_type() == "ttbar")
            {
                weight_ttbar_pt = eventWeights_HH.GetTopPtWeight(*event);
            }

            double weight_sm = 1.;
            if(args.sample_type() == "sm")
            {
                weight_sm = eventWeights_HH.GetBSMtoSMweight(*event);
            }

            tot_weight.first = tot_weight.first + (pu*mc*weight_sm);
            tot_weight.second = tot_weight.second + (pu*weight_ttbar_pt*mc*weight_sm);
		}
		return tot_weight;
	}

    void ReadThread(const std::string& treeName, const std::string& /*originalFileName*/)
    {
		if(args.sample_type() == "data")
		{
			DisabledBranches_read.erase("trigger_accepts");
			DisabledBranches_read.erase("trigger_matches");
		}
		std::shared_ptr<EventTuple> originalTuple(new EventTuple(treeName, originalFile.get(), true, DisabledBranches_read));

		tools::ProgressReporter reporter(10, std::cout, "Starting skimming...");
		const Long64_t n_entries = originalTuple->GetEntries();
        reporter.SetTotalNumberOfEvents(static_cast<size_t>(n_entries));
		for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry)
		{
			originalTuple->GetEntry(current_entry);
            reporter.Report(static_cast<size_t>(current_entry));
			EventPtr event(new Event(originalTuple->data()));
			processQueue.Push(event);
		}
		processQueue.SetAllDone();
        reporter.Report(static_cast<size_t>(n_entries), true);
    }

	void ReadThread2(const std::string& treeName)
	{
		if(args.sample_type() == "data")
		{
			DisabledBranches_read.erase("trigger_accepts");
			DisabledBranches_read.erase("trigger_matches");
		}

		tools::ProgressReporter reporter(10, std::cout, "Starting skimming...");

		double tot_entries = 0.;
		std::vector<std::shared_ptr<EventTuple>> inputTuples;
		for (size_t n = 0; n<inputFiles.size(); ++n)
		{
			std::shared_ptr<EventTuple> tempTuple(new EventTuple(treeName, inputFiles.at(n).get(), true, DisabledBranches_read));
			inputTuples.push_back(tempTuple);
			std::cout << "   Temp entries: " << tempTuple->GetEntries() << std::endl;
			tot_entries += tempTuple->GetEntries();
		}

		std::cout << "   Total entries: " << tot_entries << std::endl;
        reporter.SetTotalNumberOfEvents(static_cast<size_t>(tot_entries));
		
		for (size_t n = 0; n<inputTuples.size(); ++n)
		{
			std::cout << " Analyzing Tuple: " << inputTuples.at(n) << std::endl;
			
			const Long64_t n_entries = inputTuples.at(n)->GetEntries();
			//reporter.SetTotalNumberOfEvents(n_entries);
			for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry)
			{
			 	inputTuples.at(n)->GetEntry(current_entry);
                reporter.Report(static_cast<size_t>(current_entry));
				EventPtr event(new Event(inputTuples.at(n)->data()));
				processQueue.Push(event);
			}
		}
		
		processQueue.SetAllDone();
        reporter.Report(static_cast<size_t>(tot_entries), true);
	}

	void SaveSummaryTree()
	{
		std::cout << "Copying the summary tree..." << std::endl;
		auto original_tree = root_ext::ReadObject<TTree>(*originalFile, "summary");

		outputFile->cd();
		TTree* copied_tree = original_tree->CloneTree();
		outputFile->Write();
		delete original_tree;
		delete copied_tree;
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

    void WriteThread(const std::string& treeName, const std::string& /*outputFileName*/)
    {
		if(args.sample_type() == "data")
		{
			DisabledBranches_write.erase("trigger_accepts");
			DisabledBranches_write.erase("trigger_matches");
		}
		std::shared_ptr<EventTuple> outputTuple(new EventTuple(treeName, outputFile.get(), false, DisabledBranches_write));

        EventPtr event;
        while(writeQueue.Pop(event)) {
            (*outputTuple)() = *event;
            outputTuple->Fill();
        }

        outputTuple->Write();
    }

    bool ProcessEvent(Event& event)
    {
	
        LorentzVectorE_Float first_jet = event.jets_p4.at(0);
        LorentzVectorE_Float second_jet = event.jets_p4.at(1);
        if (event.jets_p4.size() < 2 || event.extraelec_veto == false || event.extramuon_veto == false ||
                 std::abs(first_jet.eta()) >= cuts::btag_2016::eta ||
                 std::abs(second_jet.eta()) >= cuts::btag_2016::eta) return false;

        //const EventEnergyScale es = static_cast<EventEnergyScale>(event.eventEnergyScale);
        //if(es != EventEnergyScale::Central) return false;

        //old WP for isolation: "byTightIsolationMVArun2v1DBoldDMwLT", "byVTightIsolationMVArun2v1DBoldDMwLT"
        static const std::set<std::string> tauID_Names = {
            "againstMuonTight3", "againstElectronVLooseMVA6", "againstElectronTightMVA6", "againstMuonLoose3",
            "byMediumIsolationMVArun2v1DBoldDMwLT"
        };
		
//        decltype(event.tauIDs_1) tauIDs_1, tauIDs_2;
//        for(const auto& name : tauID_Names) {
//            if(event.tauIDs_1.count(name))
//                tauIDs_1[name] = event.tauIDs_1.at(name);
//            if(event.tauIDs_2.count(name))
//                tauIDs_2[name] = event.tauIDs_2.at(name);

//        }
//        event.tauIDs_1 = tauIDs_1;
//        event.tauIDs_2 = tauIDs_2;

	
		// Event Variables
        event.n_jets = static_cast<unsigned>(event.jets_p4.size());
        event.ht_other_jets = static_cast<float>(analysis::Calculate_HT(event.jets_p4.begin()+2,event.jets_p4.end()));

		// Event Weights Variables
        event.weight_btag = eventWeights_HH.GetBtagWeight(event);
        event.weight_lepton = eventWeights_HH.GetLeptonTotalWeight(event);
        event.weight_PU = eventWeights_HH.GetPileUpWeight(event);
        event.weight_ttbar_pt = eventWeights_HH.GetTopPtWeight(event);
        event.weight_dy = eventWeights_HH.GetDY_weight(event);
        event.weight_wjets = eventWeights_HH.GetWjets_weight(event);
        event.weight_ttbar_merge = eventWeights_HH.GetTTbar_weight(event);
        event.weight_sm = eventWeights_HH.GetBSMtoSMweight(event);
//		event.shape_denominator_weight = denominator;
	
		// Jets Variables Resizing
		event.jets_csv.resize(2);
		event.jets_rawf.resize(2);
		event.jets_mva.resize(2);
		event.jets_p4.resize(2);
		event.jets_partonFlavour.resize(2);
		event.jets_hadronFlavour.resize(2);
		
        return true;
    }

private:
    Arguments args;
    EventQueue processQueue, writeQueue;
	
	std::shared_ptr<TFile> originalFile;
	std::shared_ptr<TFile> outputFile;
	std::shared_ptr<TFile> additionalFile;
	std::shared_ptr<TFile> additional2File;
	std::vector<std::shared_ptr<TFile>> inputFiles;
	
    mc_corrections::EventWeights_HH eventWeights_HH;
    std::pair<double,double> denominator;
	
	std::set<std::string> DisabledBranches_read;
	std::set<std::string> DisabledBranches_write;
};

} // namespace analysis

PROGRAM_MAIN(analysis::TupleSkimmer, Arguments)
