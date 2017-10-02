/*! Skim EventTuple.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <thread>
#include <functional>
#include <random>

#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Run/include/EntryQueue.h"
#include "AnalysisTools/Core/include/ProgressReporter.h"
#include "h-tautau/McCorrections/include/EventWeights.h"
#include "hh-bbtautau/McCorrections/include/EventWeights_HH.h"
#include "hh-bbtautau/Instruments/include/SkimmerConfigEntryReader.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "h-tautau/Analysis/include/EventLoader.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"

struct Arguments {
    REQ_ARG(std::string, cfg);
    REQ_ARG(std::string, inputPath);
    REQ_ARG(std::string, outputPath);
    REQ_ARG(std::string, jobs);
    REQ_ARG(std::string, setup_name);
    OPT_ARG(unsigned, n_threads, 1);
};

namespace analysis {
namespace tuple_skimmer {
class TupleSkimmer {
public:
    using Event = ntuple::Event;
    using EventPtr = std::shared_ptr<Event>;
    using EventTuple = ntuple::EventTuple;
    using EventQueue = run::EntryQueue<EventPtr>;
	
	using ExpressTuple = ntuple::ExpressTuple;
	using ExpressEvent = ntuple::ExpressEvent;
    using SummaryTuple = ntuple::SummaryTuple;
    using ProdSummary = ntuple::ProdSummary;

    static constexpr size_t max_queue_size = 100000;

    TupleSkimmer(const Arguments& _args) :
        args(_args), processQueue(max_queue_size), writeQueue(max_queue_size)
    {
        std::cout << "TupleSkimmer started.\nReading config... " << std::flush;
        ROOT::EnableThreadSafety();
        if(args.n_threads() > 1)
            ROOT::EnableImplicitMT(args.n_threads());

        ConfigReader configReader;
        SetupCollection setups;
        SetupEntryReader setupReader(setups);
        configReader.AddEntryReader("SETUP", setupReader, false);
        SkimJobCollection all_jobs;
        SkimJobEntryReader jobReader(all_jobs);
        configReader.AddEntryReader("JOB", jobReader, true);
        configReader.ReadConfig(args.cfg());

        if(!setups.count(args.setup_name()))
            throw exception("Tuple skimmer setup not found.");
        setup = setups.at(args.setup_name());
        setup.UpdateTauIdHashes();

        std::cout << "done.\nLoading weights... " << std::flush;
        eventWeights_HH = std::make_shared<mc_corrections::EventWeights_HH>(setup.period, setup.btag_wp);
        std::cout << "done." << std::endl;

        if(args.jobs() == "all") {
            for(const auto& job_entry : all_jobs)
                jobs.push_back(job_entry.second);
        } else {
            const auto selected_jobs = SplitValueList(args.jobs(), false, ",", true);
            if(!selected_jobs.size())
                throw exception("No jobs selected.");
            for(const auto& job_name : selected_jobs) {
                if(!all_jobs.count(job_name))
                    throw exception("Job '%1%' not found.") % job_name;
                jobs.push_back(all_jobs.at(job_name));
            }
        }
    }

    void Run()
    {
        for(const auto& job : jobs) {
            try {
                std::cout << boost::format("Skimming job %1%...") % job.name << std::endl;
                ProcessJob(job);
                std::cout << "Job " << job.name << " has been skimmed." << std::endl;
            } catch(std::exception& e) {
                std::cerr << "\nEROOR: " << e.what() << "\nJob " <<job.name << " is failed." << std::endl;
            }
        }
    }

private:
    void ProcessJob(const SkimJob& job)
    {
        std::shared_ptr<std::thread> process_thread, writer_thread;
        std::shared_ptr<ProdSummary> summary;
        unsigned desc_id = 0;
        std::map<Channel,std::mt19937_64> gen_map;
        std::shared_ptr<std::uniform_int_distribution<unsigned int>> split_distr;
        if (setup.n_splits > 0)
            split_distr = std::make_shared<std::uniform_int_distribution<unsigned>>(0, setup.n_splits-1);

        try {
            weighting_mode = job.apply_common_weights ? job.weights | setup.common_weights : job.weights;
            for(auto desc_iter = job.files.begin(); desc_iter != job.files.end(); ++desc_iter, ++desc_id) {
                if(desc_iter == job.files.begin() || !job.ProduceMergedOutput()) {
                    for (Channel channel : setup.channels){
                        gen_map[channel].seed(setup.split_seed);
                    }
                    const std::string out_name = job.ProduceMergedOutput() ? job.merged_output : desc_iter->output;
                    outputFile = root_ext::CreateRootFile(args.outputPath() + "/" + out_name);

                    processQueue.SetAllDone(false);
                    writeQueue.SetAllDone(false);
                    process_thread = std::make_shared<std::thread>(std::bind(&TupleSkimmer::ProcessThread, this));
                    writer_thread = std::make_shared<std::thread>(std::bind(&TupleSkimmer::WriteThread, this));
                    summary = std::shared_ptr<ProdSummary>();
                }
                std::cout << "\tProcessing";
                std::vector<std::shared_ptr<TFile>> inputFiles;
                for(const auto& input : desc_iter->inputs) {
                    std::cout << " " << input;
                    inputFiles.push_back(root_ext::OpenRootFile(args.inputPath() + "/" + input));
                }
                std::cout << "\n\t\textracting summary" << std::endl;
                const ProdSummary desc_summary = GetCombinedSummary(*desc_iter, inputFiles);
                if(summary)
                    ntuple::MergeProdSummaries(*summary, desc_summary);
                else
                    summary = std::make_shared<ProdSummary>(desc_summary);

                summary->file_desc_name.push_back(desc_iter->inputs.at(0));
                summary->file_desc_id.push_back(desc_id);
                summary->n_splits = setup.n_splits;
                summary->split_seed = setup.split_seed;

                for(Channel channel : setup.channels) {
                    const std::string treeName = ToString(channel);
                    for(size_t n = 0; n < desc_iter->inputs.size(); ++n) {
                        auto file = inputFiles.at(n);
                        std::cout << "\t\t" << desc_iter->inputs.at(n) << ":" << treeName << std::endl;
                        std::shared_ptr<EventTuple> tuple;
                        try {
                            tuple = ntuple::CreateEventTuple(treeName, file.get(), true, ntuple::TreeState::Full);
                        } catch(std::exception&) {
                            std::cerr << "WARNING: tree " << treeName << " not found in file '"
                                      << desc_iter->inputs.at(n) << "'." << std::endl;
                        }
                        if(!tuple) continue;
                        for(const Event& event : *tuple) {
                            auto event_ptr = std::make_shared<Event>(event);
                            event_ptr->weight_xs = desc_iter->GetCrossSectionWeight();
                            event_ptr->file_desc_id = desc_id;
                            event_ptr->split_id = split_distr ? (*split_distr)(gen_map.at(channel)) : 0;
                            processQueue.Push(event_ptr);
                        }
                    }
                }

                if(std::next(desc_iter) == job.files.end() || !job.ProduceMergedOutput()) {
                    processQueue.SetAllDone(true);
                    process_thread->join();
                    writer_thread->join();

                    if(!summary)
                        throw exception("Summary not produced for job %1%") % job.name;
                    auto summaryTuple = ntuple::CreateSummaryTuple("summary", outputFile.get(), false,
                                                                   ntuple::TreeState::Skimmed);
                    (*summaryTuple)() = *summary;
                    summaryTuple->Fill();
                    summaryTuple->Write();
                }
            }
        } catch(std::exception&) {
            processQueue.SetAllDone(true);
            if(process_thread) process_thread->join();
            if(writer_thread) writer_thread->join();
            throw;
        }
    }

    void ProcessThread()
    {
        try {
            EventPtr event, prev_full_event;
            bool prev_full_event_stored = false;
            while(processQueue.Pop(event)) {
                ntuple::StorageMode storage_mode;
                const EventEnergyScale es = static_cast<EventEnergyScale>(event->eventEnergyScale);
                const bool store_event = ProcessEvent(*event, prev_full_event, storage_mode, prev_full_event_stored);
                if(store_event)
                    writeQueue.Push(event);
                if(storage_mode.IsFull() && es == EventEnergyScale::Central) {
                    prev_full_event = event;
                    prev_full_event_stored = store_event;
                }
            }
            writeQueue.SetAllDone();
        } catch(std::exception& e) {
            std::cerr << "ERROR (ProcessThread): " << e.what() << std::endl;
            std::abort();
        }
    }

    void WriteThread()
    {
        try {
            std::map<Channel, std::shared_ptr<EventTuple>> outputTuples;

            EventPtr event;
            while(writeQueue.Pop(event)) {
                const Channel channel = static_cast<Channel>(event->channelId);
                if(!outputTuples.count(channel)) {
                    const std::string treeName = ToString(channel);
                    outputTuples[channel] = ntuple::CreateEventTuple(treeName, outputFile.get(), false,
                                                                     ntuple::TreeState::Skimmed);
                }
                (*outputTuples[channel])() = *event;
                outputTuples[channel]->Fill();
            }

            for(auto& tuple : outputTuples)
                tuple.second->Write();
        } catch(std::exception& e) {
            std::cerr << "ERROR (ProcessThread): " << e.what() << std::endl;
            std::abort();
        }
    }

    ProdSummary GetSummaryWithWeights(const std::shared_ptr<TFile>& file, double xs_weight) const
    {
        using mc_corrections::WeightType;
        using mc_corrections::WeightingMode;

        static const WeightingMode shape_weights = { WeightType::PileUp, WeightType::BSM_to_SM, WeightType::DY,
                                                   WeightType::TTbar, WeightType::Wjets};
        static const WeightingMode shape_weights_withTopPt = shape_weights | WeightingMode({WeightType::TopPt});

        auto summary_tuple = ntuple::CreateSummaryTuple("summary", file.get(), true, ntuple::TreeState::Full);
        auto summary = ntuple::MergeSummaryTuple(*summary_tuple);
        summary.totalShapeWeight = 0;
        summary.totalShapeWeight_withTopPt = 0;

        const auto mode = shape_weights & weighting_mode;
        const auto mode_withTopPt = shape_weights_withTopPt & weighting_mode;
        const bool calc_withTopPt = mode_withTopPt.count(WeightType::TopPt);
        if(mode.size() || mode_withTopPt.size()) {
            ExpressTuple all_events("all_events", file.get(), true);
            for(const auto& event : all_events) {
                summary.totalShapeWeight += eventWeights_HH->GetTotalWeight(event, mode) * xs_weight;
                if(calc_withTopPt)
                    summary.totalShapeWeight_withTopPt +=
                            eventWeights_HH->GetTotalWeight(event, mode_withTopPt) * xs_weight;
            }
        }
        return summary;
    }

    ProdSummary GetCombinedSummary(const FileDescriptor& desc, const std::vector<std::shared_ptr<TFile>>& input_files)
    {
        if(!input_files.size())
            throw exception("Input files list is empty.");
        auto file_iter = input_files.begin();
        auto summary = GetSummaryWithWeights(*file_iter++, desc.GetCrossSectionWeight());
        if(!desc.first_input_is_ref) {
            for(; file_iter != input_files.end(); ++file_iter) {
                auto other_summary = GetSummaryWithWeights(*file_iter, desc.GetCrossSectionWeight());
                ntuple::MergeProdSummaries(summary, other_summary);
            }
        }
        return summary;
    }

    void SkimTauIds(std::vector<uint32_t>& tauId_keys, std::vector<float>& tauId_values) const
    {
        std::vector<uint32_t> skimmed_keys;
        std::vector<float> skimmed_values;
        if(tauId_keys.size() != tauId_values.size())
            throw exception("Inconsistent tauId information.");
        for(size_t n = 0; n < tauId_keys.size(); ++n) {
            if(setup.tau_id_hashes.count(tauId_keys.at(n))) {
                skimmed_keys.push_back(tauId_keys.at(n));
                skimmed_values.push_back(tauId_values.at(n));
            }
        }

        tauId_keys = std::move(skimmed_keys);
        tauId_values = std::move(skimmed_values);
    }

    bool ApplyTauIdCut(const std::vector<uint32_t>& tauId_keys, const std::vector<float>& tauId_values) const
    {
        for(size_t n = 0; n < tauId_keys.size(); ++n) {
            if(setup.tau_id_cut_hashes.count(tauId_keys.at(n)) &&
                    tauId_values.at(n) < setup.tau_id_cut_hashes.at(tauId_keys.at(n)))
                return false;
        }
        return true;
    }

    bool ProcessEvent(Event& event, const std::shared_ptr<Event>& prev_event, ntuple::StorageMode& storage_mode,
                      bool prev_event_stored)
    {
        using EventPart = ntuple::StorageMode::EventPart;
        using WeightType = mc_corrections::WeightType;

        Event full_event = event;
        storage_mode = ntuple::EventLoader::Load(full_event, prev_event.get());
        const EventEnergyScale es = static_cast<EventEnergyScale>(event.eventEnergyScale);
        if(!prev_event_stored) {
            event = full_event;
            storage_mode = ntuple::StorageMode::Full();
            event.storageMode = static_cast<UInt_t>(storage_mode.Mode());
        }

        if (!setup.energy_scales.count(es) || full_event.jets_p4.size() < 2 || full_event.extraelec_veto
                || full_event.extramuon_veto
                || std::abs(full_event.jets_p4.at(0).eta()) >= cuts::btag_2016::eta
                || std::abs(full_event.jets_p4.at(1).eta()) >= cuts::btag_2016::eta) return false;

        auto bb = full_event.jets_p4.at(0) + full_event.jets_p4.at(1);
        if (setup.apply_mass_cut
                && !cuts::hh_bbtautau_2016::hh_tag::IsInsideEllipse(full_event.SVfit_p4.mass(),bb.mass())) return false;

        if (setup.apply_charge_cut && (full_event.q_1+full_event.q_2) != 0) return false;

        if(!ApplyTauIdCut(full_event.tauId_keys_1, full_event.tauId_values_1) ||
                !ApplyTauIdCut(full_event.tauId_keys_2, full_event.tauId_values_2)) return false;

        if(storage_mode.IsPresent(EventPart::FirstTauIds))
            SkimTauIds(event.tauId_keys_1, event.tauId_values_1);

        if(storage_mode.IsPresent(EventPart::SecondTauIds))
            SkimTauIds(event.tauId_keys_2, event.tauId_values_2);
	
        event.n_jets = static_cast<unsigned>(full_event.jets_p4.size());
        event.ht_other_jets = static_cast<float>(
                    Calculate_HT(full_event.jets_p4.begin() + 2, full_event.jets_p4.end()));

        event.weight_pu = weighting_mode.count(WeightType::PileUp)
                        ? eventWeights_HH->GetWeight(full_event, WeightType::PileUp) : 1;
        if(weighting_mode.count(WeightType::LeptonTrigIdIso)) {
            auto lepton_wp = eventWeights_HH->GetProviderT<mc_corrections::LeptonWeights>(WeightType::LeptonTrigIdIso);
            event.weight_lepton_trig = lepton_wp->GetTriggerWeight(full_event);
            event.weight_lepton_id_iso = lepton_wp->GetIdIsoWeight(full_event);
        } else {
            event.weight_lepton_trig = 1;
            event.weight_lepton_id_iso = 1;
        }
        if(weighting_mode.count(WeightType::BTag)) {
            auto btag_wp = eventWeights_HH->GetProviderT<mc_corrections::BTagWeight>(WeightType::BTag);
            event.weight_btag = btag_wp->Get(full_event);
            event.weight_btag_up = btag_wp->GetEx(full_event, UncertaintyScale::Up);
            event.weight_btag_down = btag_wp->GetEx(full_event, UncertaintyScale::Down);
        } else {
            event.weight_btag = 1;
            event.weight_btag_up = 1;
            event.weight_btag_down = 1;
        }
        event.weight_dy = weighting_mode.count(WeightType::DY)
                ? eventWeights_HH->GetWeight(full_event, WeightType::DY) : 1;
        event.weight_ttbar = weighting_mode.count(WeightType::TTbar)
                ? eventWeights_HH->GetWeight(full_event, WeightType::TTbar) : 1;
        event.weight_wjets = weighting_mode.count(WeightType::Wjets)
                ? eventWeights_HH->GetWeight(full_event, WeightType::Wjets) : 1;
        event.weight_bsm_to_sm = weighting_mode.count(WeightType::BSM_to_SM)
                ? eventWeights_HH->GetWeight(full_event, WeightType::BSM_to_SM) : 1;
        event.weight_top_pt = weighting_mode.count(WeightType::TopPt)
                ? eventWeights_HH->GetWeight(full_event, WeightType::TopPt) : 1;
        event.weight_total = eventWeights_HH->GetTotalWeight(full_event, weighting_mode) * event.weight_xs;

        if(storage_mode.IsPresent(EventPart::Jets)) {
            event.jets_csv.resize(2);
            event.jets_rawf.resize(2);
            event.jets_mva.resize(2);
            event.jets_p4.resize(2);
            event.jets_partonFlavour.resize(2);
            event.jets_hadronFlavour.resize(2);
        }
		
        return true;
    }

private:
    Arguments args;
    Setup setup;
    std::vector<SkimJob> jobs;
    EventQueue processQueue, writeQueue;
    std::shared_ptr<mc_corrections::EventWeights_HH> eventWeights_HH;
	std::shared_ptr<TFile> outputFile;
    mc_corrections::WeightingMode weighting_mode;


};

} // namespace tuple_skimmer
} // namespace analysis

PROGRAM_MAIN(analysis::tuple_skimmer::TupleSkimmer, Arguments)
