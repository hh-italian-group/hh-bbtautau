/*! Skim EventTuple.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <thread>
#include <functional>
#include <random>

#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Run/include/EntryQueue.h"
#include "AnalysisTools/Core/include/ProgressReporter.h"
#include "h-tautau/McCorrections/include/EventWeights.h"
#include "hh-bbtautau/McCorrections/include/EventWeights_HH.h"
#include "hh-bbtautau/Instruments/include/SkimmerConfigEntryReader.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "h-tautau/Analysis/include/EventLoader.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"
#include "hh-bbtautau/Analysis/include/AnalysisCategories.h"
#include "h-tautau/Analysis/include/MetFilters.h"

#include "h-tautau/McCorrections/include/PileUpWeight.h"
#include "h-tautau/McCorrections/include/LeptonWeights.h"
#include "h-tautau/McCorrections/include/BTagWeight.h"
#include "h-tautau/McCorrections/include/TopPtWeight.h"
#include "h-tautau/McCorrections/include/TauIdWeight.h"
#include "h-tautau/McCorrections/include/GenEventWeight.h"
#include "hh-bbtautau/McCorrections/include/HH_nonResonant_weight.h"
#include "h-tautau/Analysis/include/SignalObjectSelector.h"

struct Arguments {
    REQ_ARG(std::string, cfg);
    REQ_ARG(std::string, inputPath);
    REQ_ARG(std::string, outputPath);
    REQ_ARG(std::string, jobs);
    REQ_ARG(std::string, setup_name);
    OPT_ARG(bool, use_LLR_weights, false);
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
    using Filter = ntuple::MetFilters::Filter;

    static constexpr size_t max_queue_size = 100000;

    TupleSkimmer(const Arguments& _args) :
        args(_args), processQueue(max_queue_size), writeQueue(max_queue_size), signalObjectSelector(setup.mode)
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

        std::cout << "done.\nLoading weights... " << std::flush;
        eventWeights_HH = std::make_shared<mc_corrections::EventWeights_HH>(setup.period, setup.jet_ordering, setup.btag_wp,
                                                                            args.use_LLR_weights(),setup.applyTauId);
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
        std::set<std::string> successful_jobs, failed_jobs;
        for(const auto& job : jobs) {
            try {
                std::cout << boost::format("Skimming job %1%...") % job.name << std::endl;
                ProcessJob(job);
                std::cout << "Job " << job.name << " has been skimmed." << std::endl;
                successful_jobs.insert(job.name);
            } catch(std::exception& e) {
                std::cerr << "\nERROR: " << e.what() << "\nJob " <<job.name << " is failed." << std::endl;
                failed_jobs.insert(job.name);
            }
        }
        std::cout << "Summary:";
        if(failed_jobs.empty())
            std::cout << " all jobs has been successfully skimmed.";
        std::cout << std::endl;

        auto print_jobs = [](const std::set<std::string>& job_set, const std::string& name) {
            if(!job_set.empty()) {
                std::cout << name << " jobs:";
                for(const auto& job : job_set)
                    std::cout << " " << job;
                std::cout << std::endl;
            }
        };

        print_jobs(successful_jobs, "Successful");
        print_jobs(failed_jobs, "Failed");
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
                    outputFile = root_ext::CreateRootFile(args.outputPath() + "/" + out_name, ROOT::kLZ4, 4);

                    processQueue.SetAllDone(false);
                    writeQueue.SetAllDone(false);
                    process_thread = std::make_shared<std::thread>(std::bind(&TupleSkimmer::ProcessThread, this));
                    writer_thread = std::make_shared<std::thread>(std::bind(&TupleSkimmer::WriteThread, this));
                    summary = std::shared_ptr<ProdSummary>();

                    if(weighting_mode.count(mc_corrections::WeightType::BSM_to_SM)) {
                        std::cout << "\tPreparing EFT weights" << std::endl;
                        auto eft_weight_provider = eventWeights_HH->GetProviderT<NonResHH_EFT::WeightProvider>(
                                    mc_corrections::WeightType::BSM_to_SM);
                        auto express_tuple = ntuple::CreateExpressTuple("all_events", outputFile.get(), false,
                                                                       ntuple::TreeState::Full);

                        for(auto desc_iter_2 = job.files.begin(); desc_iter_2 != job.files.end(); ++desc_iter_2) {
                            if(desc_iter_2 != desc_iter && !job.ProduceMergedOutput()) continue;
                            for(const auto& input : desc_iter_2->inputs) {
                                auto file = root_ext::OpenRootFile(args.inputPath() + "/" + input);
                                eft_weight_provider->AddFile(*file);
                                ntuple::ExpressTuple file_events("all_events", file.get(), true);

                                using EventIdSet = std::set<EventIdentifier>;
                                EventIdSet processed_events;
                                for(const auto& event : file_events) {
                                    (*express_tuple)() = event;
                                    const EventIdentifier Id(event.run, event.lumi, event.evt);
                                    if(processed_events.count(Id)) {
                                        std::cout << "WARNING: duplicated express event " << Id << std::endl;
                                        continue;
                                    }
                                    processed_events.insert(Id);
                                    express_tuple->Fill();
                                }
                            }
                        }
                        express_tuple->Write();
                        eft_weight_provider->CreatePdfs(outputFile.get());
                    }
                }
                std::cout << "\tProcessing";
                std::vector<std::shared_ptr<TFile>> inputFiles;
                for(const auto& input : desc_iter->inputs) {
                    std::cout << " " << input;
                    inputFiles.push_back(root_ext::OpenRootFile(args.inputPath() + "/" + input));
                }
                std::cout << "\n\t\textracting summary" << std::endl;
                double weight_xs, weight_xs_withTopPt;
                const ProdSummary desc_summary = GetCombinedSummary(*desc_iter, inputFiles, weight_xs,
                                                                    weight_xs_withTopPt);
                if(summary)
                    ntuple::MergeProdSummaries(*summary, desc_summary);
                else
                    summary = std::make_shared<ProdSummary>(desc_summary);

//                for(unsigned n = 0; n < desc_iter->inputs.size() && (n == 0 || !desc_iter->first_input_is_ref); ++n) {
//                    if (desc_iter->input_is_partial.size() && desc_iter->input_is_partial.at(n) == true) continue;
//                    const auto& input = desc_iter->inputs.at(n);
//                    summary->file_desc_name.push_back(input);
//                    summary->file_desc_id.push_back(n * 1000 + desc_id);
//                }

                summary->n_splits = setup.n_splits;
                summary->split_seed = setup.split_seed;

                for(Channel channel : setup.channels) {
                    const std::string treeName = ToString(channel);

                    using FullEventId = std::tuple<EventIdentifier, EventEnergyScale>;
                    using FullEventIdSet = std::set<FullEventId>;

                    FullEventIdSet processed_events;
                    for(size_t n = 0; n < desc_iter->inputs.size(); ++n) {
                        auto file = inputFiles.at(n);
                        std::cout << "\t\t" << desc_iter->inputs.at(n) << ":" << treeName << std::endl;

                        if(!desc_iter->first_input_is_ref && (!desc_iter->input_is_partial.size() || desc_iter->input_is_partial.at(n) == false))
                            processed_events.clear();
                        if((n == 0 || !desc_iter->first_input_is_ref)
                                && (desc_iter->input_is_partial.empty() || !desc_iter->input_is_partial.at(n))
                                && weighting_mode.count(mc_corrections::WeightType::PileUp)
                                && setup.period == Period::Run2017) {
                            auto pile_up_weight = eventWeights_HH->GetProviderT<mc_corrections::PileUpWeightEx>(mc_corrections::WeightType::PileUp);

                            auto dataset_name = RemoveFileExtension(desc_iter->inputs.at(n));
                            pile_up_weight->SetActiveDataset(dataset_name);
                        }

                        std::shared_ptr<EventTuple> tuple;
                        try {
                            tuple = ntuple::CreateEventTuple(treeName, file.get(), true, ntuple::TreeState::Full);
                        } catch(std::exception&) {
                            std::cerr << "WARNING: tree " << treeName << " not found in file '"
                                      << desc_iter->inputs.at(n) << "'." << std::endl;
                        }
                        if(!tuple) continue;
                        EventIdentifier prev_event_id = EventIdentifier::Undef_event();
                        unsigned split_id = 0;

                        for(const Event& event : *tuple) {
                            const FullEventId fullId{EventIdentifier(event.run, event.lumi, event.evt),
                                                     static_cast<EventEnergyScale>(event.eventEnergyScale)};
                            if(processed_events.count(fullId)) {
                                std::cout << "WARNING: duplicated event " << std::get<EventIdentifier>(fullId) << " "
                                          << std::get<EventEnergyScale>(fullId) << std::endl;
                                continue;
                            }
                            processed_events.insert(fullId);

                            auto event_ptr = std::make_shared<Event>(event);
                            event_ptr->weight_xs = weight_xs;
                            event_ptr->weight_xs_withTopPt = weight_xs_withTopPt;
                            event_ptr->file_desc_id = desc_id;
                            event_ptr->split_id = 0;
                            event_ptr->isData = job.isData;
                            if(split_distr) {
                                const EventIdentifier event_id(event);
                                event_ptr->split_id = event_id == prev_event_id
                                                    ? split_id : (*split_distr)(gen_map.at(channel));
                                prev_event_id = event_id;
                                split_id = event_ptr->split_id;
                            }
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


    ntuple::ProdSummary GetSummaryWithWeights(const FileDescriptor& desc, const std::shared_ptr<TFile>& file,
                                              size_t file_index)
    {
        if(weighting_mode.count(mc_corrections::WeightType::PileUp) && setup.period == Period::Run2017){
            auto pile_up_weight = eventWeights_HH->GetProviderT<mc_corrections::PileUpWeightEx>(
                                                                mc_corrections::WeightType::PileUp);
            auto dataset_name = RemoveFileExtension(desc.inputs.at(file_index));
            pile_up_weight->SetActiveDataset(dataset_name);
        }
        return eventWeights_HH->GetSummaryWithWeights(file, weighting_mode);
    }

    ProdSummary GetCombinedSummary(const FileDescriptor& desc, const std::vector<std::shared_ptr<TFile>>& input_files,
                                   double& weight_xs, double& weight_xs_withTopPt)
    {
        if(!input_files.size())
            throw exception("Input files list is empty.");
        auto file_iter = input_files.begin();
        auto summary = GetSummaryWithWeights(desc, *file_iter++, 0);
        if(!desc.first_input_is_ref) {
            unsigned n=1;
            for(; file_iter != input_files.end(); ++file_iter, ++n) {
                if (desc.input_is_partial.size() && desc.input_is_partial.at(n) == true)
                    continue;
                auto other_summary = GetSummaryWithWeights(desc, *file_iter, n);
                ntuple::MergeProdSummaries(summary, other_summary);
            }
        }
        if(desc.HasCrossSection()) {
            weight_xs = desc.GetCrossSectionWeight() / summary.totalShapeWeight;
            summary.totalShapeWeight = desc.GetCrossSectionWeight();
            weight_xs_withTopPt = desc.GetCrossSectionWeight() / summary.totalShapeWeight_withTopPt;
            summary.totalShapeWeight_withTopPt = desc.GetCrossSectionWeight();
        } else {
            weight_xs = 1;
            weight_xs_withTopPt = 1;
        }
        return summary;
    }

    // bool ApplyTauIdCut(TauIdResults::BitsContainer id_bits) const
    // {
    //     const TauIdResults id_results(id_bits);
    //     for(size_t index : setup.tau_id_cut_indices) {
    //         if(!id_results.Result(index))
    //             return false;
    //     }
    //     return true;
    // }

    bool ProcessEvent(Event& event, const std::shared_ptr<Event>& prev_event, ntuple::StorageMode& storage_mode,
                      bool prev_event_stored)
    {
        // using EventPart = ntuple::StorageMode::EventPart;
        using WeightType = mc_corrections::WeightType;
        using WeightingMode = mc_corrections::WeightingMode;

        Event full_event = event;

        storage_mode = ntuple::EventLoader::Load(full_event, prev_event.get());

        boost::optional<EventInfoBase> eventInfo = CreateEventInfo(full_event,signalObjectSelector,nullptr,analysis::TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017,setup.period,setup.jet_ordering);
        if(!eventInfo.is_initialized()) return false;

        const EventEnergyScale es = static_cast<EventEnergyScale>(event.eventEnergyScale);
        if(!prev_event_stored) {
            event = full_event;
            storage_mode = ntuple::StorageMode::Full();
            event.storageMode = static_cast<UInt_t>(storage_mode.Mode());
        }

        auto event_metFilters = ntuple::MetFilters(full_event.metFilters);
        if (!event_metFilters.Pass(Filter::PrimaryVertex)  || !event_metFilters.Pass(Filter::BeamHalo) ||
            !event_metFilters.Pass(Filter::HBHE_noise)  || !event_metFilters.Pass(Filter::HBHEiso_noise) ||
            !event_metFilters.Pass(Filter::ECAL_TP) || !event_metFilters.Pass(Filter::badMuon) ||
            !event_metFilters.Pass(Filter::ecalBadCalib)) return false;

        if(event.isData && !event_metFilters.Pass(Filter::ee_badSC_noise)) return false;

        if(!setup.energy_scales.count(es) || full_event.extraelec_veto || full_event.extramuon_veto) return false;

        if(setup.apply_bb_cut && !eventInfo->HasBjetPair()) return false;

        if(setup.apply_mass_cut) {
            bool pass_mass_cut = false;
            if (!eventInfo->HasBjetPair()) return false;
            const double mbb = eventInfo->GetHiggsBB().GetMomentum().mass();
            const double mtautau = (eventInfo->GetLeg(1).GetMomentum() + eventInfo->GetLeg(2).GetMomentum()).mass();
            pass_mass_cut = pass_mass_cut
                    || cuts::hh_bbtautau_2016::hh_tag::IsInsideBoostedMassWindow(eventInfo->GetSVFitResults().momentum.mass(),mbb);

            if(setup.massWindowParams.count(SelectionCut::mh))
                pass_mass_cut = pass_mass_cut || setup.massWindowParams.at(SelectionCut::mh)
                        .IsInside(eventInfo->GetSVFitResults().momentum.mass(),mbb);

            if(setup.massWindowParams.count(SelectionCut::mhVis))
                pass_mass_cut = pass_mass_cut || setup.massWindowParams.at(SelectionCut::mhVis)
                        .IsInside(mtautau,mbb);

            if(setup.massWindowParams.count(SelectionCut::mhMET))
                pass_mass_cut = pass_mass_cut || setup.massWindowParams.at(SelectionCut::mhMET)
                        .IsInside((eventInfo->GetLeg(1).GetMomentum() + eventInfo->GetLeg(2).GetMomentum() + full_event.pfMET_p4).mass(),mbb);

            if(!pass_mass_cut) return false;
        }

        if (setup.apply_charge_cut && (eventInfo->GetLeg(1)->charge() + eventInfo->GetLeg(2)->charge()) != 0) return false;

        // if(setup.ApplyTauIdCuts()) {
        //     const Channel channel = static_cast<Channel>(event.channelId);
        //     const auto leg_types = GetChannelLegTypes(channel);
        //     if(leg_types.first == LegType::tau && !ApplyTauIdCut(full_event.tauId_flags_1)) return false;
        //     if(leg_types.second == LegType::tau && !ApplyTauIdCut(full_event.tauId_flags_2)) return false;
        // }

        if(setup.apply_kinfit){
            event.kinFit_chi2.push_back(static_cast<Float_t>(eventInfo->GetKinFitResults().chi2));
            event.kinFit_convergence.push_back(eventInfo->GetKinFitResults().convergence);
            event.kinFit_m.push_back(static_cast<Float_t>(eventInfo->GetKinFitResults().mass));
            event.kinFit_jetPairId.push_back(static_cast<unsigned>(ntuple::CombinationPairToIndex(eventInfo->GetSelectedSignalJets().selectedBjetPair, eventInfo->GetNJets())));
        }

        event.ht_other_jets = (eventInfo->HasBjetPair()) ? static_cast<Float_t>(eventInfo->GetHT(false,true)) : 0;

        event.weight_pu = weighting_mode.count(WeightType::PileUp)
                        ? eventWeights_HH->GetWeight(*eventInfo, WeightType::PileUp) : 1;
        if(weighting_mode.count(WeightType::LeptonTrigIdIso)) {
            auto lepton_wp = eventWeights_HH->GetProviderT<mc_corrections::LeptonWeights>(WeightType::LeptonTrigIdIso);
            event.weight_lepton_trig = lepton_wp->GetTriggerWeight(*eventInfo);
            event.weight_lepton_id_iso = lepton_wp->GetIdIsoWeight(*eventInfo);
        } else {
            event.weight_lepton_trig = 1;
            event.weight_lepton_id_iso = 1;
        }
        event.weight_tau_id = weighting_mode.count(WeightType::TauId)
                ? eventWeights_HH->GetWeight(*eventInfo, WeightType::TauId) : 1;
        if(weighting_mode.count(WeightType::BTag)) {
            auto btag_wp = eventWeights_HH->GetProviderT<mc_corrections::BTagWeight>(WeightType::BTag);
            event.weight_btag = btag_wp->Get(*eventInfo);
            event.weight_btag_up = btag_wp->GetEx(*eventInfo, UncertaintyScale::Up);
            event.weight_btag_down = btag_wp->GetEx(*eventInfo, UncertaintyScale::Down);
        } else {
            event.weight_btag = 1;
            event.weight_btag_up = 1;
            event.weight_btag_down = 1;
        }
        event.weight_dy = weighting_mode.count(WeightType::DY)
                ? eventWeights_HH->GetWeight(*eventInfo, WeightType::DY) : 1;
        event.weight_ttbar = weighting_mode.count(WeightType::TTbar)
                ? eventWeights_HH->GetWeight(*eventInfo, WeightType::TTbar) : 1;
        event.weight_wjets = weighting_mode.count(WeightType::Wjets)
                ? eventWeights_HH->GetWeight(*eventInfo, WeightType::Wjets) : 1;
        event.weight_bsm_to_sm = weighting_mode.count(WeightType::BSM_to_SM)
                ? eventWeights_HH->GetWeight(*eventInfo, WeightType::BSM_to_SM) : 1;

        if(weighting_mode.count(WeightType::TopPt)) {
            event.weight_top_pt = eventWeights_HH->GetWeight(*eventInfo, WeightType::TopPt);
            const auto wmode_withoutTopPt = weighting_mode - WeightingMode{WeightType::TopPt};
            event.weight_total = eventWeights_HH->GetTotalWeight(*eventInfo, wmode_withoutTopPt) * event.weight_xs;
            event.weight_total_withTopPt = eventWeights_HH->GetTotalWeight(*eventInfo, weighting_mode)
                    * event.weight_xs_withTopPt;
        } else {
            event.weight_top_pt = 1;
            event.weight_total = eventWeights_HH->GetTotalWeight(*eventInfo, weighting_mode) * event.weight_xs;
            event.weight_total_withTopPt = 0;
        }

//        if(setup.apply_bb_cut && storage_mode.IsPresent(EventPart::Jets)) {
//            event.jets_csv.resize(2);
//            event.jets_rawf.resize(2);
//            event.jets_pu_id.resize(2);
//            event.jets_p4.resize(2);
//            event.jets_hadronFlavour.resize(2);
//        }

        if(!setup.keep_genJets) {
            event.genJets_p4.clear();
            event.genJets_hadronFlavour.clear();
        }

        if(!setup.keep_genParticles){
            event.genParticles_p4.clear();
            event.genParticles_pdg.clear();
        }

        if(!setup.keep_MET_cov) {
            for(unsigned int i = 0; i < static_cast<unsigned>(event.pfMET_cov.kRows); i++) {
                for(int j = 0; j < event.pfMET_cov.kCols; j++) {
                    event.pfMET_cov[i][j] = 0;
                }
            }
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
    SignalObjectSelector signalObjectSelector;
};

} // namespace tuple_skimmer
} // namespace analysis

PROGRAM_MAIN(analysis::tuple_skimmer::TupleSkimmer, Arguments)
