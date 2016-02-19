/*! Sync Tree Producer from FlatTree.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/LightBaseFlatTreeAnalyzer.h"
#include "h-tautau/Analysis/include/SyncTree.h"

class SyncTreeProducer : public analysis::LightBaseFlatTreeAnalyzer {
public:
    typedef std::map<analysis::EventEnergyScale, ntuple::Flat> ES_toEvent_Map;
    typedef std::map<analysis::EventId, ES_toEvent_Map> EventId_ToES_Map;

    SyncTreeProducer(const std::string& inputFileName, const std::string& outputFileName)
         : LightBaseFlatTreeAnalyzer(inputFileName, outputFileName), inclusive(0), passed(0)
    {
        syncTree = std::shared_ptr<ntuple::SyncTree>(new ntuple::SyncTree("syncTree", GetOutputFile().get(), false));
        recalc_kinfit = false;
    }

    virtual ~SyncTreeProducer() { syncTree->Write(); }

protected:
    static bool PassSyncTreeSelection(const analysis::FlatEventInfo& eventInfo)
    {
        using analysis::EventRegion;
        const ntuple::Flat& event = *eventInfo.event;
        if (eventInfo.channel == analysis::Channel::MuTau){
            using namespace cuts::Htautau_Summer13::MuTau;

            return !(!event.againstMuonTight_2
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                || event.pfRelIso_1 >= muonID::pFRelIso || event.q_1 * event.q_2 == +1);
        }
        if (eventInfo.channel == analysis::Channel::ETau){
            using namespace cuts::Htautau_Summer13::ETau;

            return !(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                || event.pfRelIso_1 >= electronID::pFRelIso || event.q_1 * event.q_2 == +1);
        }
        if (eventInfo.channel == analysis::Channel::TauTau){
            using namespace cuts::Htautau_Summer13::TauTau::tauID;

//            return !(event.againstElectronLooseMVA_2 <= againstElectronLooseMVA3
//                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= byCombinedIsolationDeltaBetaCorrRaw3Hits
//                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= byCombinedIsolationDeltaBetaCorrRaw3Hits);
            if(!event.againstElectronLooseMVA_2
                    || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= BackgroundEstimation::Isolation_upperLimit
                    || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= BackgroundEstimation::Isolation_upperLimit
                    || (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= byCombinedIsolationDeltaBetaCorrRaw3Hits
                        && event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= byCombinedIsolationDeltaBetaCorrRaw3Hits))
                return false;

            const bool os = event.q_1 * event.q_2 == -1;
            const bool iso = event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                             event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < byCombinedIsolationDeltaBetaCorrRaw3Hits;

            if(eventInfo.bjet_momentums.size() < 2 || !eventInfo.fitResults.has_valid_mass)
                return false;

            using namespace cuts::massWindow;
            const double mass_tautau = eventInfo.event->m_sv_MC;
            const double mass_bb = eventInfo.Hbb.M();
            const bool inside_mass_window = mass_tautau > m_tautau_low && mass_tautau < m_tautau_high
                    && mass_bb > m_bb_low && mass_bb < m_bb_high;

            return os && iso && inside_mass_window;
        }
        throw analysis::exception("unsupported channel ") << eventInfo.channel;
    }

    virtual void AnalyzeEvent(const analysis::FlatEventInfo& eventInfo, analysis::EventCategory category) override
    {
        using analysis::EventCategory;

        if (category != EventCategory::Inclusive) return;
        ++inclusive;

        const ntuple::Flat& event = *eventInfo.event;
        if(event.evt == 2823968 && event.eventEnergyScale == 0){
            std::cout << "event: " << event.evt << std::endl;
            std::cout << "b1: " << eventInfo.bjet_momentums.at(0)<< ", csv1: " <<
                         event.csv_Bjets.at(eventInfo.selected_bjets.first) << std::endl;
            std::cout << "b2: " << eventInfo.bjet_momentums.at(1) << ", csv2: " <<
                         event.csv_Bjets.at(eventInfo.selected_bjets.second) << std::endl;
            std::cout << "msv= " << eventInfo.event->m_sv_MC << std::endl;
        }

        if (!PassSyncTreeSelection(eventInfo)) return;
        ++passed;
        if (eventInfo.eventType != ntuple::EventType::ZTT) return;
        ++passedEventType;


        const analysis::EventId eventId(event.run,event.lumi,event.evt);
        eventId_ToES_Map[eventId][static_cast<analysis::EventEnergyScale>(event.eventEnergyScale)] = event;
    }

    virtual void EndOfRun() override
    {
        std::cout << "inclusive evt = " << inclusive << ", passed = " << passed << ", passedEventType = " <<
                     passedEventType << std::endl;
        std::cout << "eventId_ToES_Map.size = " << eventId_ToES_Map.size() << std::endl;
        for (const auto& event_iter : eventId_ToES_Map){
            const ES_toEvent_Map& es_toEventMap = event_iter.second;

            if (!es_toEventMap.count(analysis::EventEnergyScale::Central)) continue;
            const ntuple::Flat& event = es_toEventMap.at(analysis::EventEnergyScale::Central);
            analysis::FlatEventInfo eventInfo(event, analysis::FlatEventInfo::BjetPair(0, 1), false);
            static const double default_value = ntuple::DefaultFillValueForSyncTree();
            syncTree->run() = event.run;
            syncTree->lumi() = event.lumi;
            syncTree->evt() = event.evt;

            syncTree->npv() = event.npv;
            syncTree->npu() = event.npu;

            syncTree->puweight() = event.puweight;
            syncTree->trigweight_1() = event.trigweight_1;
            syncTree->trigweight_2() = event.trigweight_2;
            syncTree->idweight_1() = event.idweight_1;
            syncTree->idweight_2() = event.idweight_2;
            syncTree->isoweight_1() = event.isoweight_1;
            syncTree->isoweight_2() = event.isoweight_2;
            syncTree->fakeweight() = event.fakeweight_2;
            double DYweight = 1;
            //inclusive sample
//            if (event.n_extraJets_MC == 6)
//                DYweight = 0.1941324;
//            if (event.n_extraJets_MC == 7)
//                DYweight = 0.0787854;
//            if (event.n_extraJets_MC == 8)
//                DYweight = 0.0457089;
//            if (event.n_extraJets_MC == 9)
//                DYweight = 0.0357635;
            //exclusive sample
//            DYweight = 0.1941324; //1jet
//            DYweight = 0.0787854; //2jets
//            DYweight = 0.0457089; //3jets
//            DYweight = 0.0357635; //4jets
            syncTree->DYweight() = DYweight;
//            double correct_weight = event.weight * DYweight / event.fakeweight_2;
//            if (event.decayMode_2 == ntuple::tau_id::kOneProng0PiZero)
//                correct_weight = correct_weight/0.88;
//            double fakeWeight =
//                    cuts::Htautau_Summer13::electronEtoTauFakeRateWeight::CalculateEtoTauFakeWeight(
//                        event.eta_2,
//                        ntuple::tau_id::ConvertToHadronicDecayMode(event.decayMode_2));
//            syncTree->etau_fakerate() = fakeWeight;
//            syncTree->weight() = correct_weight * fakeWeight;
            //without DYweight and decayMode weight correction
            syncTree->weight() = event.weight * DYweight;
            syncTree->embeddedWeight() = event.embeddedWeight;
            syncTree->decayModeWeight_1() = event.decayMode_1 == ntuple::tau_id::kOneProng0PiZero
                    ? cuts::Htautau_Summer13::tauCorrections::DecayModeWeight : 1;
            syncTree->decayModeWeight_2() = event.decayMode_2 == ntuple::tau_id::kOneProng0PiZero
                    ? cuts::Htautau_Summer13::tauCorrections::DecayModeWeight : 1;

            syncTree->mvis() = eventInfo.Htt.M();
            syncTree->m_sv() = event.m_sv_MC;
            syncTree->pt_sv() = event.pt_sv_MC;
            syncTree->eta_sv() = event.eta_sv_MC;
            syncTree->phi_sv() = event.phi_sv_MC;
            if (es_toEventMap.count(analysis::EventEnergyScale::TauUp)){
                const ntuple::Flat& event_up = es_toEventMap.at(analysis::EventEnergyScale::TauUp);
                syncTree->m_sv_Up() = event_up.m_sv_MC;
            }
            else
                syncTree->m_sv_Up() = default_value;
            if (es_toEventMap.count(analysis::EventEnergyScale::TauDown)){
                const ntuple::Flat& event_down = es_toEventMap.at(analysis::EventEnergyScale::TauDown);
                syncTree->m_sv_Down() = event_down.m_sv_MC;
            }
            else
                syncTree->m_sv_Down() = default_value;

            syncTree->pt_1() = event.pt_1;
            syncTree->phi_1() = event.phi_1;
            syncTree->eta_1() = event.eta_1;
            syncTree->m_1() = event.m_1;
            syncTree->q_1() = event.q_1;
            syncTree->mt_1() = event.mt_1;
            syncTree->d0_1() = event.d0_1;
            syncTree->dZ_1() = event.dZ_1;

            // MVA iso for hadronic Tau, Delta Beta for muon and electron
            syncTree->iso_1() = event.iso_1;
            syncTree->mva_1() = event.mva_1;
            syncTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1() = event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
            syncTree->againstElectronMVA3raw_1() = event.againstElectronMVA3raw_1;
            syncTree->byIsolationMVA2raw_1() = event.byIsolationMVA2raw_1;
            syncTree->againstMuonLoose2_1() = event.againstMuonLoose_1;
            syncTree->againstMuonMedium2_1() = event.againstMuonMedium_1;
            syncTree->againstMuonTight2_1() = event.againstMuonTight_1;

            syncTree->pt_2() = event.pt_2;
            syncTree->phi_2() = event.phi_2;
            syncTree->eta_2() = event.eta_2;
            syncTree->m_2() = event.m_2;
            syncTree->q_2() = event.q_2;
            syncTree->mt_2() = event.mt_2;
            syncTree->d0_2() = event.d0_2;
            syncTree->dZ_2() = event.dZ_2;

            syncTree->iso_2() = event.iso_2;
            syncTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
            syncTree->againstElectronMVA3raw_2() = event.againstElectronMVA3raw_2;
            syncTree->byIsolationMVA2raw_2() = event.byIsolationMVA2raw_2;
            syncTree->againstMuonLoose2_2() = event.againstMuonLoose_2;
            syncTree->againstMuonMedium2_2() = event.againstMuonMedium_2;
            syncTree->againstMuonTight2_2() = event.againstMuonTight_2;

            syncTree->pt_tt() = eventInfo.Htt_MET.Pt();

            syncTree->met() = event.met;
            syncTree->metphi() = event.metphi;
            syncTree->metcov00() = event.metcov00;
            syncTree->metcov01() = event.metcov01;
            syncTree->metcov10() = event.metcov10;
            syncTree->metcov11() = event.metcov11;

            syncTree->mvamet() = event.mvamet;
            syncTree->mvametphi() = event.mvametphi;
            syncTree->mvacov00() = event.mvacov00;
            syncTree->mvacov01() = event.mvacov01;
            syncTree->mvacov10() = event.mvacov10;
            syncTree->mvacov11() = event.mvacov11;

            syncTree->njets() = event.njets;
            syncTree->njetspt20() = event.njetspt20;

            syncTree->jpt_1() = default_value;
            syncTree->jeta_1() = default_value;
            syncTree->jphi_1() = default_value;
            syncTree->jptraw_1() = default_value;
            //syncTree->jptunc_1() = default_value;
            syncTree->jmva_1() = default_value;
            syncTree->jpass_1() = default_value;
            syncTree->jpt_2() = default_value;
            syncTree->jeta_2() = default_value;
            syncTree->jphi_2() = default_value;
            syncTree->jptraw_2() = default_value;
    //            syncTree->jptunc_2() = default_value;
            syncTree->jmva_2() = default_value;
            syncTree->jpass_2() = default_value;

            if (event.njets >= 1) {
                syncTree->jpt_1() = event.pt_jets.at(0);
                syncTree->jeta_1() = event.eta_jets.at(0);
                syncTree->jphi_1() = event.phi_jets.at(0);
                syncTree->jptraw_1() = event.ptraw_jets.at(0);
                //syncTree->jptunc_1();
                syncTree->jmva_1() = event.mva_jets.at(0);
                syncTree->jpass_1() = event.passPU_jets.at(0);
            }

            if (event.njets >= 2) {
                syncTree->jpt_2() = event.pt_jets.at(1);
                syncTree->jeta_2() = event.eta_jets.at(1);
                syncTree->jphi_2() = event.phi_jets.at(1);
                syncTree->jptraw_2() = event.ptraw_jets.at(1);
                //syncTree->jptunc_2();
                syncTree->jmva_2() = event.mva_jets.at(1);
                syncTree->jpass_2() = event.passPU_jets.at(1);
            }

            syncTree->nbtag() = event.nBjets_retagged;
            //syncTree->nbtag() = event.nBjets;

            syncTree->bpt_1() = default_value;
            syncTree->beta_1() = default_value;
            syncTree->bphi_1() = default_value;
            syncTree->bcsv_1() = default_value;
            syncTree->bpt_2() = default_value;
            syncTree->beta_2() = default_value;
            syncTree->bphi_2() = default_value;
            syncTree->bcsv_2() = default_value;
            syncTree->m_bb() = default_value;
            syncTree->m_ttbb() = default_value;
            syncTree->bpt_3() = default_value;
            syncTree->beta_3() = default_value;
            syncTree->bphi_3() = default_value;
            syncTree->bcsv_3() = default_value;

            if (event.nBjets >= 1 && event.csv_Bjets.at(0) > cuts::Htautau_Summer13::btag::CSVM) {
                syncTree->bpt_1() = event.pt_Bjets.at(0);
                syncTree->beta_1() = event.eta_Bjets.at(0);
                syncTree->bphi_1() = event.phi_Bjets.at(0);
                syncTree->bcsv_1() = event.csv_Bjets.at(0);
            }

            if (event.nBjets >= 2 && event.csv_Bjets.at(0) > cuts::Htautau_Summer13::btag::CSVM &&
                    event.csv_Bjets.at(1) > cuts::Htautau_Summer13::btag::CSVM) {
                syncTree->bpt_2() = event.pt_Bjets.at(1);
                syncTree->beta_2() = event.eta_Bjets.at(1);
                syncTree->bphi_2() = event.phi_Bjets.at(1);
                syncTree->bcsv_2() = event.csv_Bjets.at(1);
                syncTree->m_bb() = eventInfo.Hbb.M();
                syncTree->m_ttbb() = eventInfo.resonance.M();
            }

            if (event.nBjets >= 3 && event.csv_Bjets.at(0) > cuts::Htautau_Summer13::btag::CSVM &&
                    event.csv_Bjets.at(1) > cuts::Htautau_Summer13::btag::CSVM &&
                    event.csv_Bjets.at(2) > cuts::Htautau_Summer13::btag::CSVM){
                syncTree->bpt_3() = event.pt_Bjets.at(2);
                syncTree->beta_3() = event.eta_Bjets.at(2);
                syncTree->bphi_3() = event.phi_Bjets.at(2);
                syncTree->bcsv_3() = event.csv_Bjets.at(2);
            }

            syncTree->Fill();
        }
    }

private:
    std::shared_ptr<ntuple::SyncTree> syncTree;
    EventId_ToES_Map eventId_ToES_Map;
    unsigned inclusive, passed, passedEventType;
};
