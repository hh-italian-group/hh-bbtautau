/*! Definition of BaseFlatTreeProducer class, the base class for flat tree producers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "h-tautau/Analysis/include/FlatTree.h"

#include "BaseAnalyzer.h"

namespace analysis {

class BaseFlatTreeProducer : public BaseAnalyzer {
public:
    BaseFlatTreeProducer(const std::string& inputFileName, const std::string& outputFileName,
                         const std::string& configFileName, const std::string& _prefix = "none",
                         size_t _maxNumberOfEvents = 0,
                         std::shared_ptr<ntuple::FlatTree> _flatTree = std::shared_ptr<ntuple::FlatTree>())
        : BaseAnalyzer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents),
          flatTree(_flatTree), writeFlatTree(!flatTree)
    {
        if(!flatTree)
            flatTree = std::shared_ptr<ntuple::FlatTree>(new ntuple::FlatTree("flatTree", outputFile.get(), false));
    }

    virtual ~BaseFlatTreeProducer() override
    {
        if(writeFlatTree)
            flatTree->Write();
    }

    virtual void ProcessEvent() override
    {
        using namespace analysis;

        SelectionResults& selection = ApplyBaselineSelection();
        selection.svfitResults = sv_fit::CombinedFit({ selection.GetLeg(1), selection.GetLeg(2) },
                                                     selection.MET_with_recoil_corrections, true, true);
        selection.kinfitResults = RunKinematicFit(selection.bjets_all, *selection.higgs,
                                                  selection.MET_with_recoil_corrections);

        if(config.isMC()){
            if(config.ApplyDMweight())
                GetEventWeights().SetGenTaus(selection.GetFinalStateMC());
            GetEventWeights().CalculateSelectionDependentWeights(selection);
        }

        FillFlatTree(selection);
    }

protected:
    virtual void FillFlatTree(const SelectionResults& selection)
    {
        static const float default_value = ntuple::DefaultFloatFillValueForFlatTree();

        // Event
        flatTree->run() = event->eventInfo().run;
        flatTree->lumi() = event->eventInfo().lumis;
        flatTree->evt() = event->eventInfo().EventId;
        flatTree->eventType() = static_cast<int>(selection.eventType);
        flatTree->eventEnergyScale() = static_cast<int>(eventEnergyScale);

        flatTree->npv() = selection.vertices.size();
        if (config.ApplyPUreweight()){
            const size_t bxIndex = tools::find_index(event->eventInfo().bunchCrossing, 0);
            if(bxIndex >= event->eventInfo().bunchCrossing.size())
                throw std::runtime_error("in-time BX not found");
            flatTree->npu() = event->eventInfo().trueNInt.at(bxIndex);
        }

        // Weights
        flatTree->puweight() = GetEventWeights().GetPileUpWeight();
        flatTree->trigweight_1() = GetEventWeights().GetTriggerWeight(1);
        flatTree->trigweight_2() = GetEventWeights().GetTriggerWeight(2);
        flatTree->idweight_1() = GetEventWeights().GetIdWeight(1);
        flatTree->idweight_2() = GetEventWeights().GetIdWeight(2);
        flatTree->isoweight_1() = GetEventWeights().GetIsoWeight(1);
        flatTree->isoweight_2() = GetEventWeights().GetIsoWeight(2);
        flatTree->fakeweight_1() = GetEventWeights().GetFakeWeight(1); //
        flatTree->fakeweight_2() = GetEventWeights().GetFakeWeight(2); // e -> tau fake rate * jet -> tau fake rate
        flatTree->weight() = GetEventWeights().GetFullWeight();
        flatTree->embeddedWeight() = GetEventWeights().GetEmbeddedWeight();
        flatTree->decayModeWeight_1() = GetEventWeights().GetDecayModeWeight(1);
        flatTree->decayModeWeight_2() = GetEventWeights().GetDecayModeWeight(2);

        // HTT candidate
        flatTree->mvis() = selection.higgs->GetMomentum().M();
        flatTree->m_sv_vegas() = selection.svfitResults.fit_vegas.has_valid_mass
                ? selection.svfitResults.fit_vegas.mass : default_value;
        flatTree->m_sv_MC() = selection.svfitResults.fit_mc.has_valid_mass
                ? selection.svfitResults.fit_mc.mass : default_value;
        flatTree->pt_sv_MC() = selection.svfitResults.fit_mc.has_valid_momentum
                ? selection.svfitResults.fit_mc.momentum.Pt() : default_value;
        flatTree->eta_sv_MC() = selection.svfitResults.fit_mc.has_valid_momentum
                ? selection.svfitResults.fit_mc.momentum.Eta() : default_value;
        flatTree->phi_sv_MC() = selection.svfitResults.fit_mc.has_valid_momentum
                ? selection.svfitResults.fit_mc.momentum.Phi() : default_value;

        flatTree->DeltaR_leptons() = selection.GetLeg(1)->GetMomentum().DeltaR(selection.GetLeg(2)->GetMomentum()) ;
        flatTree->pt_tt()          = (selection.GetLeg(1)->GetMomentum() + selection.GetLeg(2)->GetMomentum()).Pt();

        // Kinematic fit
        flatTree->kinfit_bb_tt_mass() = selection.kinfitResults.mass;
        flatTree->kinfit_bb_tt_convergence() = selection.kinfitResults.convergence;
        flatTree->kinfit_bb_tt_chi2() = selection.kinfitResults.chi2;
        flatTree->kinfit_bb_tt_pull_balance() = selection.kinfitResults.pull_balance;

        // Hhh generator info candidate
        if(selection.GetFinalStateMC().resonance) {
            const TLorentzVector& momentum = selection.GetFinalStateMC().resonance->momentum;
            flatTree->pt_resonance_MC()    = momentum.Pt()  ;
            flatTree->eta_resonance_MC()   = momentum.Eta() ;
            flatTree->phi_resonance_MC()   = momentum.Phi() ;
            flatTree->mass_resonance_MC()  = momentum.M()   ;
            flatTree->pdgId_resonance_MC() = selection.GetFinalStateMC().resonance->pdg.ToInteger();
        } else {
            flatTree->pt_resonance_MC()    = default_value  ;
            flatTree->eta_resonance_MC()   = default_value  ;
            flatTree->phi_resonance_MC()   = default_value  ;
            flatTree->mass_resonance_MC()  = default_value  ;
            flatTree->pdgId_resonance_MC() = particles::NONEXISTENT.RawCode() ;
        }

        if(selection.GetFinalStateMC().Higgs_TauTau) {
            const TLorentzVector& momentum = selection.GetFinalStateMC().Higgs_TauTau->momentum;
            flatTree->pt_Htt_MC()    = momentum.Pt()  ;
            flatTree->eta_Htt_MC()   = momentum.Eta() ;
            flatTree->phi_Htt_MC()   = momentum.Phi() ;
            flatTree->mass_Htt_MC()  = momentum.M()   ;
            flatTree->pdgId_Htt_MC() = selection.GetFinalStateMC().Higgs_TauTau->pdg.ToInteger();
        } else {
            flatTree->pt_Htt_MC()    = default_value  ;
            flatTree->eta_Htt_MC()   = default_value  ;
            flatTree->phi_Htt_MC()   = default_value  ;
            flatTree->mass_Htt_MC()  = default_value  ;
            flatTree->pdgId_Htt_MC() = particles::NONEXISTENT.RawCode() ;
        }

        if(selection.GetFinalStateMC().Higgs_BB) {
            const TLorentzVector& momentum = selection.GetFinalStateMC().Higgs_BB->momentum;
            flatTree->pt_Hbb_MC()    = momentum.Pt()  ;
            flatTree->eta_Hbb_MC()   = momentum.Eta() ;
            flatTree->phi_Hbb_MC()   = momentum.Phi() ;
            flatTree->mass_Hbb_MC()  = momentum.M()   ;
            flatTree->pdgId_Hbb_MC() = selection.GetFinalStateMC().Higgs_BB->pdg.ToInteger();
        } else {
            flatTree->pt_Hbb_MC()    = default_value  ;
            flatTree->eta_Hbb_MC()   = default_value  ;
            flatTree->phi_Hbb_MC()   = default_value  ;
            flatTree->mass_Hbb_MC()  = default_value  ;
            flatTree->pdgId_Hbb_MC() = particles::NONEXISTENT.RawCode()  ;
        }

        // needs to be filles with NUP!
        // https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauTauAnalyzer.py#L51
        if (config.MaxTreeVersion() == 2)
            flatTree->n_extraJets_MC() = event->genEvent().nup;
        else
            flatTree->n_extraJets_MC() = default_value;

        // MET
        const TLorentzVector MET_momentum = MakeLorentzVectorPtEtaPhiM(selection.MET_with_recoil_corrections.pt, 0,
                                                                       selection.MET_with_recoil_corrections.phi, 0);
        flatTree->pt_tt_MET() = (selection.GetLeg(1)->GetMomentum() + selection.GetLeg(2)->GetMomentum()
                                 + MET_momentum).Pt();

        flatTree->met() = selection.pfMET.pt;
        flatTree->metphi() = selection.pfMET.phi;
        flatTree->mvamet() = MET_momentum.Pt();
        flatTree->mvametphi() = MET_momentum.Phi();
        //flatTree->pzetavis();
        //flatTree->pzetamiss();
        if(selection.pfMET.significanceMatrix.size()) {
            const TMatrixD metPFcov = ntuple::VectorToSignificanceMatrix(selection.pfMET.significanceMatrix);
            flatTree->metcov00() = metPFcov[0][0];
            flatTree->metcov01() = metPFcov[0][1];
            flatTree->metcov10() = metPFcov[1][0];
            flatTree->metcov11() = metPFcov[1][1];
        }
        const TMatrixD metMVAcov =
                ntuple::VectorToSignificanceMatrix(selection.MET_with_recoil_corrections.significanceMatrix);
        flatTree->mvacov00() = metMVAcov[0][0];
        flatTree->mvacov01() = metMVAcov[0][1];
        flatTree->mvacov10() = metMVAcov[1][0];
        flatTree->mvacov11() = metMVAcov[1][1];

        // Leg 1, lepton
        flatTree->pt_1()     = selection.GetLeg(1)->GetMomentum().Pt();
        flatTree->phi_1()    = selection.GetLeg(1)->GetMomentum().Phi();
        flatTree->eta_1()    = selection.GetLeg(1)->GetMomentum().Eta();
        flatTree->m_1()      = selection.GetLeg(1)->GetMomentum().M();
        flatTree->energy_1() = selection.GetLeg(1)->GetMomentum().E();
        flatTree->q_1()      = selection.GetLeg(1)->GetCharge();
        flatTree->mt_1()     = Calculate_MT(selection.GetLeg(1)->GetMomentum(), MET_momentum.Pt(), MET_momentum.Phi());
        flatTree->d0_1()     = Calculate_dxy(selection.GetLeg(1)->GetVertexPosition(), primaryVertex->GetPosition(),
                                             selection.GetLeg(1)->GetMomentum());
        flatTree->dZ_1()     = selection.GetLeg(1)->GetVertexPosition().Z() - primaryVertex->GetPosition().Z();

        // Leg 2, tau
        flatTree->pt_2()     = selection.GetLeg(2)->GetMomentum().Pt();
        flatTree->phi_2()    = selection.GetLeg(2)->GetMomentum().Phi();
        flatTree->eta_2()    = selection.GetLeg(2)->GetMomentum().Eta();
        flatTree->m_2()      = selection.GetLeg(2)->GetMomentum().M();
        flatTree->energy_2() = selection.GetLeg(2)->GetMomentum().E();
        flatTree->q_2()      = selection.GetLeg(2)->GetCharge();
        flatTree->mt_2()     = Calculate_MT(selection.GetLeg(2)->GetMomentum(), MET_momentum.Pt(), MET_momentum.Phi());
        flatTree->d0_2()     = Calculate_dxy(selection.GetLeg(2)->GetVertexPosition(), primaryVertex->GetPosition(),
                                             selection.GetLeg(2)->GetMomentum());
        flatTree->dZ_2()     = selection.GetLeg(2)->GetVertexPosition().Z() - primaryVertex->GetPosition().Z();

        // RM: for the three channels, mt, et, tt this leg is always a tau
        const ntuple::Tau& ntuple_tau_leg2 = selection.GetLeg(2)->GetNtupleObject<ntuple::Tau>();
        flatTree->decayMode_2()                                = ntuple_tau_leg2.decayMode;
        flatTree->iso_2()                  = ntuple_tau_leg2.byIsolationMVAraw;
        flatTree->againstElectronLooseMVA_2() = ntuple_tau_leg2.againstElectronLooseMVA3;
        flatTree->againstElectronMediumMVA_2() = ntuple_tau_leg2.againstElectronMediumMVA3;
        flatTree->againstElectronTightMVA_2() = ntuple_tau_leg2.againstElectronTightMVA3;
        flatTree->againstElectronVTightMVA_2() = ntuple_tau_leg2.againstElectronVTightMVA3;

        flatTree->againstElectronLoose_2()                     = ntuple_tau_leg2.againstElectronLoose  ;
        flatTree->againstElectronMedium_2()                    = ntuple_tau_leg2.againstElectronMedium ;
        flatTree->againstElectronTight_2()                     = ntuple_tau_leg2.againstElectronTight  ;
        flatTree->againstMuonLoose_2()                         = ntuple_tau_leg2.againstMuonLoose      ;
        flatTree->againstMuonMedium_2()                        = ntuple_tau_leg2.againstMuonMedium     ;
        flatTree->againstMuonTight_2()                         = ntuple_tau_leg2.againstMuonTight      ;
        flatTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = ntuple_tau_leg2.byCombinedIsolationDeltaBetaCorrRaw3Hits ;
        flatTree->againstElectronMVA3raw_2()                   = ntuple_tau_leg2.againstElectronMVA3raw;
        flatTree->byIsolationMVA2raw_2()                       = ntuple_tau_leg2.byIsolationMVA2raw;

        // Jets
        flatTree->njets()     = selection.jets.size();
        flatTree->njetspt20() = selection.jetsPt20.size();
        flatTree->nBjets()    = selection.bjets_all.size();
        flatTree->nBjets_retagged()    = selection.retagged_bjets.size();


        for (const CandidatePtr& jet : selection.jets){
            const ntuple::Jet& ntuple_jet = jet->GetNtupleObject<ntuple::Jet>();
            flatTree->pt_jets().push_back(jet->GetMomentum().Pt());
            flatTree->eta_jets().push_back(jet->GetMomentum().Eta());
            flatTree->phi_jets().push_back(jet->GetMomentum().Phi());
            flatTree->ptraw_jets().push_back(ntuple_jet.pt_raw);
            flatTree->ptunc_jets().push_back(ntuple_jet.pt_raw);//to put uncertainties -> put pt_raw only to compile
            flatTree->mva_jets().push_back(ntuple_jet.puIdMVA);
            flatTree->passPU_jets().push_back(ntuple::JetID_MVA::PassLooseId(ntuple_jet.puIdBits));
        }

        for (const CandidatePtr& jet : selection.bjets_all) {
            const ntuple::Jet& ntuple_jet = jet->GetNtupleObject<ntuple::Jet>();

            flatTree->pt_Bjets()      .push_back( jet->GetMomentum().Pt() );
            flatTree->eta_Bjets()     .push_back( jet->GetMomentum().Eta() );
            flatTree->phi_Bjets()     .push_back( jet->GetMomentum().Phi() );
            flatTree->energy_Bjets()  .push_back( jet->GetMomentum().E() );
            flatTree->chargedHadronEF_Bjets().push_back( ntuple_jet.chargedHadronEnergyFraction );
            flatTree->neutralHadronEF_Bjets()  .push_back( ntuple_jet.neutralHadronEnergyFraction );
            flatTree->photonEF_Bjets()         .push_back( ntuple_jet.photonEnergyFraction );
            flatTree->muonEF_Bjets()  .push_back( ntuple_jet.muonEnergyFraction );
            flatTree->electronEF_Bjets()  .push_back( ntuple_jet.electronEnergyFraction );
            flatTree->csv_Bjets()     .push_back( ntuple_jet.combinedSecondaryVertexBJetTags );
            // inspect the flavour of the gen jet
            const VisibleGenObjectVector matched_bjets_MC = FindMatchedObjects(jet->GetMomentum(),
                                                                               selection.GetFinalStateMC().b_jets,
                                                                               cuts::DeltaR_MC_Match);
            const bool isJet_MC_Bjet = matched_bjets_MC.size() != 0;
            const bool isJet_MC_Bjet_withLeptonicDecay = isJet_MC_Bjet
                    && matched_bjets_MC.at(0).finalStateChargedLeptons.size() != 0;
            flatTree->isBjet_MC_Bjet()                  .push_back( isJet_MC_Bjet );
            flatTree->isBjet_MC_Bjet_withLeptonicDecay().push_back( isJet_MC_Bjet_withLeptonicDecay );
        }

        flatTree->x_PV() = primaryVertex->GetPosition().x();
        flatTree->y_PV() = primaryVertex->GetPosition().y();
        flatTree->z_PV() = primaryVertex->GetPosition().z();
    }

protected:
    std::shared_ptr<ntuple::FlatTree> flatTree;
    bool writeFlatTree;
};
} // analysis
