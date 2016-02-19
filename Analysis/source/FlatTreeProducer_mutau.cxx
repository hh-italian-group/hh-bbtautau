/*! Generate flat-tree for Htautau analysis using looser selection.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseFlatTreeProducer.h"

namespace analysis {
struct SelectionResults_mutau : public SelectionResults {
    finalState::bbMuTaujet muTau_MC;
    CandidatePtr GetMuon() const { return higgs->GetDaughter(Candidate::Type::Muon); }
    CandidatePtr GetTau() const { return higgs->GetDaughter(Candidate::Type::Tau); }

    virtual CandidatePtr GetLeg(size_t leg_id) const override
    {
        if(leg_id == 1) return GetMuon();
        if(leg_id == 2) return GetTau();
        throw exception("Bad leg id = ") << leg_id;
    }

    virtual const finalState::bbTauTau& GetFinalStateMC() const override { return muTau_MC; }
};

class EventWeights_mutau : public EventWeights {
public:
    EventWeights_mutau(bool is_data, bool is_embedded, bool apply_pu_weight, bool _apply_DM_weight,
                       const std::string& pu_reweight_file_name,
                       double _max_available_pu, double _default_pu_weight, bool _applyJetToTauFakeRate)
        : EventWeights(is_data, is_embedded, apply_pu_weight, _apply_DM_weight, pu_reweight_file_name, _max_available_pu, _default_pu_weight),
          applyJetToTauFakeRate(_applyJetToTauFakeRate) {}

protected:
    virtual double CalculateTriggerWeight(CandidatePtr leg) override
    {
        using namespace analysis::Htautau_Summer13::trigger;
        using namespace analysis::Htautau_Summer13::trigger::Run2012ABCD::MuTau;
        typedef std::pair<bool, Candidate::Type> Key;

        static const std::map<Key, TriggerFunction> trigger_functions = {
            { { false, Candidate::Type::Muon }, &CalculateMuonWeight },
            { { false, Candidate::Type::Tau }, &CalculateTauWeight },
            { { true, Candidate::Type::Muon }, &CalculateMuonTurnOnCurveData },
            { { true, Candidate::Type::Tau }, &CalculateTauTurnOnCurveData }
        };

        const Key key(IsEmbedded(), leg->GetType());
        if(!trigger_functions.count(key))
            throw exception("Bad leg type ") << leg->GetType();
        return trigger_functions.at(key)(leg->GetMomentum());
    }

    virtual double CalculateIsoWeight(CandidatePtr leg) override
    {
        using namespace cuts::Htautau_Summer13::MuTau::muonISOscaleFactor;

        if(leg->GetType() == Candidate::Type::Tau)
            return 1;
        if(leg->GetType() != Candidate::Type::Muon)
            throw exception("Bad leg type ") << leg->GetType();

        const double mu_pt = leg->GetMomentum().Pt(), mu_eta = std::abs(leg->GetMomentum().Eta());
        if(mu_pt < pt.at(0))
            throw exception("No information about ISO. Muon pt is too small");
        const size_t pt_bin = mu_pt < pt.at(1) ? 0 : 1;
        if(mu_eta >= eta.at(2))
            throw exception("No information about ISO. Muon eta is too big");
        const size_t eta_bin = mu_eta < eta.at(0) ? 0 : ( mu_eta < eta.at(1) ? 1 : 2 );
        return scaleFactors.at(pt_bin).at(eta_bin);
    }

    virtual double CalculateIdWeight(CandidatePtr leg) override
    {
        using namespace cuts::Htautau_Summer13::MuTau::muonIDscaleFactor;

        if(leg->GetType() == Candidate::Type::Tau)
            return 1;
        if(leg->GetType() != Candidate::Type::Muon)
            throw exception("Bad leg type ") << leg->GetType();

        const double mu_pt = leg->GetMomentum().Pt(), mu_eta = std::abs(leg->GetMomentum().Eta());
        if(mu_pt < pt.at(0))
            throw std::runtime_error("No information about ID. Muon pt is too small");
        const size_t pt_bin = mu_pt < pt.at(1) ? 0 : 1;
        if(mu_eta >= eta.at(2))
            throw std::runtime_error("No information about ID. Muon eta is too big");
        const size_t eta_bin = mu_eta < eta.at(0) ? 0 : ( mu_eta < eta.at(1) ? 1 : 2 );
        return scaleFactors.at(pt_bin).at(eta_bin);
    }

    virtual double CalculateFakeWeight(CandidatePtr leg) override
    {
        using namespace cuts::Htautau_Summer13::jetToTauFakeRateWeight;

        if(leg->GetType() == Candidate::Type::Muon)
            return 1;
        if(leg->GetType() != Candidate::Type::Tau)
            throw exception("Bad leg type ") << leg->GetType();

        return applyJetToTauFakeRate ? CalculateJetToTauFakeWeight(leg->GetMomentum().Pt()) : 1;
    }

private:
    bool applyJetToTauFakeRate;
};

} // namespace analysis

class FlatTreeProducer_mutau : public virtual analysis::BaseFlatTreeProducer {
public:
    FlatTreeProducer_mutau(const std::string& inputFileName, const std::string& outputFileName,
                           const std::string& configFileName, const std::string& _prefix = "none",
                           size_t _maxNumberOfEvents = 0,
                           std::shared_ptr<ntuple::FlatTree> _flatTree = std::shared_ptr<ntuple::FlatTree>())
        : BaseFlatTreeProducer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents, _flatTree),
          baseAnaData(outputFile),
          eventWeights(!config.isMC(), config.IsEmbeddedSample(), config.ApplyPUreweight(), config.ApplyDMweight(),
                       config.PUreweight_fileName(), config.PUreweight_maxAvailablePU(),
                       config.PUreweight_defaultWeight(), config.ApplyJetToTauFakeRate())
    {
        if(config.ApplyRecoilCorrection())
            recoilCorrectionProducer_mutau = std::shared_ptr<analysis::RecoilCorrectionProducer>(
                        new analysis::RecoilCorrectionProducer(config.RecoilCorrection_fileCorrectTo_MuTau(),
                                                               config.RecoilCorrection_fileZmmData_MuTau(),
                                                               config.RecoilCorrection_fileZmmMC_MuTau()));
    }

    virtual analysis::BaseAnalyzerData& GetAnaData() override { return baseAnaData; }
    virtual analysis::EventWeights& GetEventWeights() override { return eventWeights; }
    virtual analysis::RecoilCorrectionProducer& GetRecoilCorrectionProducer() override
    {
        return *recoilCorrectionProducer_mutau;
    }

protected:
    virtual analysis::SelectionResults& ApplyBaselineSelection() override
    {
        using namespace analysis;
        using namespace cuts::Htautau_Summer13;
        using namespace cuts::Htautau_Summer13::MuTau;

        selection = SelectionResults_mutau();

        cuts::Cutter cut(&GetAnaData().Selection("event"));
        cut(true, "total");

        cut(FindAnalysisFinalState(selection.muTau_MC) || !config.RequireSpecificFinalState(), "spec_final_state");
        cut(!config.isDYEmbeddedSample() || GenFilterForZevents(selection.muTau_MC), "genFilter");

        const auto& selectedTriggerPath = config.IsEmbeddedSample()
                ? DYEmbedded::trigger::hltPaths : trigger::hltPaths;
        cut(HaveTriggerFired(selectedTriggerPath), "trigger");

        selection.vertices = CollectVertices();
        cut(selection.vertices.size(), "vertex");
        primaryVertex = selection.GetPrimaryVertex();

        const auto z_muons = CollectZmuons();
        const auto z_muon_candidates = FindCompatibleObjects(z_muons, ZmumuVeto::deltaR, Candidate::Type::Z,
                                                             "Z_mu_mu", 0);
        cut(!z_muon_candidates.size(), "z_mumu_veto");

        const auto electrons_bkg = CollectBackgroundElectrons();
        cut(!electrons_bkg.size(), "no_electron");

        const auto signal_muons = CollectSignalMuons(); //rename to signalMuons

        const auto muons_extra = CollectBackgroundMuons();
        const bool have_bkg_muon = muons_extra.size() > 1 || signal_muons.size() > 1 ||
                ( muons_extra.size() == 1 && signal_muons.size() == 1 && *muons_extra.front() != *signal_muons.front() );
        cut(!have_bkg_muon, "no_bkg_muon");

        const auto allmuons = CollectMuons();
        cut(allmuons.size(), "muon_cand");

        correctedTaus = config.ApplyTauESCorrection()
                ? ApplyTauCorrections(selection.muTau_MC.hadronic_taus,false) : event->taus();

        const auto alltaus = CollectTaus();
        cut(alltaus.size(), "tau_cand");

        const auto signaltaus = CollectSignalTaus() ;

        // First check OS, isolated higgs candidates
        // If no OS candidate, keep any higgs-ish candidates for bkg estimation (don't cut on sign nor isolation)
        auto higgses = FindCompatibleObjects(signal_muons, signaltaus, DeltaR_betweenSignalObjects, Candidate::Type::Higgs,
                                             "H_mu_tau", 0);
        //check SS Isolated higgs candidates
        if(!higgses.size())
            higgses = FindCompatibleObjects(signal_muons, signaltaus, DeltaR_betweenSignalObjects,
                                            Candidate::Type::Higgs, "H_mu_tau");

        //check OS notIsolated higgs candidates
        if(!higgses.size())
            higgses = FindCompatibleObjects(allmuons, signaltaus, DeltaR_betweenSignalObjects,
                                            Candidate::Type::Higgs, "H_mu_tau", 0);

        //check SS notIsolated higgs candidates
        if(!higgses.size())
            higgses = FindCompatibleObjects(allmuons, signaltaus, DeltaR_betweenSignalObjects,
                                            Candidate::Type::Higgs, "H_mu_tau");


        cut(higgses.size(), "mu_tau");

        const auto higgsTriggered = config.IsEmbeddedSample() ? higgses :
                                                                ApplyTriggerMatch(higgses,trigger::hltPaths,false);

        cut(higgsTriggered.size(), "trigger obj match");

        selection.higgs = SelectSemiLeptonicHiggs(higgsTriggered);
        GetAnaData().Htautau_Mass().Fill(selection.higgs->GetMomentum().M());
        selection.eventType = DoEventCategorization(selection.higgs);

        cut(!config.isDYEmbeddedSample() || selection.eventType == ntuple::EventType::ZTT, "tau match with MC truth");

        if (!config.isMC() || config.isDYEmbeddedSample()){
            selection.pfMET = mvaMetProducer.ComputePFMet(event->pfCandidates(), primaryVertex);
        }
        else
            selection.pfMET = event->metPF();

        const ntuple::MET mvaMet = ComputeMvaMet(selection.higgs, selection.vertices);

        const ntuple::MET correctedMET = config.ApplyTauESCorrection()
                ? ApplyTauCorrectionsToMVAMET(mvaMet, correctedTaus) : mvaMet;

//        std::cout << "correctedMET: " << correctedMET.pt << ", phi:" << correctedMET.phi << ", mt: " <<
//                      analysis::Calculate_MT(selection.higgs->GetDaughter(analysis::Candidate::Type::Muon)->GetMomentum(),correctedMET.pt,correctedMET.phi) << std::endl;
        const auto looseJets = CollectLooseJets();
        selection.jetsPt20 = FilterCompatibleObjects(looseJets, selection.higgs, jetID::deltaR_signalObjects);


        selection.jets = CollectJets(selection.jetsPt20);
        selection.bjets_all = CollectBJets(selection.jetsPt20, false, false);
        selection.retagged_bjets = CollectBJets(selection.jetsPt20, config.isMC(), true);


        selection.MET_with_recoil_corrections = ApplyRecoilCorrections(selection.higgs, selection.muTau_MC.Higgs_TauTau,
                                                                       selection.jets.size(), correctedMET);
//        std::cout << "recoilMET: " << selection.MET_with_recoil_corrections.pt <<
//                     ", phi:" << selection.MET_with_recoil_corrections.phi << ", mt: " <<
//                     analysis::Calculate_MT(selection.higgs->GetDaughter(analysis::Candidate::Type::Muon)->GetMomentum(),selection.MET_with_recoil_corrections.pt,selection.MET_with_recoil_corrections.phi) << std::endl;
        return selection;
    }

    virtual void SelectMuon(const analysis::CandidatePtr& muon, analysis::SelectionManager& selectionManager,
                            cuts::Cutter& cut) override
    {
        using namespace cuts::Htautau_Summer13::MuTau;
        using namespace cuts::Htautau_Summer13::MuTau::muonID;
        const ntuple::Muon& object = muon->GetNtupleObject<ntuple::Muon>();

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(isGlobalMuonPromptTight) == isGlobalMuonPromptTight, "tight");
        cut(X(isPFMuon) == isPFMuon, "PF");
        cut(X(nMatchedStations) > nMatched_Stations, "stations");
        cut(X(pixHits) > pixHits, "pix_hits");
        cut(X(trackerLayersWithMeasurement) > trackerLayersWithMeasurement, "layers");
        const double DeltaZ = std::abs(object.vz - primaryVertex->GetPosition().Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double dB_PV = analysis::Calculate_dxy(mu_vertex, primaryVertex->GetPosition(), muon->GetMomentum());
        cut(std::abs( Y(dB_PV) ) < dB, "dB");
    }

    virtual void SelectSignalMuon(const analysis::CandidatePtr& muon,
                                  analysis::SelectionManager& selectionManager, cuts::Cutter& cut) override
    {
        using namespace cuts::Htautau_Summer13::MuTau;
        using namespace cuts::Htautau_Summer13::MuTau::muonID;
        const ntuple::Muon& object = muon->GetNtupleObject<ntuple::Muon>();

        SelectMuon(muon, selectionManager, cut);
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");
    }

    virtual void SelectTau(const analysis::CandidatePtr& tau, analysis::SelectionManager& selectionManager,
                           cuts::Cutter& cut) override
    {
        using namespace cuts::Htautau_Summer13::MuTau;
        using namespace cuts::Htautau_Summer13::MuTau::tauID;
        const ntuple::Tau& object = tau->GetNtupleObject<ntuple::Tau>();

        cut(true, ">0 tau cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(againstMuonLoose) > cuts::skim::MuTau::tauID::againstMuonLoose, "vs_mu_loose");
        cut(X(againstElectronLoose) > againstElectronLoose, "vs_e_loose");
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits)
            < cuts::skim::MuTau::tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits, "relaxed_Iso3Hits");
        const double DeltaZ = std::abs(object.vz - primaryVertex->GetPosition().Z());
        cut(Y(DeltaZ)  < dz, "dz");
    }

    virtual void SelectSignalTau(const analysis::CandidatePtr& tau,
                                 analysis::SelectionManager& selectionManager, cuts::Cutter& cut) override
    {
        using namespace cuts::Htautau_Summer13::MuTau;
        using namespace cuts::Htautau_Summer13::MuTau::tauID;
        const ntuple::Tau& object = tau->GetNtupleObject<ntuple::Tau>();

        SelectTau(tau, selectionManager, cut);
        cut(X(againstMuonTight) > againstMuonTight, "vs_mu_tight");
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits) < byCombinedIsolationDeltaBetaCorrRaw3Hits, "loose_Iso3Hits");
    }

    analysis::CandidatePtrVector CollectZmuons()
    {
        const auto base_selector = [&](const analysis::CandidatePtr& candidate,
                                       analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
            { SelectZmuon(candidate, selectionManager, cut); };
        return CollectObjects<analysis::Candidate>("z_muons", base_selector, event->muons());
    }

    virtual void SelectZmuon(const analysis::CandidatePtr& muon, analysis::SelectionManager& selectionManager,
                             cuts::Cutter& cut)
    {
        using namespace cuts::Htautau_Summer13::MuTau::ZmumuVeto;
        const ntuple::Muon& object = muon->GetNtupleObject<ntuple::Muon>();

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex->GetPosition().Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = analysis::Calculate_dxy(mu_vertex, primaryVertex->GetPosition(), muon->GetMomentum());
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(isTrackerMuon) == isTrackerMuon, "trackerMuon");
        cut(X(isPFMuon) == isPFMuon, "PFMuon");
        cut(X(pfRelIso) < pfRelIso, "pFRelIso");
    }

    bool FindAnalysisFinalState(analysis::finalState::bbMuTaujet& final_state)
    {
        const bool base_result = BaseFlatTreeProducer::FindAnalysisFinalState(final_state);
        if(!base_result)
            return base_result;

        for(const analysis::VisibleGenObject& tau_MC : final_state.taus) {
//            if(tau_MC.finalStateChargedLeptons.size() == 1
//                    && (*tau_MC.finalStateChargedLeptons.begin())->pdg.Code == particles::mu)
//                final_state.muon = *tau_MC.finalStateChargedLeptons.begin();
            analysis::GenParticlePtrVector tauProducts;
            if (analysis::FindDecayProducts(*tau_MC.origin,analysis::TauMuonicDecay,tauProducts,false))
                final_state.muon = tauProducts.at(0);
//            else if(tau_MC.finalStateChargedHadrons.size() >= 1)
            else if(!analysis:: IsLeptonicTau(*tau_MC.origin)){
                final_state.tau_jet = &tau_MC;
            }
        }

        if (!final_state.muon || !final_state.tau_jet) return false;
        return true;
    }

    virtual void FillFlatTree(const analysis::SelectionResults& /*selection*/) override
    {
        static const float default_value = ntuple::DefaultFloatFillValueForFlatTree();
        static const float default_int_value = ntuple::DefaultIntegerFillValueForFlatTree();

        BaseFlatTreeProducer::FillFlatTree(selection);

        const analysis::CandidatePtr& muon = selection.GetMuon();
        const ntuple::Muon& ntuple_muon = muon->GetNtupleObject<ntuple::Muon>();
        const analysis::CandidatePtr& tau = selection.GetTau();
        const ntuple::Tau& ntuple_tau = tau->GetNtupleObject<ntuple::Tau>();

        flatTree->channel() = static_cast<int>(analysis::Channel::MuTau);
        flatTree->pfRelIso_1() = ntuple_muon.pfRelIso;
        flatTree->mva_1() = 0;
        flatTree->passid_1() = true;
        flatTree->passiso_1() = ntuple_muon.pfRelIso < cuts::Htautau_Summer13::MuTau::muonID::pFRelIso;
        flatTree->decayMode_1() = default_int_value;
        flatTree->iso_1() = ntuple_muon.pfRelIso;
        flatTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1() = default_value;
        flatTree->againstElectronLooseMVA_1() = false;
        flatTree->againstElectronMediumMVA_1() = false;
        flatTree->againstElectronTightMVA_1() = false;
        flatTree->againstElectronVTightMVA_1() = false;
        flatTree->againstElectronLooseMVA_custom_1() = false;
        flatTree->againstElectronMediumMVA_custom_1() = false;
        flatTree->againstElectronTightMVA_custom_1() = false;
        flatTree->againstElectronVTightMVA_custom_1() = false;
        flatTree->againstElectronLoose_1() = false;
        flatTree->againstElectronMedium_1() = false;
        flatTree->againstElectronTight_1() = false;
        flatTree->againstMuonLoose_1() = false;
        flatTree->againstMuonMedium_1() = false;
        flatTree->againstMuonTight_1() = false;
        flatTree->againstElectronMVA3raw_1() = default_int_value;
        flatTree->byIsolationMVA2raw_1() = default_int_value;
        flatTree->againstElectronLooseMVA_custom_2() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau, 0, true);
        flatTree->againstElectronMediumMVA_custom_2() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau, 1, true);
        flatTree->againstElectronTightMVA_custom_2() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau, 2, true);
        flatTree->againstElectronVTightMVA_custom_2() = cuts::Htautau_Summer13::customTauMVA::ComputeAntiElectronMVA3New(
                    ntuple_tau, 3, true);


        const bool muon_matched = analysis::HasMatchWithMCParticle(muon->GetMomentum(), selection.muTau_MC.muon,
                                                                   cuts::DeltaR_MC_Match);
        if(muon_matched) {
            const TLorentzVector& momentum = selection.muTau_MC.muon->momentum;
            flatTree->pt_1_MC   () = momentum.Pt()  ;
            flatTree->phi_1_MC  () = momentum.Phi() ;
            flatTree->eta_1_MC  () = momentum.Eta() ;
            flatTree->m_1_MC    () = momentum.M()   ;
            flatTree->pt_1_visible_MC   () = momentum.Pt()  ;
            flatTree->phi_1_visible_MC  () = momentum.Phi() ;
            flatTree->eta_1_visible_MC  () = momentum.Eta() ;
            flatTree->m_1_visible_MC    () = momentum.M()   ;
            flatTree->pdgId_1_MC() = selection.muTau_MC.muon->pdg.ToInteger();
        } else {
            flatTree->pt_1_MC   () = default_value ;
            flatTree->phi_1_MC  () = default_value ;
            flatTree->eta_1_MC  () = default_value ;
            flatTree->m_1_MC    () = default_value ;
            flatTree->pt_1_visible_MC   () = default_value ;
            flatTree->phi_1_visible_MC  () = default_value ;
            flatTree->eta_1_visible_MC  () = default_value ;
            flatTree->m_1_visible_MC    () = default_value ;
            flatTree->pdgId_1_MC() = particles::NONEXISTENT.RawCode();
        }

        const bool tau_matched = analysis::HasMatchWithMCObject(tau->GetMomentum(), selection.muTau_MC.tau_jet,
                                                                cuts::DeltaR_MC_Match);
        if(tau_matched) {
            const TLorentzVector& momentum = selection.muTau_MC.tau_jet->origin->momentum;
            flatTree->pt_2_MC   () = momentum.Pt()  ;
            flatTree->phi_2_MC  () = momentum.Phi() ;
            flatTree->eta_2_MC  () = momentum.Eta() ;
            flatTree->m_2_MC    () = momentum.M()   ;
            const TLorentzVector& visible_momentum = selection.muTau_MC.tau_jet->visibleMomentum;
            flatTree->pt_2_visible_MC   () = visible_momentum.Pt()  ;
            flatTree->phi_2_visible_MC  () = visible_momentum.Phi() ;
            flatTree->eta_2_visible_MC  () = visible_momentum.Eta() ;
            flatTree->m_2_visible_MC    () = visible_momentum.M()   ;
            flatTree->pdgId_2_MC() = selection.muTau_MC.tau_jet->origin->pdg.ToInteger();
        } else {
            flatTree->pt_2_MC   () = default_value ;
            flatTree->phi_2_MC  () = default_value ;
            flatTree->eta_2_MC  () = default_value ;
            flatTree->m_2_MC    () = default_value ;
            flatTree->pt_2_visible_MC   () = default_value ;
            flatTree->phi_2_visible_MC  () = default_value ;
            flatTree->eta_2_visible_MC  () = default_value ;
            flatTree->m_2_visible_MC    () = default_value ;
            flatTree->pdgId_2_MC() = particles::NONEXISTENT.RawCode();
        }

        flatTree->Fill();
    }

protected:
    analysis::BaseAnalyzerData baseAnaData;
    analysis::SelectionResults_mutau selection;
    std::shared_ptr<analysis::RecoilCorrectionProducer> recoilCorrectionProducer_mutau;
    analysis::EventWeights_mutau eventWeights;
};

#include "METPUSubtraction/interface/GBRProjectDict.cxx"
