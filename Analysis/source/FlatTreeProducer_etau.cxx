/*! Generate flat-tree for H->tautau->e_taujet analysis using looser selection.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseFlatTreeProducer.h"

namespace analysis {
struct SelectionResults_etau : public SelectionResults {
    finalState::bbETaujet eTau_MC;
    CandidatePtr GetElectron() const { return higgs->GetDaughter(Candidate::Type::Electron); }
    CandidatePtr GetTau() const { return higgs->GetDaughter(Candidate::Type::Tau); }

    virtual CandidatePtr GetLeg(size_t leg_id) const override
    {
        if(leg_id == 1) return GetElectron();
        if(leg_id == 2) return GetTau();
        throw exception("Bad leg id = ") << leg_id;
    }

    virtual const finalState::bbTauTau& GetFinalStateMC() const override { return eTau_MC; }
};

class EventWeights_etau : public EventWeights {
public:
    EventWeights_etau(bool is_data, bool is_embedded, bool apply_pu_weight, bool _apply_DM_weight,
                      const std::string& pu_reweight_file_name,
                      double _max_available_pu, double _default_pu_weight, bool _applyJetToTauFakeRate,
                      bool _applyEtoTauFakeRate)
        : EventWeights(is_data, is_embedded, apply_pu_weight, _apply_DM_weight, pu_reweight_file_name,
                       _max_available_pu, _default_pu_weight),
          applyJetToTauFakeRate(_applyJetToTauFakeRate), applyEtoTauFakeRate(_applyEtoTauFakeRate) {}

    virtual void Reset() override
    {
        EventWeights::Reset();
        genElectrons.clear();
        has_gen_electrons = false;
    }

    void SetGenElectrons(const GenEvent& genEvent)
    {
        using namespace cuts::Htautau_Summer13::DrellYannCategorization;
        static const particles::ParticleCodes electron_code = { particles::e };

        if(has_gen_electrons)
            throw exception("Gen electrons are already set.");

        genElectrons = genEvent.GetParticles(electron_code, minimal_genParticle_pt);

        has_gen_electrons = true;
    }

protected:
    virtual double CalculateTriggerWeight(CandidatePtr leg) override
    {
        using namespace analysis::Htautau_Summer13::trigger;
        using namespace analysis::Htautau_Summer13::trigger::Run2012ABCD::ETau;
        typedef std::pair<bool, Candidate::Type> Key;

        static const std::map<Key, TriggerFunction> trigger_functions = {
            { { false, Candidate::Type::Electron }, &CalculateElectronWeight },
            { { false, Candidate::Type::Tau }, &CalculateTauWeight },
            { { true, Candidate::Type::Electron }, &CalculateElectronTurnOnCurveData },
            { { true, Candidate::Type::Tau }, &CalculateTauTurnOnCurveData }
        };

        const Key key(IsEmbedded(), leg->GetType());
        if(!trigger_functions.count(key))
            throw exception("Bad leg type ") << leg->GetType();
        return trigger_functions.at(key)(leg->GetMomentum());
    }

    virtual double CalculateIsoWeight(CandidatePtr leg) override
    {
        using namespace cuts::Htautau_Summer13::ETau::electronISOscaleFactor;

        if(leg->GetType() == Candidate::Type::Tau)
            return 1;
        if(leg->GetType() != Candidate::Type::Electron)
            throw exception("Bad leg type ") << leg->GetType();

        const double ele_pt = leg->GetMomentum().Pt(), ele_eta = std::abs(leg->GetMomentum().Eta());
        if(ele_pt < pt.at(0))
            throw exception("No information about ISO. Electron pt is too small");
        const size_t pt_bin = ele_pt < pt.at(1) ? 0 : 1;
        if(ele_eta >= eta.at(1))
            throw exception("No information about ISO. Electron eta is too big");
        const size_t eta_bin = ele_eta < eta.at(0) ? 0 : 1;
        return scaleFactors.at(pt_bin).at(eta_bin);
    }

    virtual double CalculateIdWeight(CandidatePtr leg) override
    {
        using namespace cuts::Htautau_Summer13::ETau::electronIDscaleFactor;

        if(leg->GetType() == Candidate::Type::Tau)
            return 1;
        if(leg->GetType() != Candidate::Type::Electron)
            throw exception("Bad leg type ") << leg->GetType();

        const double ele_pt = leg->GetMomentum().Pt(), ele_eta = std::abs(leg->GetMomentum().Eta());
        if(ele_pt < pt.at(0))
            throw std::runtime_error("No information about ID. Electron pt is too small");
        const size_t pt_bin = ele_pt < pt.at(1) ? 0 : 1;
        if(ele_eta >= eta.at(1))
            throw std::runtime_error("No information about ID. Electron eta is too big");
        const size_t eta_bin = ele_eta < eta.at(0) ? 0 : 1;
        return scaleFactors.at(pt_bin).at(eta_bin);
    }

    virtual double CalculateFakeWeight(CandidatePtr leg) override
    {
        using namespace cuts::Htautau_Summer13::DrellYannCategorization;
        using namespace cuts::Htautau_Summer13::electronEtoTauFakeRateWeight;
        using namespace cuts::Htautau_Summer13::jetToTauFakeRateWeight;

        if(leg->GetType() == Candidate::Type::Electron)
            return 1;
        if(leg->GetType() != Candidate::Type::Tau)
            throw exception("Bad leg type ") << leg->GetType();

        double fakeEtoTauWeight = 1;
        if(applyEtoTauFakeRate) {
            if(!has_gen_electrons)
                throw exception("Gen electrons are not set.");
            if (analysis::FindMatchedParticles(leg->GetMomentum(), genElectrons, deltaR_matchGenParticle).size() > 0)
                fakeEtoTauWeight = CalculateEtoTauFakeWeight(leg->GetNtupleObject<ntuple::Tau>().eta,
                                                             ntuple::tau_id::ConvertToHadronicDecayMode(leg->GetNtupleObject<ntuple::Tau>().decayMode));
        }

        const double fakeJetToTauWeight = applyJetToTauFakeRate
                ? CalculateJetToTauFakeWeight(leg->GetMomentum().Pt()) : 1;

        return fakeEtoTauWeight * fakeJetToTauWeight;
    }

private:
    bool applyJetToTauFakeRate, applyEtoTauFakeRate;
    bool has_gen_electrons;
    GenParticleSet genElectrons;
};

} // namespace analysis

class FlatTreeProducer_etau : public virtual analysis::BaseFlatTreeProducer {
public:
    FlatTreeProducer_etau(const std::string& inputFileName, const std::string& outputFileName,
                          const std::string& configFileName, const std::string& _prefix = "none",
                          size_t _maxNumberOfEvents = 0,
                          std::shared_ptr<ntuple::FlatTree> _flatTree = std::shared_ptr<ntuple::FlatTree>())
        : BaseFlatTreeProducer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents, _flatTree),
          baseAnaData(outputFile),
          eventWeights(!config.isMC(), config.IsEmbeddedSample(), config.ApplyPUreweight(), config.ApplyDMweight(),
                       config.PUreweight_fileName(), config.PUreweight_maxAvailablePU(),
                       config.PUreweight_defaultWeight(), config.ApplyJetToTauFakeRate(), config.ApplyEtoTauFakeRate())
    {
        if(config.ApplyRecoilCorrection())
            recoilCorrectionProducer_etau = std::shared_ptr<analysis::RecoilCorrectionProducer>(
                        new analysis::RecoilCorrectionProducer(config.RecoilCorrection_fileCorrectTo_ETau(),
                                                               config.RecoilCorrection_fileZmmData_ETau(),
                                                               config.RecoilCorrection_fileZmmMC_ETau()));
    }

    virtual analysis::BaseAnalyzerData& GetAnaData() override { return baseAnaData; }
    virtual analysis::EventWeights& GetEventWeights() override { return eventWeights; }
    virtual analysis::RecoilCorrectionProducer& GetRecoilCorrectionProducer() override
    {
        return *recoilCorrectionProducer_etau;
    }

protected:
    virtual analysis::SelectionResults& ApplyBaselineSelection() override
    {
        using namespace analysis;
        using namespace cuts::Htautau_Summer13;
        using namespace cuts::Htautau_Summer13::ETau;

        selection = SelectionResults_etau();

        cuts::Cutter cut(&GetAnaData().Selection("event"));
        cut(true, "total");

        cut(FindAnalysisFinalState(selection.eTau_MC) || !config.RequireSpecificFinalState(), "spec_final_state");
        cut(!config.isDYEmbeddedSample() || GenFilterForZevents(selection.eTau_MC), "genFilter");

        const auto& selectedTriggerPath = config.IsEmbeddedSample()
                ? DYEmbedded::trigger::hltPaths : trigger::hltPaths;
        cut(HaveTriggerFired(selectedTriggerPath), "trigger");

        selection.vertices = CollectVertices();
        cut(selection.vertices.size(), "vertex");
        primaryVertex = selection.GetPrimaryVertex();

        const auto z_electrons = CollectZelectrons();
        const auto z_electron_candidates = FindCompatibleObjects(z_electrons, ZeeVeto::deltaR,Candidate::Type::Z,
                                                                 "Z_e_e", 0);
        cut(!z_electron_candidates.size(), "z_ee_veto");

        const auto muons_bkg = CollectBackgroundMuons();
        cut(!muons_bkg.size(), "no_muon");

        const auto signal_electrons = CollectSignalElectrons();

        const auto electrons_bkg = CollectBackgroundElectrons();
        const bool have_bkg_electron = electrons_bkg.size() > 1 || signal_electrons.size() > 1 ||
                ( electrons_bkg.size() == 1 && signal_electrons.size() == 1 && *electrons_bkg.front() != *signal_electrons.front() );
        cut(!have_bkg_electron, "no_bkg_electron");

        const auto allelectrons = CollectElectrons();
        cut(allelectrons.size(), "electron_cand");

        correctedTaus = config.ApplyTauESCorrection()
                ? ApplyTauCorrections(selection.eTau_MC.hadronic_taus,false) : event->taus();

        const auto alltaus = CollectTaus();
        cut(alltaus.size(), "tau_cand");

        const auto signaltaus = CollectSignalTaus() ;

        // First check OS, isolated higgs candidates
        // If no OS candidate, keep any higgs-ish candidates for bkg estimation (don't cut on sign nor isolation)
        auto higgses = FindCompatibleObjects(signal_electrons, signaltaus, DeltaR_betweenSignalObjects,
                                             Candidate::Type::Higgs, "H_e_tau", 0);

        //check SS isolated higgs candidates
        if (!higgses.size())
          higgses = FindCompatibleObjects(signal_electrons, signaltaus, DeltaR_betweenSignalObjects,
                                          Candidate::Type::Higgs, "H_e_tau");

        //check 0S antiisolated higgs candidates
        if (!higgses.size())
          higgses = FindCompatibleObjects(allelectrons, signaltaus, DeltaR_betweenSignalObjects,
                                          Candidate::Type::Higgs, "H_e_tau", 0);

        //check SS antiisolated higgs candidates
        if (!higgses.size())
          higgses = FindCompatibleObjects(allelectrons, signaltaus, DeltaR_betweenSignalObjects,
                                          Candidate::Type::Higgs, "H_e_tau");

        cut(higgses.size(), "e_tau");

        const auto higgsTriggered = config.IsEmbeddedSample() ? higgses :
                                                                ApplyTriggerMatch(higgses,trigger::hltPaths,false);

        cut(higgsTriggered.size(), "trigger obj match");


        selection.higgs = SelectSemiLeptonicHiggs(higgsTriggered);
        selection.eventType = DoEventCategorization(selection.higgs);

        cut(!config.isDYEmbeddedSample() || selection.eventType == ntuple::EventType::ZTT, "tau match with MC truth");

        if(config.ApplyEtoTauFakeRate())
            eventWeights.SetGenElectrons(genEvent);

        if (!config.isMC() || config.isDYEmbeddedSample()){
            selection.pfMET = mvaMetProducer.ComputePFMet(event->pfCandidates(), primaryVertex);
        }
        else
            selection.pfMET = event->metPF();

        const ntuple::MET mvaMet = ComputeMvaMet(selection.higgs, selection.vertices);

        const ntuple::MET correctedMET = config.ApplyTauESCorrection()
                ? ApplyTauCorrectionsToMVAMET(mvaMet, correctedTaus) : mvaMet;

//        std::cout << "correctedMET: " << correctedMET.pt << ", phi:" << correctedMET.phi << ", mt: " <<
//                      analysis::Calculate_MT(selection.higgs->
//                                             GetDaughter(analysis::Candidate::Type::Electron)->GetMomentum(),correctedMET.pt,correctedMET.phi) << std::endl;
        const auto looseJets = CollectLooseJets();
        selection.jetsPt20 = FilterCompatibleObjects(looseJets, selection.higgs, jetID::deltaR_signalObjects);


        selection.jets = CollectJets(selection.jetsPt20);
        selection.bjets_all = CollectBJets(selection.jetsPt20, false, false);
        selection.retagged_bjets = CollectBJets(selection.jetsPt20, config.isMC(), true);


        selection.MET_with_recoil_corrections = ApplyRecoilCorrections(selection.higgs, selection.eTau_MC.Higgs_TauTau,
                                                                       selection.jets.size(), correctedMET);
//        std::cout << "recoilMET: " << selection.MET_with_recoil_corrections.pt <<
//                     ", phi:" << selection.MET_with_recoil_corrections.phi << ", mt: " <<
//                     analysis::Calculate_MT(selection.higgs->GetDaughter(analysis::Candidate::Type::Electron)->
//                                            GetMomentum(),selection.MET_with_recoil_corrections.pt,selection.MET_with_recoil_corrections.phi) << std::endl;
        return selection;
    }

    virtual void SelectElectron(const analysis::CandidatePtr& electron, analysis::SelectionManager& selectionManager,
                                cuts::Cutter& cut) override
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::electronID;
        const ntuple::Electron& object = electron->GetNtupleObject<ntuple::Electron>();

        cut(true, ">0 ele cand");
        cut(X(pt) > pt, "pt");
        const double eta = std::abs( X(eta) );
        cut(eta < eta_high, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex->GetPosition().Z());
        cut(Y(DeltaZ)  < dz, "dz");
        cut(X(missingHits) < missingHits, "missingHits");
        cut(X(hasMatchedConversion) == hasMatchedConversion, "conversion");
        const TVector3 ele_vertex(object.vx, object.vy, object.vz);
        // same as dB
        const double d0_PV = analysis::Calculate_dxy(ele_vertex, primaryVertex->GetPosition(), electron->GetMomentum());
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig) > MVApogNonTrig[eta_index], "mva");
    }

    virtual void SelectSignalElectron(const analysis::CandidatePtr& electron,
                                      analysis::SelectionManager& selectionManager, cuts::Cutter& cut) override
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::electronID;
        const ntuple::Electron& object = electron->GetNtupleObject<ntuple::Electron>();

        SelectElectron(electron, selectionManager, cut);
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");
    }

    virtual void SelectTau(const analysis::CandidatePtr& tau, analysis::SelectionManager& selectionManager,
                           cuts::Cutter& cut) override
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::tauID;
        using namespace cuts::Htautau_Summer13::customTauMVA;
        const ntuple::Tau& object = tau->GetNtupleObject<ntuple::Tau>();

        cut(true, ">0 tau cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        cut(X(decayModeFinding) > decayModeFinding, "decay_mode");
        cut(X(againstMuonLoose) > againstMuonLoose, "vs_mu_loose");
        const bool againstElectron =  ComputeAntiElectronMVA3New(object, againstElectronMVA3_customWP_id, true);
        cut(Y(againstElectron), "vs_e_mediumMVA");
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits) <
            cuts::skim::ETau::tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits, "relaxed_Iso3Hits");
        const double DeltaZ = std::abs(object.vz - primaryVertex->GetPosition().Z());
        cut(Y(DeltaZ)  < dz, "dz");
    }

    virtual void SelectSignalTau(const analysis::CandidatePtr& tau, analysis::SelectionManager& selectionManager,
                                 cuts::Cutter& cut) override
    {
        using namespace cuts::Htautau_Summer13::ETau;
        using namespace cuts::Htautau_Summer13::ETau::tauID;
        using namespace cuts::Htautau_Summer13::customTauMVA;
        const ntuple::Tau& object = tau->GetNtupleObject<ntuple::Tau>();

        SelectTau(tau, selectionManager, cut);
        cut(X(byCombinedIsolationDeltaBetaCorrRaw3Hits) < byCombinedIsolationDeltaBetaCorrRaw3Hits, "looseIso3Hits");
    }

    analysis::CandidatePtrVector CollectZelectrons()
    {
        const auto base_selector = [&](const analysis::CandidatePtr& candidate,
                analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
            { SelectZelectron(candidate, selectionManager, cut); };
        return CollectObjects<analysis::Candidate>("z_electrons", base_selector, event->electrons());
    }

    virtual void SelectZelectron(const analysis::CandidatePtr& electron, analysis::SelectionManager& selectionManager,
                                 cuts::Cutter& cut)
    {
        using namespace cuts::Htautau_Summer13::ETau::ZeeVeto;
        const ntuple::Electron& object = electron->GetNtupleObject<ntuple::Electron>();

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex->GetPosition().Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 ele_vertex(object.vx, object.vy, object.vz);
        // same as dB
        const double d0_PV = analysis::Calculate_dxy(ele_vertex, primaryVertex->GetPosition(), electron->GetMomentum());
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(pfRelIso) < pfRelIso, "pFRelIso");
        const size_t eta_index = std::abs(object.eta) <= barrel_eta_high ? barrel_index : endcap_index;
        cut(X(sigmaIEtaIEta) < sigma_ieta_ieta[eta_index], "sigmaIetaIeta");
        cut(X(deltaEtaTrkSC) < delta_eta[eta_index], "deltaEtaSC");
        cut(X(deltaPhiTrkSC) < delta_phi[eta_index], "deltaPhiSC");
        cut(X(hcalOverEcal) < HoverE[eta_index], "HoverE");
    }

    bool FindAnalysisFinalState(analysis::finalState::bbETaujet& final_state)
    {
        const bool base_result = BaseFlatTreeProducer::FindAnalysisFinalState(final_state);
        if(!base_result)
            return base_result;

        for (const analysis::VisibleGenObject& tau_MC : final_state.taus) {
//            if(tau_MC.finalStateChargedLeptons.size() == 1
//                    && (*tau_MC.finalStateChargedLeptons.begin())->pdg.Code == particles::e)
//            final_state.electron = *tau_MC.finalStateChargedLeptons.begin();
            analysis::GenParticlePtrVector tauProducts;
            if (analysis::FindDecayProducts(*tau_MC.origin,analysis::TauElectronDecay,tauProducts,false))
                final_state.electron = tauProducts.at(0);

//            else if(tau_MC.finalStateChargedHadrons.size() >= 1)
            else if(!analysis::IsLeptonicTau(*tau_MC.origin)){
                final_state.tau_jet = &tau_MC;
            }
        }

        if (!final_state.electron || !final_state.tau_jet) return false;
        return true;
    }

    virtual void FillFlatTree(const analysis::SelectionResults& /*selection*/) override
    {
        static const float default_value = ntuple::DefaultFloatFillValueForFlatTree();
        static const float default_int_value = ntuple::DefaultIntegerFillValueForFlatTree();

        const analysis::CandidatePtr& electron = selection.GetElectron();
        const ntuple::Electron& ntuple_electron = electron->GetNtupleObject<ntuple::Electron>();
        const analysis::CandidatePtr& tau = selection.GetTau();
        const ntuple::Tau& ntuple_tau = tau->GetNtupleObject<ntuple::Tau>();

        BaseFlatTreeProducer::FillFlatTree(selection);

        flatTree->channel() = static_cast<int>(analysis::Channel::ETau);
        flatTree->pfRelIso_1() = ntuple_electron.pfRelIso;
        flatTree->mva_1() = ntuple_electron.mvaPOGNonTrig;
        flatTree->passid_1() = true;
        flatTree->passiso_1() = ntuple_electron.pfRelIso < cuts::Htautau_Summer13::ETau::electronID::pFRelIso;
        flatTree->decayMode_1() = default_int_value;
        flatTree->iso_1() = ntuple_electron.pfRelIso;
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


        const bool electron_matched = analysis::HasMatchWithMCParticle(electron->GetMomentum(),
                                                                       selection.eTau_MC.electron,
                                                                       cuts::DeltaR_MC_Match);
        if(electron_matched) {
            const TLorentzVector& momentum = selection.eTau_MC.electron->momentum;
            flatTree->pt_1_MC   () = momentum.Pt()  ;
            flatTree->phi_1_MC  () = momentum.Phi() ;
            flatTree->eta_1_MC  () = momentum.Eta() ;
            flatTree->m_1_MC    () = momentum.M()   ;
            flatTree->pt_1_visible_MC   () = momentum.Pt()  ;
            flatTree->phi_1_visible_MC  () = momentum.Phi() ;
            flatTree->eta_1_visible_MC  () = momentum.Eta() ;
            flatTree->m_1_visible_MC    () = momentum.M()   ;
            flatTree->pdgId_1_MC() = selection.eTau_MC.electron->pdg.ToInteger();
        } else {
            flatTree->pt_1_MC   () = default_value ;
            flatTree->phi_1_MC  () = default_value ;
            flatTree->eta_1_MC  () = default_value ;
            flatTree->m_1_MC    () = default_value ;
            flatTree->pt_1_visible_MC   () = default_value ;
            flatTree->phi_1_visible_MC  () = default_value ;
            flatTree->eta_1_visible_MC  () = default_value ;
            flatTree->m_1_visible_MC    () = default_value ;
            flatTree->pdgId_1_MC() = default_int_value;
        }

        const bool tau_matched = analysis::HasMatchWithMCObject(tau->GetMomentum(), selection.eTau_MC.tau_jet,
                                                                cuts::DeltaR_MC_Match);
        if(tau_matched) {
            const TLorentzVector& momentum = selection.eTau_MC.tau_jet->origin->momentum;
            flatTree->pt_2_MC   () = momentum.Pt()  ;
            flatTree->phi_2_MC  () = momentum.Phi() ;
            flatTree->eta_2_MC  () = momentum.Eta() ;
            flatTree->m_2_MC    () = momentum.M()   ;
            const TLorentzVector& visible_momentum = selection.eTau_MC.tau_jet->visibleMomentum;
            flatTree->pt_2_visible_MC   () = visible_momentum.Pt()  ;
            flatTree->phi_2_visible_MC  () = visible_momentum.Phi() ;
            flatTree->eta_2_visible_MC  () = visible_momentum.Eta() ;
            flatTree->m_2_visible_MC    () = visible_momentum.M()   ;
            flatTree->pdgId_2_MC() = selection.eTau_MC.tau_jet->origin->pdg.ToInteger();
        } else {
            flatTree->pt_2_MC   () = default_value ;
            flatTree->phi_2_MC  () = default_value ;
            flatTree->eta_2_MC  () = default_value ;
            flatTree->m_2_MC    () = default_value ;
            flatTree->pt_2_visible_MC   () = default_value ;
            flatTree->phi_2_visible_MC  () = default_value ;
            flatTree->eta_2_visible_MC  () = default_value ;
            flatTree->m_2_visible_MC    () = default_value ;
            flatTree->pdgId_2_MC() = default_int_value;
        }

        flatTree->Fill();
    }

protected:
    analysis::BaseAnalyzerData baseAnaData;
    analysis::SelectionResults_etau selection;
    std::shared_ptr<analysis::RecoilCorrectionProducer> recoilCorrectionProducer_etau;
    analysis::EventWeights_etau eventWeights;
};

#include "METPUSubtraction/interface/GBRProjectDict.cxx"
