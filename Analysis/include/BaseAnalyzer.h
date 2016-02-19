/*! Definition of BaseAnalyzer class which is the base class for all X->HH->bbTauTau and H->tautau analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <iomanip>
#include <functional>
#include <string>
#include <iostream>

#include "h-tautau/Analysis/include/TreeExtractor.h"
#include "AnalysisTools/Core/include/CutTools.h"
#include "h-tautau/Analysis/include/GenParticle.h"
#include "h-tautau/Analysis/include/MCfinalState.h"
#include "h-tautau/Analysis/include/RunReport.h"
#include "AnalysisTools/Core/include/Tools.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/Analysis/include/FlatTree.h"
#include "AnalysisTools/Core/include/ProgressReporter.h"


#include "h-tautau/Analysis/include/Htautau_Summer13.h"
#include "h-tautau/Analysis/include/Htautau_TriggerEfficiency.h"
#include "BTagWeight.h"
#include "Config.h"
#include "custom_cuts.h"
#include "h-tautau/Analysis/include/MvaMet.h"
#include "h-tautau/Analysis/include/RecoilCorrection.h"
#include "h-tautau/Analysis/include/JetEnergyUncertainty.h"
#include "EventWeights.h"


namespace analysis {

class BaseAnalyzerData : public root_ext::AnalyzerData {
public:
    BaseAnalyzerData(std::shared_ptr<TFile> outputFile) : AnalyzerData(outputFile) {}

    SELECTION_ENTRY(Selection)

    TH1D_ENTRY_FIX(N_objects, 1, 500, -0.5)
    TH1D_ENTRY(Mass, 3000, 0.0, 3000.0)
    TH1D_ENTRY(Htautau_Mass, 60, 0.0, 300.0)
};

class BaseAnalyzer {
public:
    BaseAnalyzer(const std::string& inputFileName, const std::string& outputFileName, const std::string& configFileName,
                 const std::string& _prefix = "none", size_t _maxNumberOfEvents = 0)
        : config(configFileName),
          outputFile(root_ext::CreateRootFile(outputFileName)),
          anaDataBeforeCut(outputFile, "before_cut"), anaDataAfterCut(outputFile, "after_cut"),
          anaDataFinalSelection(outputFile, "final_selection"),
          maxNumberOfEvents(_maxNumberOfEvents),
          mvaMetProducer(config.MvaMet_dZcut(), config.MvaMet_inputFileNameU(), config.MvaMet_inputFileNameDPhi(),
                         config.MvaMet_inputFileNameCovU1(), config.MvaMet_inputFileNameCovU2())
    {
        if ( _prefix != "external" ){
            progressReporter = std::shared_ptr<tools::ProgressReporter>(
                        new tools::ProgressReporter(config.ReportInterval(), std::cout));
            treeExtractor = std::shared_ptr<TreeExtractor>(
                        new TreeExtractor(_prefix == "none" ? "" : _prefix, inputFileName, config.extractMCtruth(),
                                          config.MaxTreeVersion()));

        }
        TH1::SetDefaultSumw2();
        TH1::AddDirectory(kFALSE);
        TH2::AddDirectory(kFALSE);
        if(config.EstimateJetEnergyUncertainties()) {
            jetEnergyUncertaintyCorrector = std::shared_ptr<JetEnergyUncertaintyCorrector>(
                        new JetEnergyUncertaintyCorrector(config.JetEnergyUncertainties_inputFile(),
                                                          config.JetEnergyUncertainties_inputSection()));
        }
    }

    virtual ~BaseAnalyzer() {}

    virtual void Run()
    {
        size_t n = 0;
        auto _event = std::shared_ptr<EventDescriptor>(new EventDescriptor());
        if (!treeExtractor || !progressReporter)
            throw exception("treeExtractor not initialized");
        for(; ( !maxNumberOfEvents || n < maxNumberOfEvents ) && treeExtractor->ExtractNext(*_event); ++n) {
            progressReporter->Report(n);
//            std::cout << "event = " << _event->eventId().eventId << std::endl;
            if(config.RunSingleEvent() && _event->eventId().eventId != config.SingleEventId()) continue;
            ProcessEventWithEnergyUncertainties(_event);
            if(config.RunSingleEvent()) break;
        }
        progressReporter->Report(n, true);
    }

    void ProcessEventWithEnergyUncertainties(std::shared_ptr<const EventDescriptor> _event)
    {
        TryProcessEvent(_event, EventEnergyScale::Central);
        if(config.EstimateTauEnergyUncertainties()) {
            TryProcessEvent(_event, EventEnergyScale::TauUp);
            TryProcessEvent(_event, EventEnergyScale::TauDown);
        }
        if(config.EstimateJetEnergyUncertainties()) {
            TryProcessEvent(_event, EventEnergyScale::JetUp);
            TryProcessEvent(_event, EventEnergyScale::JetDown);
        }
        if(config.EstimateBtagEfficiencyUncertainties()) {
            TryProcessEvent(_event, EventEnergyScale::BtagEfficiencyUp);
            TryProcessEvent(_event, EventEnergyScale::BtagEfficiencyDown);
        }
        if(config.EstimateBtagFakeUncertainties()) {
            TryProcessEvent(_event, EventEnergyScale::BtagFakeUp);
            TryProcessEvent(_event, EventEnergyScale::BtagFakeDown);
        }
    }

private:
    void TryProcessEvent(std::shared_ptr<const EventDescriptor> _event, EventEnergyScale energyScale)
    {
        eventEnergyScale = energyScale;
        scaledTaus = _event->taus();
        scaledJets = _event->jets();
        if(energyScale == EventEnergyScale::TauUp || energyScale == EventEnergyScale::TauDown) {
            const double sign = energyScale == EventEnergyScale::TauUp ? +1 : -1;
            const double sf = 1.0 + sign * cuts::Htautau_Summer13::tauCorrections::energyUncertainty;
            for(ntuple::Tau& tau : scaledTaus) {
                const TLorentzVector momentum = MakeLorentzVectorPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass);
                const TLorentzVector scaled_momentum = momentum * sf;
                tau.pt = scaled_momentum.Pt();
                tau.eta = scaled_momentum.Eta();
                tau.phi = scaled_momentum.Phi();
                tau.mass = scaled_momentum.M();
            }
        } else if(energyScale == EventEnergyScale::JetUp || energyScale == EventEnergyScale::JetDown) {
            const bool scale_up = energyScale == EventEnergyScale::JetUp;
            jetEnergyUncertaintyCorrector->ApplyCorrection(scaledJets, scale_up);
        }

        event = _event;
        GetEventWeights().Reset();
        if(config.ApplyPUreweight())
            GetEventWeights().CalculatePuWeight(event->eventInfo());
        if(config.IsEmbeddedSample())
            GetEventWeights().SetEmbeddedWeight(event->genEvent().embeddedWeight);

        try {
            ProcessEvent();
        } catch(cuts::cut_failed&){}
        GetAnaData().Selection("event").fill_selection(GetEventWeights().GetPartialWeight());
    }

protected:
    virtual BaseAnalyzerData& GetAnaData() = 0;
    virtual RecoilCorrectionProducer& GetRecoilCorrectionProducer() = 0;
    virtual SelectionResults& ApplyBaselineSelection() = 0;
    virtual EventWeights& GetEventWeights() = 0;
    virtual void ProcessEvent() = 0;

    const ntuple::TauVector& GetNtupleTaus() const { return scaledTaus; }
    const ntuple::JetVector& GetNtupleJets() const { return scaledJets; }

    bool HaveTriggerFired(const std::set<std::string>& interestinghltPaths) const
    {
        for (const ntuple::Trigger& trigger : event->triggers()){
            for (size_t n = 0; HaveTriggerMatched(trigger.hltpaths, interestinghltPaths, n); ++n){
                if (trigger.hltresults.at(n) == 1 && trigger.hltprescales.at(n) == 1)
                    return true;
            }
        }
        return false;
    }

    std::vector<std::string> CollectPathsForTriggerFired(const std::set<std::string>& interestinghltPaths)
    {
        std::vector<std::string> firedPaths;
        for (const ntuple::Trigger& trigger : event->triggers()){
            for (size_t n = 0; HaveTriggerMatched(trigger.hltpaths, interestinghltPaths, n); ++n){
                if (trigger.hltresults.at(n) == 1 && trigger.hltprescales.at(n) == 1){
                    std::string firedPath = trigger.hltpaths.at(n);
                    firedPaths.push_back(firedPath);
                }
            }
        }
        return firedPaths;
    }

    template<typename ObjectType, typename BaseSelectorType, typename NtupleObjectType,
             typename ObjectPtrType = std::shared_ptr<const ObjectType>,
             typename Comparitor = std::less<ObjectPtrType> >
    std::vector<ObjectPtrType> CollectObjects(const std::string& selection_label, const BaseSelectorType& base_selector,
                                              const std::vector<NtupleObjectType>& ntuple_objects,
                                              Comparitor comparitor = Comparitor())
    {
        cuts::ObjectSelector& objectSelector = GetAnaData().Selection(selection_label);
        SelectionManager selectionManager(anaDataBeforeCut, selection_label, GetEventWeights().GetPartialWeight());

        const auto selector = [&](size_t id) -> ObjectPtrType {
            ObjectPtrType candidate(new ObjectType(ntuple_objects.at(id)));
            cuts::Cutter cut(&objectSelector);
            base_selector(candidate, selectionManager, cut);
            return candidate;
        };

        const auto selected = objectSelector.collect_objects<ObjectPtrType>(GetEventWeights().GetPartialWeight(),
                                                                            ntuple_objects.size(), selector,
                                                                            comparitor);
        SelectionManager selectionManager_afterCut(anaDataAfterCut, selection_label,
                                                   GetEventWeights().GetPartialWeight());
        for(const auto& candidate : selected) {
            cuts::Cutter cut(nullptr);
            base_selector(candidate, selectionManager_afterCut, cut);
        }
        GetAnaData().N_objects(selection_label).Fill(selected.size(), GetEventWeights().GetPartialWeight());
        GetAnaData().N_objects(selection_label + "_ntuple").Fill(ntuple_objects.size(),
                                                                 GetEventWeights().GetPartialWeight());
        return selected;
    }

    template<typename BaseSelectorMethod, typename NtupleObjectType>
    CandidatePtrVector CollectCandidateObjects(const std::string& selection_label, BaseSelectorMethod selector_method,
                                               const std::vector<NtupleObjectType>& ntuple_objects)
    {
        const auto base_selector = [&](const CandidatePtr& candidate, SelectionManager& selectionManager,
                                       cuts::Cutter& cut)
            { (this->*selector_method)(candidate, selectionManager, cut); };
        return CollectObjects<Candidate>(selection_label, base_selector, ntuple_objects);
    }

    CandidatePtrVector CollectMuons()
    {
        return CollectCandidateObjects("muons", &BaseAnalyzer::SelectMuon, event->muons());
    }

    CandidatePtrVector CollectSignalMuons()
    {
        return CollectCandidateObjects("muons_sgn", &BaseAnalyzer::SelectSignalMuon, event->muons());
    }

    virtual void SelectMuon(const CandidatePtr& candidate, SelectionManager& selectionManager, cuts::Cutter& cut)
    {
        throw std::runtime_error("Muon selection for signal not implemented");
    }

    virtual void SelectSignalMuon(const CandidatePtr& candidate, SelectionManager& selectionManager, cuts::Cutter& cut)
    {
        throw std::runtime_error("Signal muon selection for signal not implemented");
    }

    CandidatePtrVector CollectBackgroundMuons()
    {
        return CollectCandidateObjects("muons_bkg", &BaseAnalyzer::SelectBackgroundMuon, event->muons());
    }

    virtual void SelectBackgroundMuon(const CandidatePtr& muon, SelectionManager& selectionManager, cuts::Cutter& cut)
    {
        using namespace cuts::Htautau_Summer13::muonVeto;
        const ntuple::Muon& object = muon->GetNtupleObject<ntuple::Muon>();

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex->GetPosition().Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = Calculate_dxy(mu_vertex, primaryVertex->GetPosition(), muon->GetMomentum());
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(isGlobalMuonPromptTight) == isGlobalMuonPromptTight, "tight");
        cut(X(isPFMuon) == isPFMuon, "PF");
        cut(X(nMatchedStations) > nMatched_Stations, "stations");
        cut(X(pixHits) > pixHits, "pix_hits");
        cut(X(trackerLayersWithMeasurement) > trackerLayersWithMeasurement, "layers");
        cut(X(pfRelIso) < pfRelIso, "pFRelIso");
    }

    CandidatePtrVector CollectTaus()
    {
        return CollectCandidateObjects("taus", &BaseAnalyzer::SelectTau, correctedTaus);
    }

    CandidatePtrVector CollectSignalTaus()
    {
        return CollectCandidateObjects("taus_sgn", &BaseAnalyzer::SelectSignalTau, correctedTaus);
    }

    virtual void SelectTau(const CandidatePtr& candidate, SelectionManager& selectionManager, cuts::Cutter& cut)
    {
        throw std::runtime_error("Tau selection for signal not implemented");
    }

    virtual void SelectSignalTau(const CandidatePtr& candidate, SelectionManager& selectionManager, cuts::Cutter& cut)
    {
        throw std::runtime_error("Signal tau selection for signal not implemented");
    }

    CandidatePtrVector CollectElectrons()
    {
        return CollectCandidateObjects("electrons", &BaseAnalyzer::SelectElectron, event->electrons());
    }

    CandidatePtrVector CollectSignalElectrons()
    {
        return CollectCandidateObjects("electrons_sgn", &BaseAnalyzer::SelectSignalElectron, event->electrons());
    }

    virtual void SelectElectron(const CandidatePtr& candidate, SelectionManager& selectionManager, cuts::Cutter& cut)
    {
        throw std::runtime_error("Electron selection for signal not implemented");
    }

    virtual void SelectSignalElectron(const CandidatePtr& candidate, SelectionManager& selectionManager,
                                      cuts::Cutter& cut)
    {
        throw std::runtime_error("Electron selection for signal not implemented");
    }

    CandidatePtrVector CollectBackgroundElectrons()
    {
        return CollectCandidateObjects("electrons_bkg", &BaseAnalyzer::SelectBackgroundElectron, event->electrons());
    }

    virtual void SelectBackgroundElectron(const CandidatePtr& electron, SelectionManager& selectionManager,
                                          cuts::Cutter& cut)
    {
        using namespace cuts::Htautau_Summer13::electronVeto;
        const ntuple::Electron& object = electron->GetNtupleObject<ntuple::Electron>();

        cut(true, ">0 ele cand");
        cut(X(pt) > pt, "pt");
        const double eta = std::abs( X(eta) );
        cut(eta < eta_high, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex->GetPosition().Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 ele_vertex(object.vx, object.vy, object.vz);
        // same as dB
        const double d0_PV = Calculate_dxy(ele_vertex, primaryVertex->GetPosition(), electron->GetMomentum());
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");
        const size_t pt_index = object.pt < ref_pt ? 0 : 1;
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig) > MVApogNonTrig[pt_index][eta_index], "mva");
        cut(X(missingHits) < missingHits, "missingHits");
        cut(X(hasMatchedConversion) == hasMatchedConversion, "conversion");
    }

    CandidatePtrVector CollectBJets(const CandidatePtrVector& looseJets, bool doReTag, bool applyCsvCut)
    {
        using namespace cuts::Htautau_Summer13::btag;

        static const std::map<EventEnergyScale, std::pair<int,int>> btag_modes_map = {
            { EventEnergyScale::Central, {0,0} }, { EventEnergyScale::TauUp, {0,0} },
            { EventEnergyScale::TauDown, {0,0} }, { EventEnergyScale::JetUp, {0,0} },
            { EventEnergyScale::JetDown, {0,0} }, { EventEnergyScale::BtagEfficiencyUp, {2,0} },
            { EventEnergyScale::BtagEfficiencyDown, {1,0} }, { EventEnergyScale::BtagFakeUp, {0,2} },
            { EventEnergyScale::BtagFakeDown, {0,1} }
        };

        CandidatePtrVector bjets;
        const std::pair<int,int> btag_pair = btag_modes_map.at(eventEnergyScale);
        const int btag_mode = btag_pair.first;
        const int bfake_mode = btag_pair.second;
        for(const CandidatePtr& looseJetCandidate : looseJets) {
            const ntuple::Jet& looseJet = looseJetCandidate->GetNtupleObject<ntuple::Jet>();
            if(looseJet.pt <= pt || std::abs(looseJet.eta) >= eta) continue;
            if(doReTag && !btag::ReTag(looseJet, btag::payload::EPS13, btag::tagger::CSVM, btag_mode, bfake_mode, CSV))
                continue;
            else if(!doReTag && applyCsvCut && looseJet.combinedSecondaryVertexBJetTags <= CSV)
                continue;

            bjets.push_back(looseJetCandidate);
        }

        const auto bjetsSelector = [&] (const CandidatePtr& first, const CandidatePtr& second) -> bool
        {
            const ntuple::Jet& first_bjet = first->GetNtupleObject<ntuple::Jet>();
            const ntuple::Jet& second_bjet = second->GetNtupleObject<ntuple::Jet>();

            return first_bjet.combinedSecondaryVertexBJetTags > second_bjet.combinedSecondaryVertexBJetTags;
        };

        std::sort(bjets.begin(), bjets.end(), bjetsSelector);
        return bjets;
    }

    CandidatePtrVector CollectLooseJets()
    {
        return CollectCandidateObjects("loose_jets", &BaseAnalyzer::SelectLooseJet, GetNtupleJets());
    }

    virtual void SelectLooseJet(const CandidatePtr& jet, SelectionManager& selectionManager, cuts::Cutter& cut)
    {
        using namespace cuts::Htautau_Summer13::jetID;
        const ntuple::Jet& object = jet->GetNtupleObject<ntuple::Jet>();

        cut(true, ">0 jet cand");
        cut(X(pt) > pt_loose, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const bool passLooseID = passPFLooseId(object);
        cut(Y(passLooseID) == pfLooseID, "pfLooseID");
        const bool passPUlooseID = ntuple::JetID_MVA::PassLooseId(object.puIdBits);
        cut(Y(passPUlooseID) == puLooseID, "puLooseID");
    }

    CandidatePtrVector CollectJets(const CandidatePtrVector& looseJets)
    {
        using namespace cuts::Htautau_Summer13::jetID;
        CandidatePtrVector jets;
        for(const CandidatePtr& looseJet : looseJets) {
            if(looseJet->GetMomentum().Pt() > pt)
                jets.push_back(looseJet);
        }

        const auto jetsSelector = [&] (const CandidatePtr& first_jet, const CandidatePtr& second_jet) -> bool
        {
            const double first_pt = first_jet->GetMomentum().Pt();
            const double second_pt = second_jet->GetMomentum().Pt();
            return first_pt > second_pt;
        };

        std::sort(jets.begin(), jets.end(), jetsSelector);

        return jets;
    }

    VertexPtrVector CollectVertices()
    {
        const auto base_selector = [&](const VertexPtr& vertex, SelectionManager& selectionManager,
                                       cuts::Cutter& cut)
            { SelectVertex(vertex, selectionManager, cut); };
        const auto vertex_comparitor = [&](const VertexPtr& first, const VertexPtr& second) -> bool
            { return *first < *second; };
        return CollectObjects<Vertex>("vertices", base_selector, event->vertices(), vertex_comparitor);
    }

    void SelectVertex(const VertexPtr& vertex, SelectionManager& selectionManager, cuts::Cutter& cut)
    {
        using namespace cuts::Htautau_Summer13::vertex;
        const ntuple::Vertex& object = vertex->GetNtupleObject();

        cut(true, ">0 vertex");
        cut(X(ndf) > ndf, "ndf");
        cut(std::abs( X(z) ) < z, "z");
        const double r_vertex = std::sqrt(object.x*object.x+object.y*object.y);
        cut(std::abs( Y(r_vertex) ) < r, "r");
    }

    CandidatePtrVector FindCompatibleObjects(const CandidatePtrVector& objects1, const CandidatePtrVector& objects2,
                                          double minDeltaR, Candidate::Type type, const std::string& hist_name,
                                          int expectedCharge = Candidate::UnknownCharge())
    {
        CandidatePtrVector result;
        for(const CandidatePtr& object1 : objects1) {
            for(const CandidatePtr& object2 : objects2) {
                if(object2->GetMomentum().DeltaR(object1->GetMomentum()) > minDeltaR) {
                    const CandidatePtr candidate(new Candidate(type, object1, object2));
                    if (expectedCharge != Candidate::UnknownCharge() && candidate->GetCharge() != expectedCharge)
                        continue;
                    result.push_back(candidate);
                    GetAnaData().Mass(hist_name).Fill(candidate->GetMomentum().M(),
                                                      GetEventWeights().GetPartialWeight());
                }
            }
        }
        GetAnaData().N_objects(hist_name).Fill(result.size(), GetEventWeights().GetPartialWeight());
        return result;
    }


    CandidatePtrVector FindCompatibleObjects(const CandidatePtrVector& objects, double minDeltaR, Candidate::Type type,
                                          const std::string& hist_name, int expectedCharge = Candidate::UnknownCharge())
    {
        CandidatePtrVector result;
        for (unsigned n = 0; n < objects.size(); ++n){
            for (unsigned k = n+1; k < objects.size(); ++k){
//                std::cout << "first tau momentum " << objects.at(n).momentum << std::endl;
//                std::cout << "second tau momentum " << objects.at(k).momentum << std::endl;
//                std::cout << "DeltaR " << objects.at(n).momentum.DeltaR(objects.at(k).momentum) << std::endl;
//                std::cout << "first tau charge " << objects.at(n).charge << std::endl;
//                std::cout << "second tau charge " << objects.at(k).charge << std::endl;
                if(objects.at(n)->GetMomentum().DeltaR(objects.at(k)->GetMomentum()) > minDeltaR) {
                    const CandidatePtr candidate(new Candidate(type, objects.at(n), objects.at(k)));
                    if (expectedCharge != Candidate::UnknownCharge() && candidate->GetCharge() != expectedCharge )
                        continue;
                    result.push_back(candidate);
                    GetAnaData().Mass(hist_name).Fill(candidate->GetMomentum().M(),
                                                      GetEventWeights().GetPartialWeight());
                }
            }
        }
        GetAnaData().N_objects(hist_name).Fill(result.size(), GetEventWeights().GetPartialWeight());
        return result;
    }

    CandidatePtrVector FilterCompatibleObjects(const CandidatePtrVector& objectsToFilter,
                                               const CandidatePtr& referenceObject, double minDeltaR)
    {
        CandidatePtrVector result;
        for(const CandidatePtr& filterObject : objectsToFilter) {
            bool allDaughterPassed = true;
            for (const CandidatePtr& daughter : referenceObject->GetFinalStateDaughters()){
                if(filterObject->GetMomentum().DeltaR(daughter->GetMomentum()) <= minDeltaR) {
                    allDaughterPassed = false;
                    break;
                }
            }
            if (allDaughterPassed) result.push_back(filterObject);
        }
        return result;
    }

    ntuple::TauVector ApplyTauCorrections(const VisibleGenObjectVector& hadronic_taus, bool useLegacyCorrections)
    {
        using namespace cuts::Htautau_Summer13::tauCorrections;

        ntuple::TauVector correctedTaus;

        for(const ntuple::Tau& tau : GetNtupleTaus()) {
            TLorentzVector momentum;
            momentum.SetPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass);

            const bool hasMCmatch = FindMatchedObjects(momentum, hadronic_taus, deltaR).size() != 0;
            const double scaleFactor = MomentumScaleFactor(hasMCmatch, momentum.Pt(),
                                   ntuple::tau_id::ConvertToHadronicDecayMode(tau.decayMode), useLegacyCorrections);
            const TLorentzVector correctedMomentum = momentum * scaleFactor;
            ntuple::Tau correctedTau(tau);
            correctedTau.pt = correctedMomentum.Pt();
            correctedTau.eta = correctedMomentum.Eta();
            correctedTau.phi = correctedMomentum.Phi();
            correctedTau.mass = correctedMomentum.M();
            correctedTaus.push_back(correctedTau);
        }
        return correctedTaus;
    }

    ntuple::MET ApplyTauCorrectionsToMVAMET(const ntuple::MET& metMVA, const ntuple::TauVector& correctedTaus)
    {
        TLorentzVector sumCorrectedTaus, sumTaus;
        for(const ntuple::Tau& tau : GetNtupleTaus()) {
            TLorentzVector momentum;
            momentum.SetPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass);
            sumTaus += momentum;
        }

        for (const ntuple::Tau& correctedTau : correctedTaus) {
            TLorentzVector correctedMomentum;
            correctedMomentum.SetPtEtaPhiM(correctedTau.pt, correctedTau.eta, correctedTau.phi,correctedTau.mass);
            sumCorrectedTaus += correctedMomentum;
        }

        TLorentzVector met, metCorrected;
        met.SetPtEtaPhiM(metMVA.pt, 0, metMVA.phi, 0.);
        metCorrected = met + sumTaus - sumCorrectedTaus;
        ntuple::MET correctedMET = metMVA;
        correctedMET.pt = metCorrected.Pt();
        correctedMET.phi = metCorrected.Phi();
        return correctedMET;
    }

    CandidatePtrVector ApplyTriggerMatch(const CandidatePtrVector& higgses, const std::set<std::string>& hltPaths,
                                         bool useStandardTriggerMatch)
    {
        CandidatePtrVector triggeredHiggses;
        for (const auto& higgs : higgses){
            if(!useStandardTriggerMatch && HaveTriggerMatched(event->triggerObjects(), hltPaths, *higgs,
                                                                        cuts::Htautau_Summer13::DeltaR_triggerMatch))
                triggeredHiggses.push_back(higgs);
            if (useStandardTriggerMatch && HaveTriggerMatched(hltPaths, *higgs))
                triggeredHiggses.push_back(higgs);
        }
        return triggeredHiggses;
    }

    ntuple::MET ApplyRecoilCorrections(const CandidatePtr& higgs, const GenParticle* resonance,
                                       const size_t njets, const ntuple::MET& correctedMET)
    {
        if (config.ApplyRecoilCorrection()){
            if(!resonance && config.ApplyRecoilCorrectionForW())
                resonance = FindWboson();
            if(!resonance && config.ApplyRecoilCorrectionForZ()){
                GenParticlePtrVector ZProducts;
                bool ztt;
                resonance = FindZboson(ZProducts,ztt);
            }
            if(resonance)
                return GetRecoilCorrectionProducer().ApplyCorrection(correctedMET, higgs->GetMomentum(),
                                                                     resonance->momentum, njets);
        }
        return correctedMET;
    }

    const GenParticle* FindWboson()
    {
        static const particles::ParticleCodes Wcode = { particles::W_plus };
        static const particles::ParticleCodes WDecay_tau = { particles::tau, particles::nu_tau };
        static const particles::ParticleCodes WDecay_electron = { particles::e, particles::nu_e };
        static const particles::ParticleCodes WDecay_muon = { particles::mu, particles::nu_mu };

        const GenParticleSet Wparticles_all = genEvent.GetParticles(Wcode);

        GenParticleSet Wparticles;
        for(const GenParticle* w : Wparticles_all) {
            if(w->mothers.size() == 1) {
                const GenParticle* mother = w->mothers.at(0);
                if(mother->pdg.Code == particles::W_plus && mother->status == particles::HardInteractionProduct)
                    Wparticles.insert(w);
            }
         }

        if (Wparticles.size() == 0) return nullptr;

        if (Wparticles.size() > 1)
            throw exception("more than 1 W in the event");

        const GenParticle* Wboson = *Wparticles.begin();
        while(Wboson->daughters.size() == 1 && Wboson->daughters.front()->pdg.Code == particles::W_plus)
            Wboson = Wboson->daughters.front();

        GenParticlePtrVector WProducts;
        if(FindDecayProducts(*Wboson, WDecay_tau, WProducts,true) ||
                FindDecayProducts(*Wboson, WDecay_electron, WProducts,true) ||
                FindDecayProducts(*Wboson, WDecay_muon, WProducts,true))
            return Wboson;

        throw exception("not leptonic W decay");

    }

    const GenParticle* FindZboson(GenParticlePtrVector& ZProducts, bool& ztt)
    {
        static const particles::ParticleCodes Zcode = { particles::Z };
        static const particles::ParticleCodes ZDecay_electrons = { particles::e, particles::e };
        static const particles::ParticleCodes ZDecay_muons = { particles::mu, particles::mu };
        static const particles::ParticleCodes ZDecay_taus = { particles::tau, particles::tau };

        const GenParticleSet Zparticles_all = genEvent.GetParticles(Zcode);

        GenParticleSet Zparticles;
        for(const GenParticle* z : Zparticles_all) {
            const bool is_hard_interaction_z = z->mothers.size() == 1 && z->mothers.front()->pdg.Code == particles::Z
                    && z->mothers.front()->status == particles::HardInteractionProduct;
            const bool is_pp_z = z->mothers.size() == 2 && z->mothers.at(0)->pdg.Code == particles::p
                    && z->mothers.at(1)->pdg.Code == particles::p;
            if(is_hard_interaction_z || is_pp_z)
                Zparticles.insert(z);
         }

        if (Zparticles.size() > 1 || Zparticles.size() == 0)
            throw exception("not 1 Z per event");

        const GenParticle* Z_mc = *Zparticles.begin();
        while(Z_mc->daughters.size() == 1 && Z_mc->daughters.front()->pdg.Code == particles::Z)
            Z_mc = Z_mc->daughters.front();

        ztt = FindDecayProducts(*Z_mc, ZDecay_taus, ZProducts, true);
        if (ztt || FindDecayProducts(*Z_mc, ZDecay_electrons, ZProducts, true)
                 || FindDecayProducts(*Z_mc, ZDecay_muons, ZProducts, true))
            return Z_mc;
        throw exception("not leptonic Z decay");
    }

    CandidatePtr SelectSemiLeptonicHiggs(const CandidatePtrVector& higgses)
    {
        if(!higgses.size())
            throw std::runtime_error("no available higgs candidate to select");
        const auto higgsSelector = [&] (const CandidatePtr& first, const CandidatePtr& second) -> bool
        {
            const double first_Pt1 = first->GetDaughters().at(0)->GetMomentum().Pt();
            const double first_Pt2 = first->GetDaughters().at(1)->GetMomentum().Pt();
            const double first_sumPt = first_Pt1 + first_Pt2;
            const double second_Pt1 = second->GetDaughters().at(0)->GetMomentum().Pt();
            const double second_Pt2 = second->GetDaughters().at(1)->GetMomentum().Pt();
            const double second_sumPt = second_Pt1 + second_Pt2;

            return first_sumPt < second_sumPt;
        };
        return *std::max_element(higgses.begin(), higgses.end(), higgsSelector) ;
    }

    bool FindAnalysisFinalState(finalState::bbTauTau& final_state)
    {
        static const particles::ParticleCodes resonanceCodes = { particles::MSSM_H, particles::MSSM_A };
        static const particles::ParticleCodes resonanceDecay_1 = { particles::Higgs, particles::Higgs };
        static const particles::ParticleCodes resonanceDecay_2 = { particles::Z, particles::Higgs };
        static const particles::ParticleCodes SM_ResonanceCodes = { particles::Higgs, particles::Z,
                                                                    particles::MSSM_H, particles::MSSM_A };
        static const particles::ParticleCodes SM_ResonanceDecay_1 = { particles::tau, particles::tau };
        static const particles::ParticleCodes SM_ResonanceDecay_2 = { particles::b, particles::b };
        static const particles::ParticleCodes2D HiggsDecays = { SM_ResonanceDecay_1, SM_ResonanceDecay_2 };

        genEvent.Initialize(event->genParticles());
        final_state.Reset();

        const GenParticleSet resonances = genEvent.GetParticles(resonanceCodes);

        if (resonances.size() > 1)
            throw exception("more than 1 resonance per event");

        if (resonances.size() == 1) {
            final_state.resonance = *resonances.begin();

            bool doubleHiggsSignal = true;
            GenParticlePtrVector HiggsBosons;
            if(!FindDecayProducts(*final_state.resonance, resonanceDecay_1, HiggsBosons) &&
                    !FindDecayProducts(*final_state.resonance, resonanceDecay_2, HiggsBosons))
                doubleHiggsSignal = false;

            if(doubleHiggsSignal) {
                GenParticleVector2D HiggsDecayProducts;
                GenParticleIndexVector HiggsIndexes;
                GenParticlePtrVector Higgs_ToTauTau_Product;
                GenParticlePtrVector Higgs_ToBB_Product;
                const bool HH_bbtautau = FindDecayProducts2D(HiggsBosons, HiggsDecays, HiggsDecayProducts, HiggsIndexes);
                if(HH_bbtautau){
                    Higgs_ToTauTau_Product = HiggsDecayProducts.at(0);
                    Higgs_ToBB_Product = HiggsDecayProducts.at(1);
                    final_state.Higgs_TauTau = HiggsBosons.at(HiggsIndexes.at(0));
                    final_state.Higgs_BB = HiggsBosons.at(HiggsIndexes.at(1));
                } else if (FindDecayProducts(*HiggsBosons.at(1),SM_ResonanceDecay_1,Higgs_ToTauTau_Product)){
                    final_state.Higgs_TauTau = HiggsBosons.at(1);
                } else if (FindDecayProducts(*HiggsBosons.at(1),SM_ResonanceDecay_2,Higgs_ToBB_Product)){
                    final_state.Higgs_BB = HiggsBosons.at(1);
                } else
                    throw exception("Nor HH-> bbtautau, nor A->Zh->lltautau, nor A->Zh->llbb");

                for(const GenParticle* tau : Higgs_ToTauTau_Product) {
                    const VisibleGenObject tau_products(tau);
                    final_state.taus.push_back(tau_products);
//                    if(tau_products.finalStateChargedHadrons.size() != 0)
                    if(!IsLeptonicTau(*tau)){

                        final_state.hadronic_taus.push_back(tau_products);
                    }
                }
                for(const GenParticle* b : Higgs_ToBB_Product)
                    final_state.b_jets.push_back(VisibleGenObject(b));

                return HH_bbtautau;
            }
        }

        if(config.ExpectedOneNonSMResonance())
            throw exception("Non-SM resonance not found.");

        //search H->bb, H->tautau
        const GenParticleSet SM_particles = genEvent.GetParticles(SM_ResonanceCodes);

        GenParticlePtrVector SM_ResonanceToTauTau_products;
        GenParticlePtrVector SM_ResonanceToBB_products;

        for (const GenParticle* SM_particle : SM_particles){
            GenParticlePtrVector resonanceDecayProducts;
            if(FindDecayProducts(*SM_particle, SM_ResonanceDecay_1,resonanceDecayProducts)){
                if(!final_state.Higgs_TauTau || (final_state.Higgs_TauTau->pdg.Code != particles::Higgs
                                             && SM_particle->pdg.Code == particles::Higgs)) {
                    final_state.Higgs_TauTau = SM_particle;
                    SM_ResonanceToTauTau_products = resonanceDecayProducts;
                } else if((final_state.Higgs_TauTau->pdg.Code == particles::Higgs
                           && SM_particle->pdg.Code == particles::Higgs)
                          || (final_state.Higgs_TauTau->pdg.Code != particles::Higgs
                              && SM_particle->pdg.Code != particles::Higgs)) {
                    throw exception("more than one SM resonance to tautau per event");
                }
            }
            else if (FindDecayProducts(*SM_particle, SM_ResonanceDecay_2,resonanceDecayProducts)){
                if(!final_state.Higgs_BB || (final_state.Higgs_BB->pdg.Code != particles::Higgs
                                             && SM_particle->pdg.Code == particles::Higgs)) {
                    final_state.Higgs_BB = SM_particle;
                    SM_ResonanceToBB_products = resonanceDecayProducts;
                } else if((final_state.Higgs_BB->pdg.Code == particles::Higgs
                           && SM_particle->pdg.Code == particles::Higgs)
                          || (final_state.Higgs_BB->pdg.Code != particles::Higgs
                              && SM_particle->pdg.Code != particles::Higgs)) {
                    throw exception("more than one SM resonance to bb per event");
                }
            }
        }

        for(const GenParticle* tau : SM_ResonanceToTauTau_products) {
            const VisibleGenObject tau_products(tau);
            final_state.taus.push_back(tau_products);
            if(!IsLeptonicTau(*tau))
                final_state.hadronic_taus.push_back(tau_products);
        }

        for(const GenParticle* b : SM_ResonanceToBB_products)
            final_state.b_jets.push_back(VisibleGenObject(b));

        if (!final_state.Higgs_TauTau && !final_state.Higgs_BB) {
            if(config.ExpectedAtLeastOneSMResonanceToTauTauOrToBB())
                throw exception("SM resonance to tautau or to bb not found.");
            return false;
        }

        return true;
    }

    ntuple::EventType DoEventCategorization(const CandidatePtr& higgs)
    {
        using namespace cuts::Htautau_Summer13::DrellYannCategorization;
        if(!config.DoZEventCategorization())
            return ntuple::EventType::Unknown;

        GenParticlePtrVector ZProducts;
        bool ztt;
        FindZboson(ZProducts,ztt);

        static const particles::ParticleCodes light_lepton_codes = { particles::e, particles::mu };

        const GenParticleSet light_leptons = genEvent.GetParticles(light_lepton_codes, minimal_genParticle_pt);
        const CandidatePtrVector hadronic_taus = higgs->GetDaughters(Candidate::Type::Tau);

        size_t n_hadronic_matches = 0, n_leptonic_matches = 0;
        for(const CandidatePtr& reco_tau : hadronic_taus) {

            for(const GenParticle* gen_product : ZProducts) {
                const VisibleGenObject visible_gen_object(gen_product);
//                 std::cout <<  "GenVisibleTau: " << visible_gen_object.visibleMomentum <<
//                              "; NofLeptons: " << visible_gen_object.finalStateChargedLeptons.size() <<
//                             "; GenTauOrigin: " << visible_gen_object.origin->momentum <<
//                               "; pdg: " << visible_gen_object.origin->pdg.Code.Name() << ", status= " <<
//                               visible_gen_object.origin->status << std::endl;
                if(gen_product->pdg.Code != particles::tau || IsLeptonicTau(*gen_product) ||
                        visible_gen_object.visibleMomentum.Pt() <= minimal_visible_momentum) continue;
                if(HasMatchWithMCObject(reco_tau->GetMomentum(), &visible_gen_object, deltaR_matchGenParticle, true)) {
                    ++n_hadronic_matches;
                    break;
                }
            }

            for(const GenParticle* gen_product : light_leptons) {
                if(HasMatchWithMCParticle(reco_tau->GetMomentum(), gen_product, deltaR_matchGenParticle)) {
                    ++n_leptonic_matches;
                    break;
                }
            }
        }
        //genEvent.Print();

        if(ztt && n_hadronic_matches == hadronic_taus.size()) return ntuple::EventType::ZTT;
        if(n_leptonic_matches) return ztt ? ntuple::EventType::ZTT_L : ntuple::EventType::ZL;
        return ntuple::EventType::ZJ;
    }

    bool GenFilterForZevents(const finalState::bbTauTau& final_state)
    {
        using namespace cuts::Htautau_Summer13::DYEmbedded;
        if (final_state.taus.size() != 2)
            throw exception("not 2 taus in the event at Gen Level");
        const GenParticle* firstTau = final_state.taus.at(0).origin;
        const GenParticle* secondTau = final_state.taus.at(1).origin;
        if ((firstTau->momentum + secondTau->momentum).M() > invariantMassCut) return true;
        return false;
    }

    kinematic_fit::four_body::FitResults RunKinematicFit(const CandidatePtrVector& bjets,
                                                         const Candidate& higgs_to_taus, const ntuple::MET& met)
    {
        using namespace kinematic_fit::four_body;
        using namespace cuts::Htautau_Summer13;

        if(bjets.size() < 2)
            return FitResults();

        const TLorentzVector met_momentum = MakeLorentzVectorPtEtaPhiM(met.pt, 0, met.phi, 0);
        const TMatrix met_cov = ntuple::VectorToSignificanceMatrix(met.significanceMatrix);

        const FitInput input(bjets.at(0)->GetMomentum(), bjets.at(1)->GetMomentum(),
                             higgs_to_taus.GetDaughters().at(0)->GetMomentum(),
                             higgs_to_taus.GetDaughters().at(1)->GetMomentum(),
                             met_momentum, met_cov);
        return Fit(input);
    }

    ntuple::MET ComputeMvaMet(const CandidatePtr& higgs, const VertexPtrVector& goodVertices)
    {
        CandidatePtrVector originalDaughters;
        for(const CandidatePtr& daughter : higgs->GetDaughters()) {
            CandidatePtr original = daughter;
            if(daughter->GetType() == Candidate::Type::Tau) {
                const ntuple::Tau* ntuple_tau = &daughter->GetNtupleObject<ntuple::Tau>();
                const auto position = ntuple_tau - &(*correctedTaus.begin());
                original = CandidatePtr(new Candidate(GetNtupleTaus().at(position)));
            }
            originalDaughters.push_back(original);
        }
        CandidatePtr originalHiggs(new Candidate(higgs->GetType(), originalDaughters.at(0), originalDaughters.at(1)));
        return mvaMetProducer.ComputeMvaMet(originalHiggs, event->pfCandidates(), GetNtupleJets(), primaryVertex,
                                            goodVertices);
    }

protected:
    Config config;
    std::shared_ptr<tools::ProgressReporter> progressReporter;
    std::shared_ptr<const EventDescriptor> event;
    std::shared_ptr<TreeExtractor> treeExtractor;
    std::shared_ptr<TFile> outputFile;
    root_ext::AnalyzerData anaDataBeforeCut, anaDataAfterCut, anaDataFinalSelection;
    size_t maxNumberOfEvents;
    GenEvent genEvent;
    VertexPtr primaryVertex;
    MvaMetProducer mvaMetProducer;
    ntuple::TauVector correctedTaus;
    EventEnergyScale eventEnergyScale;
    std::shared_ptr<JetEnergyUncertaintyCorrector> jetEnergyUncertaintyCorrector;

private:
    ntuple::TauVector scaledTaus;
    ntuple::JetVector scaledJets;
};

} // analysis

namespace make_tools {
template<typename T>
struct Factory {
    static T* Make(int argc, char *argv[])
    {
        std::cout << "Command line: ";
        for(int n = 0; n < argc; ++n)
            std::cout << argv[n] << " ";
        std::cout << std::endl;
        if(argc < 4 || argc > 6)
            throw std::runtime_error("Invalid number of command line arguments.");

        int n = 0;
        const std::string inputFileName = argv[++n];
        const std::string outputFileName = argv[++n];
        const std::string configFileName = argv[++n];
        if(argc <= n)
            return new T(inputFileName, outputFileName, configFileName);

        const std::string prefix = argv[++n];
        if(argc <= n)
            return new T(inputFileName, outputFileName, configFileName, prefix);

        char c;
        size_t maxNumberOfEvents;
        std::istringstream ss_nEvents(argv[++n]);
        ss_nEvents >> c >> maxNumberOfEvents;
        if(c != '@')
            throw std::runtime_error("Bad command line format.");

        return new T(inputFileName, outputFileName, configFileName, prefix, maxNumberOfEvents);
    }
};
} // make_tools
