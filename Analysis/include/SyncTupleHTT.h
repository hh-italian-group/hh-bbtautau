/*! Definition of a tuple with all event information that is required at the analysis level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "h-tautau/Core/include/TauIdResults.h"
#include "hh-bbtautau/Analysis/include/MvaReader.h"

#define LVAR(type, name, pref) VAR(type, name##_##pref)
#define JVAR(type, name, suff, pref) VAR(type, suff##name##_##pref)

#define LEG_DATA(pref) \
    LVAR(Float_t, pt, pref) \
    LVAR(Float_t, pt_tau_ES_up, pref) \
    LVAR(Float_t, pt_tau_ES_down, pref) \
    LVAR(Float_t, phi, pref) \
    LVAR(Float_t, eta, pref) \
    LVAR(Float_t, E, pref) \
    LVAR(Float_t, m, pref) \
    LVAR(Float_t, q, pref) \
    LVAR(Float_t, d0, pref) /* is dxy between leading track and first PV */ \
    LVAR(Float_t, dZ, pref) /* dZ between leading track and first PV */ \
    LVAR(Float_t, mt, pref) /* Use MVAMet sqrt(2 * l_pt * met_pt * (1 - cos( d_phi(l, met) )) */ \
    LVAR(Float_t, pfmt, pref) /* As above but using PF Met */ \
    LVAR(Float_t, puppimt, pref) /* As above but using Puppi Met */ \
    LVAR(Float_t, iso, pref) /* iso */ \
    LVAR(Float_t, id_e_mva_nt_loose, pref) /* Non-triggering electron ID MVA score */ \
    LVAR(Float_t, gen_match, pref) /* Type of gen particle matched to object (see HiggsToTauTauWorking2016#MC_Matching) */ \
    LVAR(Float_t, againstElectronLooseMVA6, pref) \
    LVAR(Float_t, againstElectronMediumMVA6, pref) \
    LVAR(Float_t, againstElectronTightMVA6, pref) \
    LVAR(Float_t, againstElectronVLooseMVA6, pref) \
    LVAR(Float_t, againstElectronVTightMVA6, pref) \
    LVAR(Float_t, againstMuonLoose3, pref) \
    LVAR(Float_t, againstMuonTight3, pref) \
    LVAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw3Hits, pref) \
    LVAR(Float_t, byIsolationMVArun2v1DBoldDMwLTraw, pref) \
    LVAR(Float_t, byIsolationMVArun2017v2DBoldDMwLTraw2017, pref) \
    LVAR(Float_t, chargedIsoPtSum, pref) \
    LVAR(Float_t, decayModeFindingOldDMs, pref) \
    LVAR(Float_t, neutralIsoPtSum, pref) \
    LVAR(Float_t, puCorrPtSum, pref) \
    LVAR(Float_t, trigweight, pref) \
    LVAR(Float_t, idisoweight, pref) \
    /**/

#define JET_DATA(suff, pref) \
    JVAR(Float_t, pt, suff, pref) \
    JVAR(Float_t, eta, suff, pref) \
    JVAR(Float_t, phi, suff, pref) \
    JVAR(Float_t, rawf, suff, pref) /* factor to be applied to the jet p4 to obtain its uncorrected p4 */ \
    JVAR(Float_t, csv, suff, pref) \
    JVAR(Float_t, deepcsv, suff, pref) \
    JVAR(Float_t, deepflavour, suff, pref) \
    JVAR(Float_t, pu_id_raw, suff, pref) \
    JVAR(Float_t, resolution, suff, pref) /* Jet energy resolution in percentage */ \
    JVAR(Float_t, pt_tau_ES_up, suff, pref) \
    JVAR(Float_t, pt_tau_ES_down, suff, pref) \
    JVAR(Float_t, pt_jet_ES_up, suff, pref) \
    JVAR(Float_t, pt_jet_ES_down, suff, pref) \
    /**/

#define SYNC_DATA() \
    VAR(UInt_t, run) /* Run */ \
    VAR(UInt_t, lumi) /* Lumi */ \
    VAR(ULong64_t, evt) /* Evt */ \
    VAR(UInt_t, sampleId) /* sample id */ \
    /* Event Variables */ \
    VAR(Int_t, npv) /* Number of offline primary vertices */ \
    VAR(Float_t, npu) /* Number of in-time pu interactions added to the event */ \
    VAR(Float_t, rho) /* Use fixedGridRhoFastjetAll */ \
    /* Signal leptons */ \
    LEG_DATA(1) /* Leg 1 (leading tau for tt, electon for et,em muon for mt) */ \
    LEG_DATA(2) /* Leg 2 (trailing tau for tt, tau for et, mt, muon for em) */ \
    /* di-tau system */ \
    VAR(Float_t, pt_tt) /* like HIG-13-004 (p4_l1+p4_l2+p4_MET)_T, use PF met */  \
    VAR(Float_t, mt_tot) /* Use MVA MET (see HiggsToTauTauWorking2016#Synchronisation_Ntuple) */  \
    VAR(Float_t, m_vis) \
    VAR(Float_t, m_sv) /* using MarkovChain MC integration svfit.mass() */ \
    VAR(Float_t, m_sv_tau_ES_up) /*  */ \
    VAR(Float_t, m_sv_tau_ES_down) /*  */ \
    VAR(Float_t, m_sv_jet_ES_up) /*  */ \
    VAR(Float_t, m_sv_jet_ES_down) /*  */ \
    VAR(Float_t, mt_sv) /* using MarkovChain MC integration svfit.transverseMass() */ \
    /* MT2 */ \
    VAR(Double_t, mt2) /* MT2 */ \
    VAR(Double_t, mt2_tau_ES_up) /*  */ \
    VAR(Double_t, mt2_tau_ES_down) /*  */ \
    VAR(Double_t, mt2_jet_ES_up) /*  */ \
    VAR(Double_t, mt2_jet_ES_down) /*  */ \
    /* MET */ \
    VAR(Float_t, met) /* type 1 corrected PF MET */  \
    VAR(Float_t, met_tau_ES_up) /* */  \
    VAR(Float_t, met_tau_ES_down) /* */  \
    VAR(Float_t, met_jet_ES_up) /* */  \
    VAR(Float_t, met_jet_ES_down) /* */  \
    VAR(Float_t, metphi) /* type 1 corrected PF MET phi */  \
    VAR(Float_t, metphi_tau_ES_up) /* */  \
    VAR(Float_t, metphi_tau_ES_down) /* */  \
    VAR(Float_t, metphi_jet_ES_up) /* */  \
    VAR(Float_t, metphi_jet_ES_down) /* */  \
    VAR(Float_t, puppimet) /* Puppi corrrected MET */  \
    VAR(Float_t, puppimetphi) /* Puppi corrected MET phi */  \
    VAR(Float_t, mvamet) /* mva met */  \
    VAR(Float_t, mvametphi) /* mva met phi */  \
    VAR(Float_t, pzetavis) /* see HiggsToTauTauWorking2016#Synchronisation_Ntuple */  \
    VAR(Float_t, pzetamiss) /* use MVA met, see HiggsToTauTauWorking2016#Synchronisation_Ntuple */  \
    VAR(Float_t, pfpzetamiss) /* As above but using pf met */  \
    VAR(Float_t, puppipzetamiss) /* As above but using puppi met */  \
    VAR(Float_t, mvacov00) /* mva met */  \
    VAR(Float_t, mvacov01) /* mva met */  \
    VAR(Float_t, mvacov10) /* mva met */  \
    VAR(Float_t, mvacov11) /* mva met */  \
    VAR(Float_t, metcov00) /* pf met */  \
    VAR(Float_t, metcov01) /* pf met */  \
    VAR(Float_t, metcov10) /* pf met */  \
    VAR(Float_t, metcov11) /* pf met */  \
    /* VBF system (Only fill if njetspt20>=2) */ \
    VAR(Float_t, mjj) /* (jet_1->vector()+jet_2->vector() ).M() */ \
    VAR(Float_t, jdeta) /* delta eta between leading jet and subleading jet */ \
    VAR(Float_t, njetingap) /* Number of jets passing pfJetID and pt > 30 GeV, in pseudorapidity gap between the jets */ \
    VAR(Float_t, njetingap20) /* Number of jets passing pfJetID and pt > 20 GeV, in pseudorapidity gap between the jets */ \
    VAR(Float_t, jdphi) /* delta phi between leading jet and subleading jet */ \
    /* additional jets */ \
    VAR(Int_t, njets) /* pt>30 and abs(eta)<4.7 */ \
    VAR(Int_t, njetspt20) /* pt>20 and abs(eta)<4.7 */ \
    VAR(Int_t, njets_vbf) /* pt>10 and no eta cut for vbf selection */ \
    VAR(Bool_t, isVBF) /* Event is vbf if there are 2 additional jets which satisfy the VBF selection */ \
    JET_DATA(j, 1) /* leading jet sorted by pt (Fill only if corrected jet pt > 20 GeV) */ \
    JET_DATA(j, 2) /* trailing jet sorted by pt (Fill only if corrected jet pt>20 GeV) */ \
    JET_DATA(j, vbf_1) /* leading jet sorted by pt (Fill only if corrected jet pt > 20 GeV) */ \
    JET_DATA(j, vbf_2) /* trailing jet sorted by pt (Fill only if corrected jet pt>20 GeV) */ \
    JET_DATA(b, 1) /* leading b-jet sorted by pt (Fill only if corrected b-jet pt>20 GeV) */ \
    JET_DATA(b, 2) /* leading b-jet sorted by pt (Fill only if corrected b-jet pt>20 GeV) */ \
    /* Extra lepton vetos */ \
    VAR(Bool_t, extraelec_veto) /* Event is vetoed by the extra electron veto if true */ \
    VAR(Bool_t, extramuon_veto) /* Event is vetoed by the extra muon veto if true */ \
    VAR(Float_t, puweight) \
    VAR(Float_t, leptonidisoWeight) \
    VAR(Float_t, leptontrigWeight) \
    VAR(Float_t, final_weight) \
    VAR(Float_t, dy_weight) \
    /* hh->bbtautau part */ \
    VAR(Float_t, shapeWeight) /* genWeight * puWeight * genEventSpec */ \
    VAR(Float_t, topWeight) /* gen top pt weight for TTbar */ \
    VAR(Float_t, btagWeight) /* b tag weight */ \
    VAR(Int_t, lhe_n_partons) \
    VAR(Int_t, lhe_n_b_partons) \
    VAR(Float_t, lhe_HT) \
    VAR(Int_t, nbjets) /* pt>30 and abs(eta)<2.4 */ \
    JET_DATA(bjet_, 1) /* leading b-jet sorted by csv (Fill only if corrected b-jet pt>20 GeV) */ \
    JET_DATA(bjet_, 2) /* leading b-jet sorted by csv (Fill only if corrected b-jet pt>20 GeV) */ \
    VAR(Double_t, ht_other_jets) /* Ht of all jets in the event except the first 2 jets */\
    /* mva_score */ \
    VAR(Double_t, mva_score_nonRes_kl1) /* mva_score */\
    VAR(Double_t, mva_score_lm_320) /* mva_score */\
    VAR(Double_t, mva_score_mm_400) /* mva_score */\
    VAR(Double_t, mva_score_hm_650) /* mva_score */\
    VAR(Double_t, mva_score_nonRes_kl1_tau_ES_up) /* mva_score */\
    VAR(Double_t, mva_score_lm_320_tau_ES_up) /* mva_score */\
    VAR(Double_t, mva_score_mm_400_tau_ES_up) /* mva_score */\
    VAR(Double_t, mva_score_hm_650_tau_ES_up) /* mva_score */\
    VAR(Double_t, mva_score_nonRes_kl1_tau_ES_down) /* mva_score */\
    VAR(Double_t, mva_score_lm_320_tau_ES_down) /* mva_score */\
    VAR(Double_t, mva_score_mm_400_tau_ES_down) /* mva_score */\
    VAR(Double_t, mva_score_hm_650_tau_ES_down) /* mva_score */\
    VAR(Double_t, mva_score_nonRes_kl1_jet_ES_up) /* mva_score */\
    VAR(Double_t, mva_score_lm_320_jet_ES_up) /* mva_score */\
    VAR(Double_t, mva_score_mm_400_jet_ES_up) /* mva_score */\
    VAR(Double_t, mva_score_hm_650_jet_ES_up) /* mva_score */\
    VAR(Double_t, mva_score_nonRes_kl1_jet_ES_down) /* mva_score */\
    VAR(Double_t, mva_score_lm_320_jet_ES_down) /* mva_score */\
    VAR(Double_t, mva_score_mm_400_jet_ES_down) /* mva_score */\
    VAR(Double_t, mva_score_hm_650_jet_ES_down) /* mva_score */\
    /* m_kinFit */ \
    VAR(Float_t, m_kinfit) \
    VAR(Float_t, m_kinfit_tau_ES_up) \
    VAR(Float_t, m_kinfit_tau_ES_down) \
    VAR(Float_t, m_kinfit_jet_ES_up) \
    VAR(Float_t, m_kinfit_jet_ES_down) \
    VAR(Int_t, kinfit_convergence) \
    VAR(Float_t, deltaR_ll) \
    VAR(UInt_t, nFatJets) \
    VAR(Int_t, hasFatJet) \
    VAR(Float_t, fatJet_pt) \
    VAR(Float_t, fatJet_eta) \
    VAR(Float_t, fatJet_phi) \
    VAR(Float_t, fatJet_energy) \
    VAR(Float_t, fatJet_m_filtered) \
    VAR(Float_t, fatJet_m_trimmed) \
    VAR(Float_t, fatJet_m_softDrop) \
    VAR(Int_t, fatJet_n_subjets) \
    VAR(Float_t, fatJet_n_subjettiness_tau1) \
    VAR(Float_t, fatJet_n_subjettiness_tau2) \
    VAR(Float_t, fatJet_n_subjettiness_tau3) \
    VAR(UInt_t, genJets_nTotal) \
    VAR(UInt_t, genJets_nStored) \
    VAR(UInt_t, genJets_nStored_hadronFlavour_b) \
    VAR(UInt_t, genJets_nStored_hadronFlavour_c) \
    VAR(UInt_t, jets_nTotal_hadronFlavour_b) \
    VAR(UInt_t, jets_nTotal_hadronFlavour_c) \
    VAR(UInt_t, jets_nSelected_hadronFlavour_b) \
    VAR(UInt_t, jets_nSelected_hadronFlavour_c) \
    /* mva variable */ \
    VAR(Double_t, pt_hbb) \
    VAR(Double_t, pt_l1l2) \
    VAR(Double_t, pt_htautau) \
    VAR(Double_t, pt_htautau_sv) \
    VAR(Double_t, pt_MET) \
    VAR(Double_t, HT_otherjets) \
    VAR(Double_t, p_zeta) \
    VAR(Double_t, p_zetavisible) \
    VAR(Double_t, dphi_l1l2) \
    VAR(Double_t, abs_dphi_b1b2) \
    VAR(Double_t, dphi_b1b2) \
    VAR(Double_t, dphi_l1MET) \
    VAR(Double_t, abs_dphi_METhtautau_sv) \
    VAR(Double_t, dphi_METhtautau_sv) \
    VAR(Double_t, dphi_hbbMET) \
    VAR(Double_t, abs_dphi_hbbhatutau_sv) \
    VAR(Double_t, abs_deta_b1b2) \
    VAR(Double_t, abs_deta_l2MET) \
    VAR(Double_t, abs_deta_hbbMET) \
    VAR(Double_t, dR_l1l2) \
    VAR(Double_t, dR_hbbMET) \
    VAR(Double_t, dR_hbbhtautau) \
    VAR(Double_t, dR_l1l2Pt_htautau) \
    VAR(Double_t, dR_l1l2Pt_htautau_sv) \
    VAR(Double_t, MT_l1) \
    VAR(Double_t, MT_htautau) \
    VAR(Double_t, MT_htautau_sv) \
    VAR(Double_t, MT_tot) \
    VAR(Double_t, MT2) \
    VAR(Double_t, mass_top1) \
    VAR(Double_t, mass_X) \
    VAR(Double_t, mass_H) \
    VAR(Double_t, mass_H_sv) \
    VAR(Double_t, mass_H_vis) \
    VAR(Double_t, mass_H_kinfit_chi2) \
    VAR(Double_t, phi_sv) \
    VAR(Double_t, phi_1_sv) \
    VAR(Double_t, phi_2_sv) \
    VAR(Double_t, costheta_METhtautau_sv) \
    VAR(Double_t, costheta_METhbb) \
    VAR(Double_t, costheta_b1hbb) \
    VAR(Double_t, costheta_htautau_svhhMET) \
    VAR(double, prefiringweight) \
    VAR(double, prefiringweightup) \
    VAR(double, prefiringweightdown) \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(htt_sync, SyncEvent, SyncTuple, SYNC_DATA, "events")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(htt_sync, SyncTuple, SYNC_DATA)
#undef VAR
#undef SYNC_DATA
#undef LEG_DATA
#undef LVAR
#undef JET_DATA
#undef JVAR

namespace htt_sync {


void FillSyncTuple(analysis::EventInfo& event, htt_sync::SyncTuple& sync, analysis::Period run_period,
                   bool apply_svFit,
                   double weight,
                   //double dy_weight,
                   analysis::mva_study::MvaReader* mva_reader = nullptr,
                   analysis::EventInfo* event_tau_up = nullptr,
                   analysis::EventInfo* event_tau_down = nullptr,
                   analysis::EventInfo* event_jet_up = nullptr,
                   analysis::EventInfo* event_jet_down = nullptr);
}

#undef COND_VAL
#undef COND_VAL_INT
