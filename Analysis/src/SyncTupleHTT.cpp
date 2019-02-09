/*! Definition of a tuple with all event information that is required at the analysis level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "hh-bbtautau/Analysis/include/SyncTupleHTT.h"

#define COND_VAL(cond, val) cond ? static_cast<float>(val) : default_value
#define COND_VAL_INT(cond, val) cond ? static_cast<int>(val) : default_int_value

namespace htt_sync {

void FillSyncTuple(analysis::EventInfoBase& event, htt_sync::SyncTuple& sync, analysis::Period run_period,
                   analysis::mva_study::MvaReader* mva_reader,
                   analysis::EventInfoBase* event_tau_up,
                   analysis::EventInfoBase* event_tau_down,
                   analysis::EventInfoBase* event_jet_up,
                   analysis::EventInfoBase* event_jet_down)
{

    static constexpr float default_value = std::numeric_limits<float>::lowest();
    static constexpr int default_int_value = std::numeric_limits<int>::lowest();

    sync().run = event->run;
    sync().lumi = event->lumi;
    sync().evt = event->evt;
    sync().sampleId = event->file_desc_id;
    sync().npv = event->npv;
    sync().npu = event->npu;

    sync().pt_1 = event->p4_1.Pt();
    sync().pt_tau_ES_up_1 = COND_VAL(event_tau_up, (*event_tau_up)->p4_1.Pt());
    sync().pt_tau_ES_down_1 = COND_VAL(event_tau_down, (*event_tau_down)->p4_1.Pt());
    sync().phi_1 = event->p4_1.Phi();
    sync().eta_1 = event->p4_1.Eta();
    sync().m_1 = event->p4_1.mass();
    sync().q_1 = event->q_1;
    sync().d0_1 = event->dxy_1;
    sync().dZ_1 = event->dz_1;
    sync().pfmt_1 = static_cast<float>(analysis::Calculate_MT(event->p4_1, event->pfMET_p4));
    sync().iso_1 =  event->iso_1;
    sync().gen_match_1 = event->gen_match_1;

    const analysis::TauIdResults tauId_1(event->tauId_flags_1);
    sync().againstElectronLooseMVA6_1 = tauId_1.Result("againstElectronLooseMVA6");
    sync().againstElectronMediumMVA6_1 = tauId_1.Result("againstElectronMediumMVA6");
    sync().againstElectronTightMVA6_1 = tauId_1.Result("againstElectronTightMVA6");
    sync().againstElectronVLooseMVA6_1 = tauId_1.Result("againstElectronVLooseMVA6");
    sync().againstElectronVTightMVA6_1 = tauId_1.Result("againstElectronVTightMVA6");
    sync().againstMuonLoose3_1 = tauId_1.Result("againstMuonLoose3");
    sync().againstMuonTight3_1 = tauId_1.Result("againstMuonTight3");
    sync().byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = event->tauId_byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
    sync().byIsolationMVArun2v1DBoldDMwLTraw_1 = event->tauId_byIsolationMVArun2v1DBoldDMwLTraw_1;
    sync().byIsolationMVArun2017v2DBoldDMwLTraw2017_1 = event->tauId_byIsolationMVArun2017v2DBoldDMwLTraw2017_1;

    sync().pt_2 = event->p4_2.Pt();
    sync().pt_tau_ES_up_2 = COND_VAL(event_tau_up, (*event_tau_up)->p4_2.Pt());
    sync().pt_tau_ES_down_2 = COND_VAL(event_tau_down, (*event_tau_down)->p4_2.Pt());
    sync().phi_2 = event->p4_2.Phi();
    sync().eta_2 = event->p4_2.Eta();
    sync().m_2 = event->p4_2.mass();
    sync().q_2 = event->q_2;
    sync().d0_2 = event->dxy_2;
    sync().dZ_2 = event->dz_2;
    sync().pfmt_2 = static_cast<float>(analysis::Calculate_MT(event->p4_2, event->pfMET_p4));
    sync().iso_2 =  event->iso_2;
    sync().gen_match_2 = event->gen_match_2;

    const analysis::TauIdResults tauId_2(event->tauId_flags_2);
    sync().againstElectronLooseMVA6_2 = tauId_2.Result("againstElectronLooseMVA6");
    sync().againstElectronMediumMVA6_2 = tauId_2.Result("againstElectronMediumMVA6");
    sync().againstElectronTightMVA6_2 = tauId_2.Result("againstElectronTightMVA6");
    sync().againstElectronVLooseMVA6_2 = tauId_2.Result("againstElectronVLooseMVA6");
    sync().againstElectronVTightMVA6_2 = tauId_2.Result("againstElectronVTightMVA6");
    sync().againstMuonLoose3_2 = tauId_2.Result("againstMuonLoose3");
    sync().againstMuonTight3_2 = tauId_2.Result("againstMuonTight3");
    sync().byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = event->tauId_byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
    sync().byIsolationMVArun2v1DBoldDMwLTraw_2 = event->tauId_byIsolationMVArun2v1DBoldDMwLTraw_2;
    sync().byIsolationMVArun2017v2DBoldDMwLTraw2017_2 = event->tauId_byIsolationMVArun2017v2DBoldDMwLTraw2017_2;

    sync().pt_tt = (event->p4_1 + event->p4_2 + event->pfMET_p4).Pt();
    sync().m_vis = (event->p4_1 + event->p4_2).M();
    sync().m_sv = event->SVfit_p4.M();
    sync().m_sv_tau_ES_up = COND_VAL(event_tau_up, (*event_tau_up)->SVfit_p4.M());
    sync().m_sv_tau_ES_down = COND_VAL(event_tau_down, (*event_tau_down)->SVfit_p4.M());
    sync().m_sv_jet_ES_up = COND_VAL(event_jet_up, (*event_jet_up)->SVfit_p4.M());
    sync().m_sv_jet_ES_down = COND_VAL(event_jet_down, (*event_jet_down)->SVfit_p4.M());
    sync().mt_sv = event->SVfit_mt;

    sync().met = event->pfMET_p4.Pt();
    sync().met_tau_ES_up = COND_VAL(event_tau_up, (*event_tau_up)->pfMET_p4.Pt());
    sync().met_tau_ES_down = COND_VAL(event_tau_down, (*event_tau_down)->pfMET_p4.Pt());
    sync().met_jet_ES_up = COND_VAL(event_jet_up, event_jet_up->GetMET().GetMomentum().Pt());
    sync().met_jet_ES_down = COND_VAL(event_jet_down, event_jet_down->GetMET().GetMomentum().Pt());

    sync().mt2 = COND_VAL(event.HasBjetPair(), event.GetMT2());
    sync().mt2_tau_ES_up = COND_VAL(event_tau_up && event_tau_up->HasBjetPair(), event_tau_up->GetMT2());
    sync().mt2_tau_ES_down = COND_VAL(event_tau_down && event_tau_down->HasBjetPair(), event_tau_down->GetMT2());
    sync().mt2_jet_ES_up = COND_VAL(event_jet_up && event_jet_up->HasBjetPair(), event_jet_up->GetMT2());
    sync().mt2_jet_ES_down = COND_VAL(event_jet_down && event_jet_down->HasBjetPair(), event_jet_down->GetMT2());

    sync().metphi = static_cast<float>(TVector2::Phi_0_2pi(event->pfMET_p4.Phi()));
    sync().metphi_tau_ES_up = COND_VAL(event_tau_up, TVector2::Phi_0_2pi((*event_tau_up)->pfMET_p4.Phi()));
    sync().metphi_tau_ES_down = COND_VAL(event_tau_down, TVector2::Phi_0_2pi((*event_tau_down)->pfMET_p4.Phi()));
    sync().metphi_jet_ES_up = COND_VAL(event_jet_up,
                                       TVector2::Phi_0_2pi(event_jet_up->GetMET().GetMomentum().Phi()));
    sync().metphi_jet_ES_down = COND_VAL(event_jet_down,
                                         TVector2::Phi_0_2pi(event_jet_down->GetMET().GetMomentum().Phi()));
    sync().pzetavis = static_cast<float>(analysis::Calculate_visiblePzeta(event->p4_1, event->p4_2));

    sync().metcov00 = static_cast<float>(event->pfMET_cov[0][0]);
    sync().metcov01 = static_cast<float>(event->pfMET_cov[0][1]);
    sync().metcov10 = static_cast<float>(event->pfMET_cov[1][0]);
    sync().metcov11 = static_cast<float>(event->pfMET_cov[1][1]);

    analysis::EventInfoBase::JetCollection jets_pt20;
    analysis::EventInfoBase::JetCollection jets_pt30;

    auto select_jets = [&](analysis::EventInfoBase* event_info) {
        jets_pt20.clear();
        jets_pt30.clear();
        if(!event_info) return;

        if (run_period == analysis::Period::Run2016) {
            jets_pt20 = event_info->SelectJets(20, 4.7,
                                               analysis::JetOrdering::Pt);
            jets_pt30 = event_info->SelectJets(30, 4.7,
                                               analysis::JetOrdering::Pt);
        }

        if (run_period == analysis::Period::Run2017) {
            jets_pt20 = event_info->SelectJets(20, 4.7,
                                               analysis::JetOrdering::Pt);
            jets_pt30 = event_info->SelectJets(30, std::numeric_limits<double>::max(),
                                               analysis::JetOrdering::Pt);
        }

    };

    select_jets(&event);
    sync().mjj = COND_VAL(jets_pt20.size() >= 2, (jets_pt20.at(0).GetMomentum()
               + jets_pt20.at(1).GetMomentum()).M());
    sync().jdeta = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(0).GetMomentum().Eta()
                 - jets_pt20.at(1).GetMomentum().Eta());
    sync().jdphi = COND_VAL(jets_pt20.size() >= 2, TVector2::Phi_mpi_pi(jets_pt20.at(0).GetMomentum().Phi()
                 - jets_pt20.at(1).GetMomentum().Phi()));

    sync().njets = static_cast<int>(jets_pt30.size());
    sync().njetspt20 = static_cast<int>(jets_pt20.size());

    sync().jpt_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0).GetMomentum().Pt());
    sync().jeta_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0).GetMomentum().Eta());
    sync().jrawf_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0).GetMomentum().Phi());
    sync().jrawf_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0)->rawf());
    sync().jpt_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1).GetMomentum().Pt());
    sync().jeta_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1).GetMomentum().Eta());
    sync().jphi_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1).GetMomentum().Phi());
    sync().jrawf_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1)->rawf());

    sync().njets_vbf = static_cast<int>(jets_pt30.size());
    sync().isVBF = event.HasVBFjetPair();
    sync().jpt_vbf_1 = COND_VAL(event.HasVBFjetPair(), event.GetVBFJet(1).GetMomentum().Pt());
    sync().jeta_vbf_1 = COND_VAL(event.HasVBFjetPair(), event.GetVBFJet(1).GetMomentum().Eta());
    sync().jphi_vbf_1 = COND_VAL(event.HasVBFjetPair(), event.GetVBFJet(1).GetMomentum().Phi());
    sync().jpt_vbf_2 = COND_VAL(event.HasVBFjetPair(), event.GetVBFJet(2).GetMomentum().Pt());
    sync().jeta_vbf_2 = COND_VAL(event.HasVBFjetPair(), event.GetVBFJet(2).GetMomentum().Eta());
    sync().jphi_vbf_2 = COND_VAL(event.HasVBFjetPair(), event.GetVBFJet(2).GetMomentum().Phi());

    sync().extramuon_veto = event->extramuon_veto;
    sync().extraelec_veto = event->extraelec_veto;
    sync().nbjets = static_cast<int>(event.GetSelectedSignalJets().n_bjets);
    sync().bjet_pt_1 = COND_VAL(event.HasBjetPair(), event.GetBJet(1).GetMomentum().Pt());
    sync().bjet_eta_1 = COND_VAL(event.HasBjetPair(), event.GetBJet(1).GetMomentum().Eta());
    sync().bjet_phi_1 = COND_VAL(event.HasBjetPair(), event.GetBJet(1).GetMomentum().Phi());
    sync().bjet_rawf_1 = COND_VAL(event.HasBjetPair(), event.GetBJet(1)->rawf());
    sync().bjet_csv_1 = COND_VAL(event.HasBjetPair(), event.GetBJet(1)->csv());
    sync().bjet_deepcsv_1 = COND_VAL(event.HasBjetPair(), event.GetBJet(1)->deepcsv() >= 0
            ? event.GetBJet(1)->deepcsv() : -2);
    sync().bjet_resolution_1 = COND_VAL(event.HasBjetPair(),
                                        event.GetBJet(1)->resolution() * event.GetBJet(1).GetMomentum().E());
    sync().bjet_pt_2 = COND_VAL(event.HasBjetPair(), event.GetBJet(2).GetMomentum().Pt());
    sync().bjet_eta_2 = COND_VAL(event.HasBjetPair(), event.GetBJet(2).GetMomentum().Eta());
    sync().bjet_phi_2 = COND_VAL(event.HasBjetPair(), event.GetBJet(2).GetMomentum().Phi());
    sync().bjet_rawf_2 = COND_VAL(event.HasBjetPair(), event.GetBJet(2)->rawf());
    sync().bjet_csv_2 = COND_VAL(event.HasBjetPair(), event.GetBJet(2)->csv());
    sync().bjet_deepcsv_2 = COND_VAL(event.HasBjetPair(), event.GetBJet(2)->deepcsv() >= 0
            ? event.GetBJet(2)->deepcsv() : -2);
    sync().bjet_resolution_2 = COND_VAL(event.HasBjetPair(),
                                        event.GetBJet(2)->resolution() * event.GetBJet(2).GetMomentum().E());
    sync().ht_other_jets = event.GetHT(false, true);

    sync().kinfit_convergence = COND_VAL_INT(event.HasBjetPair() , event.GetKinFitResults().convergence);
    sync().m_kinfit = COND_VAL(event.HasBjetPair() && event.GetKinFitResults().HasValidMass(),
                               event.GetKinFitResults().mass);
    sync().m_kinfit_tau_ES_up = COND_VAL(event_tau_up && event_tau_up->HasBjetPair() &&
                                         event_tau_up->GetKinFitResults().HasValidMass(),
                                         event_tau_up->GetKinFitResults().mass);
    sync().m_kinfit_tau_ES_down = COND_VAL(event_tau_down && event_tau_down->HasBjetPair() &&
                                           event_tau_down->GetKinFitResults().HasValidMass(),
                                           event_tau_down->GetKinFitResults().mass);
    sync().m_kinfit_jet_ES_up = COND_VAL(event_jet_up && event_jet_up->HasBjetPair() &&
                                         event_jet_up->GetKinFitResults().HasValidMass(),
                                         event_jet_up->GetKinFitResults().mass);
    sync().m_kinfit_jet_ES_down = COND_VAL(event_jet_down && event_jet_down->HasBjetPair() &&
                                           event_jet_down->GetKinFitResults().HasValidMass(),
                                           event_jet_down->GetKinFitResults().mass);


    sync().mva_score_nonRes_kl1 = COND_VAL(mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_smANkin_BSMklscan", 125, 1}, &event));
    sync().mva_score_lm_320 = COND_VAL(mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_lmANkin", 320, 0}, &event));
    sync().mva_score_mm_400 = COND_VAL(mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_mmANkin", 400, 0}, &event));
    sync().mva_score_hm_650 = COND_VAL(mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_hmANkin", 650, 0}, &event));

    sync().deltaR_ll = ROOT::Math::VectorUtil::DeltaR(event->p4_1, event->p4_2);

    sync().nFatJets = static_cast<unsigned>(event.GetFatJets().size());
    const auto fatJet = event.SelectFatJet(30, 0.4);
    sync().hasFatJet = event.HasBjetPair() ? fatJet != nullptr : -1;
    sync().fatJet_pt = COND_VAL(fatJet, fatJet->GetMomentum().Pt());
    sync().fatJet_eta = COND_VAL(fatJet, fatJet->GetMomentum().Eta());
    sync().fatJet_phi = COND_VAL(fatJet, fatJet->GetMomentum().Phi());
    sync().fatJet_energy = COND_VAL(fatJet, fatJet->GetMomentum().E());
    sync().fatJet_m_softDrop = COND_VAL(fatJet, (*fatJet)->m(ntuple::TupleFatJet::MassType::SoftDrop));
    sync().fatJet_n_subjettiness_tau1 = COND_VAL(fatJet, (*fatJet)->jettiness(1));
    sync().fatJet_n_subjettiness_tau2 = COND_VAL(fatJet, (*fatJet)->jettiness(2));
    sync().fatJet_n_subjettiness_tau3 = COND_VAL(fatJet, (*fatJet)->jettiness(3));
    sync().fatJet_n_subjets = COND_VAL_INT(fatJet, (*fatJet)->subJets().size());

    sync().topWeight = static_cast<Float_t>(event->weight_top_pt);
    sync().shapeWeight = static_cast<Float_t>(event->weight_pu * event->weight_bsm_to_sm * event->weight_dy
                       * event->weight_ttbar * event->weight_wjets * event->weight_xs);
    sync().btagWeight = static_cast<Float_t>(event->weight_btag);
    sync().puweight = static_cast<Float_t>(event->weight_pu);
    sync().leptonidisoWeight = static_cast<Float_t>(event->weight_lepton_id_iso);
    sync().leptontrigWeight = static_cast<Float_t>(event->weight_lepton_trig);

    sync().lhe_n_b_partons = static_cast<int>(event->lhe_n_b_partons);
    sync().lhe_n_partons = static_cast<int>(event->lhe_n_partons);
    sync().lhe_HT = event->lhe_HT;

    sync().genJets_nTotal = event->genJets_nTotal;
    sync().genJets_nStored = static_cast<unsigned>(event->genJets_p4.size());
    sync().genJets_nStored_hadronFlavour_b = std::min<unsigned>(2, static_cast<unsigned>(
                std::count(event->genJets_hadronFlavour.begin(), event->genJets_hadronFlavour.end(), 5)));
    sync().genJets_nStored_hadronFlavour_c = static_cast<unsigned>(
                std::count(event->genJets_hadronFlavour.begin(), event->genJets_hadronFlavour.end(), 4));
    sync().jets_nTotal_hadronFlavour_b = event->jets_nTotal_hadronFlavour_b;
    sync().jets_nTotal_hadronFlavour_c = event->jets_nTotal_hadronFlavour_c;
    sync().jets_nSelected_hadronFlavour_b = 0;
    sync().jets_nSelected_hadronFlavour_c = 0;
    for(const auto& jet : event.GetJets()) {
        if(jet.GetMomentum().pt() <= 20) continue;
        if(jet->hadronFlavour() == 5) ++sync().jets_nSelected_hadronFlavour_b;
        if(jet->hadronFlavour() == 4) ++sync().jets_nSelected_hadronFlavour_c;
    }

    //mva variables
    if(event.HasBjetPair()){
        using namespace ROOT::Math::VectorUtil;
        const auto& Htt = event.GetHiggsTTMomentum(false);
        const auto& Htt_sv = event.GetHiggsTTMomentum(true);
        const auto& t1 = event.GetLeg(1).GetMomentum();
        const auto& t2 = event.GetLeg(2).GetMomentum();

        const auto& Hbb = event.GetHiggsBB().GetMomentum();
        const auto& b1 = event.GetHiggsBB().GetFirstDaughter().GetMomentum();
        const auto& b2 = event.GetHiggsBB().GetSecondDaughter().GetMomentum();

        const auto& met = event.GetMET().GetMomentum();

        sync().pt_hbb = Hbb.Pt();
        sync().pt_l1l2 = (t1+t2).Pt();
        sync().pt_htautau = (Htt+met).Pt();
        sync().pt_htautau_sv = Htt_sv.Pt();
        sync().p_zeta = analysis::Calculate_Pzeta(t1, t2,  met);
        sync().p_zetavisible = analysis::Calculate_visiblePzeta(t1, t2);
        sync().dphi_l1l2 = DeltaPhi(t1, t2);
        sync().abs_dphi_b1b2 = std::abs(DeltaPhi(b1, b2));
        sync().dphi_b1b2 = DeltaPhi(b1, b2);
        sync().dphi_l1MET = DeltaPhi(t1, met);
        sync().abs_dphi_METhtautau_sv = std::abs(DeltaPhi(Htt_sv, met));
        sync().dphi_METhtautau_sv = DeltaPhi(Htt_sv, met);
        sync().dphi_hbbMET = DeltaPhi(Hbb, met);
        sync().abs_dphi_hbbhatutau_sv = std::abs(DeltaPhi(Hbb, Htt_sv));
        sync().abs_deta_b1b2 = std::abs(b1.eta() - b2.eta());
        sync().abs_deta_l2MET = std::abs(t2.eta()-met.eta());
        sync().abs_deta_hbbMET = std::abs(Hbb.eta()-met.eta());
        sync().dR_l1l2 = DeltaR(t1, t2);
        sync().dR_hbbMET = DeltaR(Hbb, met);
        sync().dR_hbbhtautau = DeltaR(Hbb, Htt+met);
        sync().dR_l1l2Pt_htautau = DeltaR(t1, t2)*(Htt+met).Pt();
        sync().dR_l1l2Pt_htautau_sv = DeltaR(t1, t2)*Htt_sv.Pt();
        sync().MT_l1 = analysis::Calculate_MT(t1,met);
        sync().MT_htautau = analysis::Calculate_MT(Htt+met, met);
        sync().MT_htautau_sv = analysis::Calculate_MT(Htt_sv, met);
        sync().MT_tot = analysis::Calculate_TotalMT(t1, t2,met);
        sync().MT2 = event.GetMT2();
        sync().mass_top1 = analysis::four_bodies::Calculate_topPairMasses(t1, t2, b1, b2, met).first;
        sync().mass_X = analysis::four_bodies::Calculate_MX(t1, t2, b1, b2, met);
        sync().mass_H = InvariantMass(Hbb, Htt+met);
        sync().mass_H_sv = InvariantMass(Hbb, Htt_sv);
        sync().mass_H_vis = InvariantMass(Hbb, t1+t2);
        sync().mass_H_kinfit_chi2 = event.GetKinFitResults().chi2;
        sync().phi_sv = analysis::four_bodies::Calculate_phi(t1, t2, b1, b2, Htt_sv, Hbb);
        sync().phi_1_sv = analysis::four_bodies::Calculate_phi1(t1, t2, Htt_sv, Hbb);
        sync().phi_2_sv = analysis::four_bodies::Calculate_phi1(b1, b2, Htt_sv, Hbb);
        sync().costheta_METhtautau_sv = analysis::four_bodies::Calculate_cosTheta_2bodies(met, Htt_sv);
        sync().costheta_METhbb = analysis::four_bodies::Calculate_cosTheta_2bodies(met, Hbb);
        sync().costheta_b1hbb = analysis::four_bodies::Calculate_cosTheta_2bodies(b1, Hbb);
        sync().costheta_htautau_svhhMET = analysis::four_bodies::Calculate_cosTheta_2bodies(Htt_sv,
                                          event.GetResonanceMomentum(false,true));
    }

    select_jets(event_tau_up);
    sync().jpt_tau_ES_up_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0).GetMomentum().Pt());
    sync().jpt_tau_ES_up_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1).GetMomentum().Pt());
    sync().bjet_pt_tau_ES_up_1 = COND_VAL(event_tau_up && event_tau_up->HasBjetPair(),
                                          event_tau_up->GetBJet(1).GetMomentum().Pt());
    sync().bjet_pt_tau_ES_up_2 = COND_VAL(event_tau_up && event_tau_up->HasBjetPair(),
                                          event_tau_up->GetBJet(2).GetMomentum().Pt());
    sync().jpt_tau_ES_up_vbf_1 = COND_VAL(event_tau_up && event_tau_up->HasVBFjetPair(),
                                          event_tau_up->GetVBFJet(1).GetMomentum().Pt());
    sync().jpt_tau_ES_up_vbf_2 = COND_VAL(event_tau_up && event_tau_up->HasVBFjetPair(),
                                          event_tau_up->GetVBFJet(2).GetMomentum().Pt());
    sync().mva_score_nonRes_kl1_tau_ES_up = COND_VAL(event_tau_up && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_smANkin_BSMklscan", 125, 1}, event_tau_up));
    sync().mva_score_lm_320_tau_ES_up = COND_VAL(event_tau_up && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_lmANkin", 320, 0}, event_tau_up));
    sync().mva_score_mm_400_tau_ES_up = COND_VAL(event_tau_up && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_mmANkin", 400, 0}, event_tau_up));
    sync().mva_score_hm_650_tau_ES_up = COND_VAL(event_tau_up && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_hmANkin", 650, 0}, event_tau_up));

    select_jets(event_tau_down);
    sync().jpt_tau_ES_down_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0).GetMomentum().Pt());
    sync().jpt_tau_ES_down_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1).GetMomentum().Pt());
    sync().bjet_pt_tau_ES_down_1 = COND_VAL(event_tau_down && event_tau_down->HasBjetPair(),
                                            event_tau_down->GetBJet(1).GetMomentum().Pt());
    sync().bjet_pt_tau_ES_down_2 = COND_VAL(event_tau_down && event_tau_down->HasBjetPair(),
                                            event_tau_down->GetBJet(2).GetMomentum().Pt());
    sync().jpt_tau_ES_down_vbf_1 = COND_VAL(event_tau_down && event_tau_down->HasVBFjetPair(),
                                            event_tau_down->GetVBFJet(1).GetMomentum().Pt());
    sync().jpt_tau_ES_down_vbf_2 = COND_VAL(event_tau_down && event_tau_down->HasVBFjetPair(),
                                            event_tau_down->GetVBFJet(2).GetMomentum().Pt());
    sync().mva_score_nonRes_kl1_tau_ES_down = COND_VAL(event_tau_down && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_smANkin_BSMklscan", 125, 1}, event_tau_down));
    sync().mva_score_lm_320_tau_ES_down = COND_VAL(event_tau_down && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_lmANkin", 320, 0}, event_tau_down));
    sync().mva_score_mm_400_tau_ES_down = COND_VAL(event_tau_down && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_mmANkin", 400, 0}, event_tau_down));
    sync().mva_score_hm_650_tau_ES_down = COND_VAL(event_tau_down && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_hmANkin", 650, 0}, event_tau_down));

    select_jets(event_jet_up);
    sync().jpt_jet_ES_up_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0).GetMomentum().Pt());
    sync().jpt_jet_ES_up_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1).GetMomentum().Pt());
    sync().bjet_pt_jet_ES_up_1 = COND_VAL(event_jet_up && event_jet_up->HasBjetPair(),
                                          event_jet_up->GetBJet(1).GetMomentum().Pt());
    sync().bjet_pt_jet_ES_up_2 = COND_VAL(event_jet_up && event_jet_up->HasBjetPair(),
                                          event_jet_up->GetBJet(2).GetMomentum().Pt());
    sync().jpt_jet_ES_up_vbf_1 = COND_VAL(event_jet_up && event_jet_up->HasVBFjetPair(),
                                          event_jet_up->GetVBFJet(1).GetMomentum().Pt());
    sync().jpt_jet_ES_up_vbf_2 = COND_VAL(event_jet_up && event_jet_up->HasVBFjetPair(),
                                          event_jet_up->GetVBFJet(2).GetMomentum().Pt());
    sync().mva_score_nonRes_kl1_jet_ES_up = COND_VAL(event_jet_up && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_smANkin_BSMklscan", 125, 1}, event_jet_up));
    sync().mva_score_lm_320_jet_ES_up = COND_VAL(event_jet_up && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_lmANkin", 320, 0}, event_jet_up));
    sync().mva_score_mm_400_jet_ES_up = COND_VAL(event_jet_up && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_mmANkin", 400, 0}, event_jet_up));
    sync().mva_score_hm_650_jet_ES_up = COND_VAL(event_jet_up && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_hmANkin", 650, 0}, event_jet_up));

    select_jets(event_jet_down);
    sync().jpt_jet_ES_down_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0).GetMomentum().Pt());
    sync().jpt_jet_ES_down_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1).GetMomentum().Pt());
    sync().bjet_pt_jet_ES_down_1 = COND_VAL(event_jet_down && event_jet_down->HasBjetPair(),
                                            event_jet_down->GetBJet(1).GetMomentum().Pt());
    sync().bjet_pt_jet_ES_down_2 = COND_VAL(event_jet_down && event_jet_down->HasBjetPair(),
                                            event_jet_down->GetBJet(2).GetMomentum().Pt());
    sync().jpt_jet_ES_down_vbf_1 = COND_VAL(event_jet_down && event_jet_down->HasVBFjetPair(),
                                            event_jet_down->GetVBFJet(1).GetMomentum().Pt());
    sync().jpt_jet_ES_down_vbf_2 = COND_VAL(event_jet_down && event_jet_down->HasVBFjetPair(),
                                            event_jet_down->GetVBFJet(2).GetMomentum().Pt());
    sync().mva_score_nonRes_kl1_jet_ES_down = COND_VAL(event_jet_down && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_smANkin_BSMklscan", 125, 1}, event_jet_down));
    sync().mva_score_lm_320_jet_ES_down = COND_VAL(event_jet_down && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_lmANkin", 320, 0}, event_jet_down));
    sync().mva_score_mm_400_jet_ES_down = COND_VAL(event_jet_down && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_mmANkin", 400, 0}, event_jet_down));
    sync().mva_score_hm_650_jet_ES_down = COND_VAL(event_jet_down && mva_reader, mva_reader->Evaluate(analysis::mva_study::MvaReader::MvaKey{"mva_hmANkin", 650, 0}, event_jet_down));

    sync.Fill();
}
}

#undef COND_VAL
#undef COND_VAL_INT
