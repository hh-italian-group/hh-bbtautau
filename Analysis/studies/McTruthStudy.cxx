/*! Study of base analysis object at level of MC truth.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/LightBaseBigTreeAnalyzer.h"

class McTruthStudyData : public analysis::BaseAnalyzerData {
public:
    McTruthStudyData(std::shared_ptr<TFile> outputFile) : BaseAnalyzerData(outputFile) {}

    TH1D_ENTRY(MC_visible_pt, 20, 0, 200)
    TH1D_ENTRY(MC_visible_eta, 20, 0, 5)
    TH1D_ENTRY(MC_visible_phi, 32, 0, 3.2)
    TH1D_ENTRY(MC_visible_energy, 50, 7000, 8000)

    TH1D_ENTRY(MC_invisible_pt, 20, 0, 200)
    TH1D_ENTRY(MC_invisible_eta, 20, 0, 5)
    TH1D_ENTRY(MC_invisible_phi, 32, 0, 3.2)
    TH1D_ENTRY(MC_invisible_energy, 40, 0, 400)

    TH1D_ENTRY(MC_underlying_invisible_pt, 20, 0, 20)
    TH1D_ENTRY(MC_underlying_invisible_eta, 20, 0, 5)
    TH1D_ENTRY(MC_underlying_invisible_phi, 32, 0, 3.2)
    TH1D_ENTRY(MC_underlying_invisible_energy, 20, 0, 20)

    TH2D_ENTRY(MVA_MET_inaccuracy_pt, 20, 0, 200, 20, -1, 1)
    TH2D_ENTRY(MVA_MET_inaccuracy_pt_diff, 20, 0, 200, 10, 0, 1)
    TH1D_ENTRY(MVA_MET_inaccuracy_phi, 32, 0, 3.2)

    TH2D_ENTRY(MSSM_H_inaccuracy_pt, 20, 0, 200, 20, -1, 1)
    TH2D_ENTRY(MSSM_H_inaccuracy_pt_diff, 20, 0, 200, 10, 0, 1)
    TH1D_ENTRY(MSSM_H_inaccuracy_phi, 32, 0, 3.2)

    TH2D_ENTRY(MSSM_H_MET_inaccuracy_pt, 20, 0, 200, 20, -1, 1)
    TH2D_ENTRY(MSSM_H_MET_inaccuracy_pt_diff, 20, 0, 200, 10, 0, 1)
    TH1D_ENTRY(MSSM_H_MET_inaccuracy_phi, 32, 0, 3.2)

    TH2D_ENTRY(MSSM_H_MET_obs_inaccuracy_pt, 20, 0, 200, 20, -1, 1)
    TH2D_ENTRY(MSSM_H_MET_obs_inaccuracy_pt_diff, 20, 0, 200, 10, 0, 1)
    TH1D_ENTRY(MSSM_H_MET_obs_inaccuracy_phi, 32, 0, 3.2)

    TH2D_ENTRY(MSSM_H_SV_inaccuracy_pt, 20, 0, 200, 20, -1, 1)
    TH2D_ENTRY(MSSM_H_SV_inaccuracy_pt_diff, 20, 0, 200, 10, 0, 1)
    TH1D_ENTRY(MSSM_H_SV_inaccuracy_phi, 32, 0, 3.2)

    TH2D_ENTRY(MSSM_H_SV_obs_inaccuracy_pt, 20, 0, 200, 20, -1, 1)
    TH2D_ENTRY(MSSM_H_SV_obs_inaccuracy_pt_diff, 20, 0, 200, 10, 0, 1)
    TH1D_ENTRY(MSSM_H_SV_obs_inaccuracy_phi, 32, 0, 3.2)

    TH1D_ENTRY(MC_invisible_tau_ratio, 11, 0, 1.1)

    TH1D_ENTRY(DeltaRmin1_visible, 600, 0, 3)
    TH1D_ENTRY(DeltaRmin2_visible, 600, 0, 3)
    TH1D_ENTRY(DeltaRmin1_original, 600, 0, 3)
    TH1D_ENTRY(DeltaRmin2_original, 600, 0, 3)
    TH1D_ENTRY(DeltaRbjets_MC, 600, 0, 6)
    TH1D_ENTRY(deltaRmin_MC, 600, 0, 6)
    TH1D_ENTRY(deltaRmax_visible, 600, 0, 6)
    TH1D_ENTRY(deltaRmax_original, 600, 0, 6)
    TH1D_ENTRY(deltaPtMax, 250, 0, 5)
    TH1D_ENTRY(deltaPtMax_vis, 250, 0, 5)
    TH1D_ENTRY(MinPtBjetsMC, 20, 0, 200)
    TH1D_ENTRY(MassBB_MC, 30, 0, 300)
    TH1D_ENTRY(MassBB_MCvis, 30, 0, 300)
    TH1D_ENTRY(goodElectronsFromZee, 5, -0.5, 4.5)
    TH1D_ENTRY(hardElectronsZee, 5, -0.5, 4.5)
    TH1D_ENTRY(hardMuonsZee, 5, -0.5, 4.5)
    TH1D_ENTRY(goodMuonsFromZmm, 5, -0.5, 4.5)
    TH1D_ENTRY(hardElectronsZmm, 5, -0.5, 4.5)
    TH1D_ENTRY(hardMuonsZmm, 5, -0.5, 4.5)
};


class McTruthStudy : public analysis::LightBaseBigTreeAnalyzer {
public:
    McTruthStudy(const std::string& channelName, const std::string& inputFileName, const std::string& outputFileName,
                 const std::string& configFileName, const std::string& _prefix = "none", size_t _maxNumberOfEvents = 0)
        : BaseFlatTreeProducer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents),
          LightBaseBigTreeAnalyzer(channelName, inputFileName, outputFileName, configFileName, _prefix,
                                   _maxNumberOfEvents),
          anaData(*outputFile)
    {
    }

protected:
    virtual void AnalyzeSelection(const analysis::SelectionResults& selection) override
    {
        using namespace analysis;

        if(selection.GetFinalStateMC().b_jets.size() < 2 || selection.GetFinalStateMC().taus.size() < 2
                || selection.bjets_all.size() < 2) return;

        const GenParticle* MSSM_H_MC = selection.GetFinalStateMC().resonance->GetHardInteractionOrigin();
        const GenParticle* H_TT_MC = selection.GetFinalStateMC().Higgs_TauTau->GetHardInteractionOrigin();
//        const GenParticle* H_BB_MC = selection.GetFinalStateMC().Higgs_BB->GetHardInteractionOrigin();

        const VisibleGenObject& tau_1_visible = selection.GetFinalStateMC().taus.at(0);
        const VisibleGenObject& tau_2_visible = selection.GetFinalStateMC().taus.at(1);
        const GenParticle* tau_1_MC = tau_1_visible.origin->GetHardInteractionOrigin();
        const GenParticle* tau_2_MC = tau_2_visible.origin->GetHardInteractionOrigin();

        const VisibleGenObject& bjet_1_visible = selection.GetFinalStateMC().b_jets.at(0);
        const VisibleGenObject& bjet_2_visible = selection.GetFinalStateMC().b_jets.at(1);
        const GenParticle* bjet_1_MC = bjet_1_visible.origin->GetHardInteractionOrigin();
        const GenParticle* bjet_2_MC = bjet_2_visible.origin->GetHardInteractionOrigin();

        const TLorentzVector& tau_1_RECO = selection.GetLeg1().momentum;
        const TLorentzVector& tau_2_RECO = selection.GetLeg2().momentum;
        const TLorentzVector& bjet_1_RECO = selection.bjets_all.at(0).momentum;
        const TLorentzVector& bjet_2_RECO = selection.bjets_all.at(1).momentum;
        const TLorentzVector& MVA_MET = MakeLorentzVectorPtEtaPhiM(selection.MET_with_recoil_corrections.pt, 0,
                                                                  selection.MET_with_recoil_corrections.phi, 0);

        unsigned n_lept_bjet = 0;
        for(const VisibleGenObject& bjet : selection.GetFinalStateMC().b_jets) {
            if(bjet.finalStateChargedLeptons.size() != 0)
                ++n_lept_bjet;
        }

        TLorentzVector MC_visible, MC_invisible;
        for(const GenParticle& genParticle : genEvent.genParticles) {
            if(genParticle.status != particles::FinalStateParticle) continue;
            if(particles::neutrinos.count(genParticle.pdg.Code))
                MC_invisible += genParticle.momentum;
            else
                MC_visible += genParticle.momentum;
        }

        anaData.MC_visible_pt().Fill(MC_visible.Pt());
        anaData.MC_visible_eta().Fill(std::abs(MC_visible.Eta()));
        anaData.MC_visible_phi().Fill(std::abs(MC_visible.Phi()));
        anaData.MC_visible_energy().Fill(MC_visible.E());

        anaData.MC_invisible_pt().Fill(MC_invisible.Pt());
        anaData.MC_invisible_eta().Fill(std::abs(MC_invisible.Eta()));
        anaData.MC_invisible_phi().Fill(std::abs(MC_invisible.Phi()));
        anaData.MC_invisible_energy().Fill(MC_invisible.E());

        std::ostringstream ss_name;
        ss_name << "n_lept_b_" << n_lept_bjet;
        const std::string inc_name = ss_name.str();
        if(MC_invisible.Pt() < 30)
            ss_name << "_low_inv";
        else
            ss_name << "_high_inv";
        const std::string cat_name = ss_name.str();


        const VisibleGenObject MSSM_H_genVisible(MSSM_H_MC);
        const GenParticle* genProton = *genEvent.primaryParticles.begin();
        const VisibleGenObject underlying_visible(genProton, MSSM_H_genVisible.particlesProcessed);
        anaData.MC_underlying_invisible_pt().Fill(underlying_visible.invisibleMomentum.Pt());
        anaData.MC_underlying_invisible_eta().Fill(std::abs(underlying_visible.invisibleMomentum.Eta()));
        anaData.MC_underlying_invisible_phi().Fill(std::abs(underlying_visible.invisibleMomentum.Phi()));
        anaData.MC_underlying_invisible_energy().Fill(underlying_visible.invisibleMomentum.E());

        const TLorentzVector MVA_MET_inaccuracy = MC_invisible - MVA_MET;
        const double MVA_MET_inaccuracy_pt = (MC_invisible.Pt() - MVA_MET.Pt()) / MC_invisible.Pt();
        const double MVA_MET_inaccuracy_pt_diff = MVA_MET_inaccuracy.Pt() / MC_invisible.Pt();
        anaData.MVA_MET_inaccuracy_pt().Fill(MC_invisible.Pt(), MVA_MET_inaccuracy_pt);
        anaData.MVA_MET_inaccuracy_pt(inc_name).Fill(MC_invisible.Pt(), MVA_MET_inaccuracy_pt);
        anaData.MVA_MET_inaccuracy_pt_diff().Fill(MC_invisible.Pt(), MVA_MET_inaccuracy_pt_diff);
        anaData.MVA_MET_inaccuracy_pt_diff(inc_name).Fill(MC_invisible.Pt(), MVA_MET_inaccuracy_pt_diff);
        anaData.MVA_MET_inaccuracy_phi().Fill(std::abs(MC_invisible.DeltaPhi(MVA_MET)));
        anaData.MVA_MET_inaccuracy_phi(inc_name).Fill(std::abs(MC_invisible.DeltaPhi(MVA_MET)));

        const TLorentzVector MSSM_H_MC_visible = MSSM_H_MC->momentum - MSSM_H_genVisible.invisibleMomentum;
        const TLorentzVector MSSM_H_RECO = tau_1_RECO + tau_2_RECO + bjet_1_RECO + bjet_2_RECO;
        const TLorentzVector MSSM_H_inaccuracy = MSSM_H_MC_visible - MSSM_H_RECO;
        const double MSSM_H_inaccuracy_pt = (MSSM_H_MC_visible.Pt() - MSSM_H_RECO.Pt()) / MSSM_H_MC_visible.Pt();
        const double MSSM_H_inaccuracy_pt_diff = MSSM_H_inaccuracy.Pt() / MSSM_H_MC_visible.Pt();
        anaData.MSSM_H_inaccuracy_pt().Fill(MSSM_H_MC_visible.Pt(), MSSM_H_inaccuracy_pt);
        anaData.MSSM_H_inaccuracy_pt(inc_name).Fill(MSSM_H_MC_visible.Pt(), MSSM_H_inaccuracy_pt);
        anaData.MSSM_H_inaccuracy_pt(cat_name).Fill(MSSM_H_MC_visible.Pt(), MSSM_H_inaccuracy_pt);
        anaData.MSSM_H_inaccuracy_pt_diff().Fill(MSSM_H_MC_visible.Pt(), MSSM_H_inaccuracy_pt_diff);
        anaData.MSSM_H_inaccuracy_pt_diff(inc_name).Fill(MSSM_H_MC_visible.Pt(), MSSM_H_inaccuracy_pt_diff);
        anaData.MSSM_H_inaccuracy_pt_diff(cat_name).Fill(MSSM_H_MC_visible.Pt(), MSSM_H_inaccuracy_pt_diff);
        anaData.MSSM_H_inaccuracy_phi().Fill(std::abs(MSSM_H_MC_visible.DeltaPhi(MSSM_H_RECO)));
        anaData.MSSM_H_inaccuracy_phi(inc_name).Fill(std::abs(MSSM_H_MC_visible.DeltaPhi(MSSM_H_RECO)));
        anaData.MSSM_H_inaccuracy_phi(cat_name).Fill(std::abs(MSSM_H_MC_visible.DeltaPhi(MSSM_H_RECO)));

        const TLorentzVector MSSM_H_MET_RECO = MSSM_H_RECO + MVA_MET;
        const TLorentzVector MSSM_H_MET_inaccuracy = MSSM_H_MC->momentum - MSSM_H_MET_RECO;
        const double MSSM_H_MET_inaccuracy_pt = (MSSM_H_MC->momentum.Pt() - MSSM_H_MET_RECO.Pt())
                / MSSM_H_MC->momentum.Pt();
        const double MSSM_H_MET_inaccuracy_pt_diff = MSSM_H_MET_inaccuracy.Pt() / MSSM_H_MC->momentum.Pt();
        anaData.MSSM_H_MET_inaccuracy_pt().Fill(MSSM_H_MC->momentum.Pt(), MSSM_H_MET_inaccuracy_pt);
        anaData.MSSM_H_MET_inaccuracy_pt(inc_name).Fill(MSSM_H_MC->momentum.Pt(), MSSM_H_MET_inaccuracy_pt);
        anaData.MSSM_H_MET_inaccuracy_pt(cat_name).Fill(MSSM_H_MC->momentum.Pt(), MSSM_H_MET_inaccuracy_pt);
        anaData.MSSM_H_MET_inaccuracy_pt_diff().Fill(MSSM_H_MC->momentum.Pt(), MSSM_H_MET_inaccuracy_pt_diff);
        anaData.MSSM_H_MET_inaccuracy_pt_diff(inc_name).Fill(MSSM_H_MC->momentum.Pt(), MSSM_H_MET_inaccuracy_pt_diff);
        anaData.MSSM_H_MET_inaccuracy_pt_diff(cat_name).Fill(MSSM_H_MC->momentum.Pt(), MSSM_H_MET_inaccuracy_pt_diff);
        anaData.MSSM_H_MET_inaccuracy_phi().Fill(std::abs(MSSM_H_MC->momentum.DeltaPhi(MSSM_H_MET_RECO)));
        anaData.MSSM_H_MET_inaccuracy_phi(inc_name).Fill(std::abs(MSSM_H_MC->momentum.DeltaPhi(MSSM_H_MET_RECO)));
        anaData.MSSM_H_MET_inaccuracy_phi(cat_name).Fill(std::abs(MSSM_H_MC->momentum.DeltaPhi(MSSM_H_MET_RECO)));

        const sv_fit::FitResults sv_fit_result = sv_fit::Fit(sv_fit::FitAlgorithm::MarkovChain, selection.higgs,
                                                             selection.MET_with_recoil_corrections, 1);
        if(!sv_fit_result.has_valid_momentum)
            throw analysis::exception("svFit not converged");
        const TLorentzVector MSSM_H_SV = sv_fit_result.momentum + bjet_1_RECO + bjet_2_RECO;
        const TLorentzVector MSSM_H_SV_inaccuracy = MSSM_H_MC->momentum - MSSM_H_SV;
        const double MSSM_H_SV_inaccuracy_pt = (MSSM_H_MC->momentum.Pt() - MSSM_H_SV.Pt()) / MSSM_H_MC->momentum.Pt();
        const double MSSM_H_SV_inaccuracy_pt_diff = MSSM_H_SV_inaccuracy.Pt() / MSSM_H_MC->momentum.Pt();
        anaData.MSSM_H_SV_inaccuracy_pt().Fill(MSSM_H_MC->momentum.Pt(), MSSM_H_SV_inaccuracy_pt);
        anaData.MSSM_H_SV_inaccuracy_pt(inc_name).Fill(MSSM_H_MC->momentum.Pt(), MSSM_H_SV_inaccuracy_pt);
        anaData.MSSM_H_SV_inaccuracy_pt(cat_name).Fill(MSSM_H_MC->momentum.Pt(), MSSM_H_SV_inaccuracy_pt);
        anaData.MSSM_H_SV_inaccuracy_pt_diff().Fill(MSSM_H_MC->momentum.Pt(), MSSM_H_SV_inaccuracy_pt_diff);
        anaData.MSSM_H_SV_inaccuracy_pt_diff(inc_name).Fill(MSSM_H_MC->momentum.Pt(), MSSM_H_SV_inaccuracy_pt_diff);
        anaData.MSSM_H_SV_inaccuracy_pt_diff(cat_name).Fill(MSSM_H_MC->momentum.Pt(), MSSM_H_SV_inaccuracy_pt_diff);
        anaData.MSSM_H_SV_inaccuracy_phi().Fill(std::abs(MSSM_H_MC->momentum.DeltaPhi(MSSM_H_SV)));
        anaData.MSSM_H_SV_inaccuracy_phi(inc_name).Fill(std::abs(MSSM_H_MC->momentum.DeltaPhi(MSSM_H_SV)));
        anaData.MSSM_H_SV_inaccuracy_phi(cat_name).Fill(std::abs(MSSM_H_MC->momentum.DeltaPhi(MSSM_H_SV)));

        const VisibleGenObject H_TT_visible(H_TT_MC);
        const double MC_invisible_tau_ratio = H_TT_visible.invisibleMomentum.Pt() / MC_invisible.Pt();
        anaData.MC_invisible_tau_ratio().Fill(MC_invisible_tau_ratio);
        anaData.MC_invisible_tau_ratio(inc_name).Fill(MC_invisible_tau_ratio);
        anaData.MC_invisible_tau_ratio(cat_name).Fill(MC_invisible_tau_ratio);

        const TLorentzVector Observable_MC = tau_1_visible.origin->momentum
                + tau_2_visible.origin->momentum + bjet_1_visible.origin->momentum + bjet_2_visible.origin->momentum;
        const TLorentzVector MSSM_H_MET_obs_inaccuracy = Observable_MC - MSSM_H_MET_RECO;
        const double MSSM_H_MET_obs_inaccuracy_pt = (Observable_MC.Pt() - MSSM_H_MET_RECO.Pt()) / Observable_MC.Pt();
        const double MSSM_H_MET_obs_inaccuracy_pt_diff = MSSM_H_MET_obs_inaccuracy.Pt() / Observable_MC.Pt();
        anaData.MSSM_H_MET_obs_inaccuracy_pt().Fill(Observable_MC.Pt(), MSSM_H_MET_obs_inaccuracy_pt);
        anaData.MSSM_H_MET_obs_inaccuracy_pt(inc_name).Fill(Observable_MC.Pt(), MSSM_H_MET_obs_inaccuracy_pt);
        anaData.MSSM_H_MET_obs_inaccuracy_pt(cat_name).Fill(Observable_MC.Pt(), MSSM_H_MET_obs_inaccuracy_pt);
        anaData.MSSM_H_MET_obs_inaccuracy_pt_diff().Fill(Observable_MC.Pt(), MSSM_H_MET_obs_inaccuracy_pt_diff);
        anaData.MSSM_H_MET_obs_inaccuracy_pt_diff(inc_name).Fill(Observable_MC.Pt(), MSSM_H_MET_obs_inaccuracy_pt_diff);
        anaData.MSSM_H_MET_obs_inaccuracy_pt_diff(cat_name).Fill(Observable_MC.Pt(), MSSM_H_MET_obs_inaccuracy_pt_diff);
        anaData.MSSM_H_MET_obs_inaccuracy_phi().Fill(std::abs(Observable_MC.DeltaPhi(MSSM_H_MET_RECO)));
        anaData.MSSM_H_MET_obs_inaccuracy_phi(inc_name).Fill(std::abs(Observable_MC.DeltaPhi(MSSM_H_MET_RECO)));
        anaData.MSSM_H_MET_obs_inaccuracy_phi(cat_name).Fill(std::abs(Observable_MC.DeltaPhi(MSSM_H_MET_RECO)));

        const TLorentzVector MSSM_H_SV_obs_inaccuracy = Observable_MC - MSSM_H_SV;
        const double MSSM_H_SV_obs_inaccuracy_pt = (Observable_MC.Pt() - MSSM_H_SV.Pt()) / Observable_MC.Pt();
        const double MSSM_H_SV_obs_inaccuracy_pt_diff = MSSM_H_SV_obs_inaccuracy.Pt() / Observable_MC.Pt();
        anaData.MSSM_H_SV_obs_inaccuracy_pt().Fill(Observable_MC.Pt(), MSSM_H_SV_obs_inaccuracy_pt);
        anaData.MSSM_H_SV_obs_inaccuracy_pt(inc_name).Fill(Observable_MC.Pt(), MSSM_H_SV_obs_inaccuracy_pt);
        anaData.MSSM_H_SV_obs_inaccuracy_pt(cat_name).Fill(Observable_MC.Pt(), MSSM_H_SV_obs_inaccuracy_pt);
        anaData.MSSM_H_SV_obs_inaccuracy_pt_diff().Fill(Observable_MC.Pt(), MSSM_H_SV_obs_inaccuracy_pt_diff);
        anaData.MSSM_H_SV_obs_inaccuracy_pt_diff(inc_name).Fill(Observable_MC.Pt(), MSSM_H_SV_obs_inaccuracy_pt_diff);
        anaData.MSSM_H_SV_obs_inaccuracy_pt_diff(cat_name).Fill(Observable_MC.Pt(), MSSM_H_SV_obs_inaccuracy_pt_diff);
        anaData.MSSM_H_SV_obs_inaccuracy_phi().Fill(std::abs(Observable_MC.DeltaPhi(MSSM_H_SV)));
        anaData.MSSM_H_SV_obs_inaccuracy_phi(inc_name).Fill(std::abs(Observable_MC.DeltaPhi(MSSM_H_SV)));
        anaData.MSSM_H_SV_obs_inaccuracy_phi(cat_name).Fill(std::abs(Observable_MC.DeltaPhi(MSSM_H_SV)));

        if (bjet_1_MC->momentum.Pt() > 20 && bjet_2_MC->momentum.Pt() > 20 &&
            std::abs(bjet_1_MC->momentum.Eta()) < 2.4 && std::abs(bjet_2_MC->momentum.Eta()) < 2.4 ) {

            const double deltaRmin_firstCouple = std::min(bjet_1_MC->momentum.DeltaR(tau_1_MC->momentum),
                                                          bjet_1_MC->momentum.DeltaR(tau_2_MC->momentum));
            const double deltaRmin_secondCouple = std::min(bjet_2_MC->momentum.DeltaR(tau_1_MC->momentum),
                                                           bjet_2_MC->momentum.DeltaR(tau_2_MC->momentum));
            const double deltaRmin_MC = std::min(deltaRmin_firstCouple, deltaRmin_secondCouple);

            anaData.deltaRmin_MC().Fill(deltaRmin_MC);
            anaData.DeltaRbjets_MC().Fill(bjet_1_MC->momentum.DeltaR(bjet_2_MC->momentum));
            anaData.MinPtBjetsMC().Fill(std::min(bjet_1_MC->momentum.Pt(), bjet_2_MC->momentum.Pt()));
            const TLorentzVector bb_MC = bjet_1_MC->momentum + bjet_2_MC->momentum;
            const TLorentzVector bb_MC_visible = bjet_1_visible.visibleMomentum + bjet_2_visible.visibleMomentum;

            anaData.MassBB_MC().Fill(bb_MC.M());
            anaData.MassBB_MCvis().Fill(bb_MC_visible.M());

            double deltaRmin1_original = std::numeric_limits<double>::max();
            double deltaRmin2_original = std::numeric_limits<double>::max();
            double deltaRmin1_visible = std::numeric_limits<double>::max();
            double deltaRmin2_visible = std::numeric_limits<double>::max();
            unsigned index_bjet1 = 0;
            unsigned index_bjet2 = 0;
            unsigned index_bjet1_vis = 0;
            unsigned index_bjet2_vis = 0;
            for (unsigned i = 0; i < selection.bjets_all.size(); ++i){
                const Candidate& bjet = selection.bjets_all.at(i);
//                double deltaPt1 = std::abs(bjet.momentum.Pt() - bjet1_MC->momentum.Pt())/bjet1_MC->momentum.Pt();
//                double deltaPt2 = std::abs(bjet.momentum.Pt() - bjet2_MC->momentum.Pt())/bjet2_MC->momentum.Pt();
                if (bjet.momentum.DeltaR(bjet_1_MC->momentum) < deltaRmin1_original /*&& deltaPt1 < 0.4*/){
                    deltaRmin1_original = bjet.momentum.DeltaR(bjet_1_MC->momentum);
                    index_bjet1 = i;
                }

                if (bjet.momentum.DeltaR(bjet_2_MC->momentum) < deltaRmin2_original /*&& deltaPt2 < 0.4*/){
                    deltaRmin2_original = bjet.momentum.DeltaR(bjet_2_MC->momentum);
                    index_bjet2 = i;
                }

                if (bjet.momentum.DeltaR(bjet_1_visible.visibleMomentum) < deltaRmin1_visible){
                    deltaRmin1_visible = bjet.momentum.DeltaR(bjet_1_visible.visibleMomentum);
                    index_bjet1_vis = i;
                }

                if (bjet.momentum.DeltaR(bjet_2_visible.visibleMomentum) < deltaRmin2_visible){
                    deltaRmin2_visible = bjet.momentum.DeltaR(bjet_2_visible.visibleMomentum);
                    index_bjet2_vis = i;
                }
            }

            anaData.DeltaRmin1_original().Fill(deltaRmin1_original);
            anaData.DeltaRmin2_original().Fill(deltaRmin2_original);

            double deltaRmax_original = std::max(deltaRmin1_original,deltaRmin2_original);
            anaData.deltaRmax_original().Fill(deltaRmax_original);

            anaData.DeltaRmin1_visible().Fill(deltaRmin1_visible);
            anaData.DeltaRmin2_visible().Fill(deltaRmin2_visible);

            double deltaRmax_visible = std::max(deltaRmin1_visible,deltaRmin2_visible);
            anaData.deltaRmax_visible().Fill(deltaRmax_visible);

            double deltaPtMax;
            if (deltaRmax_original == std::numeric_limits<double>::max()){
                deltaPtMax = std::numeric_limits<double>::max();
            }
            else {
                const Candidate& selectedBjets1 = selection.bjets_all.at(index_bjet1);
                const Candidate& selectedBjets2 = selection.bjets_all.at(index_bjet2);
                double deltaPt1 = std::abs(selectedBjets1.momentum.Pt() - bjet_1_MC->momentum.Pt());
                double deltaPt2 = std::abs(selectedBjets2.momentum.Pt() - bjet_2_MC->momentum.Pt());
                deltaPtMax = std::max(deltaPt1/bjet_1_MC->momentum.Pt(), deltaPt2/bjet_2_MC->momentum.Pt());
            }
            anaData.deltaPtMax().Fill(deltaPtMax);

            double deltaPtMax_vis;
            if (deltaRmax_visible == std::numeric_limits<double>::max()){
                deltaPtMax_vis = std::numeric_limits<double>::max();
            }
            else {
                const Candidate& selectedBjets1 = selection.bjets_all.at(index_bjet1_vis);
                const Candidate& selectedBjets2 = selection.bjets_all.at(index_bjet2_vis);
                double deltaPt1 = std::abs(selectedBjets1.momentum.Pt() - bjet_1_visible.visibleMomentum.Pt());
                double deltaPt2 = std::abs(selectedBjets2.momentum.Pt() - bjet_2_visible.visibleMomentum.Pt());
                deltaPtMax_vis =
                        std::max(deltaPt1/bjet_1_visible.visibleMomentum.Pt(),
                                 deltaPt2/bjet_2_visible.visibleMomentum.Pt());
            }
            anaData.deltaPtMax_vis().Fill(deltaPtMax_vis);
        }
    }

private:
    McTruthStudyData anaData;
};
