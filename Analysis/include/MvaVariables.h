/*! Definition of MvaVariables.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <random>
#include "h-tautau/Analysis/include/EventTuple.h"
#include "hh-bbtautau/Analysis/include/MT2.h"

namespace analysis {

#define VAR(name, formula) if(IsEnabled(name)) SetValue(name, formula)
class MvaVariables {
public:
    using VarNameSet = std::set<std::string>;

    MvaVariables(bool _split_training_testing, uint_fast32_t seed, const VarNameSet& _enabled_vars,
                 const VarNameSet& _disabled_vars) :
        split_training_testing(_split_training_testing), gen(seed), test_vs_training(0, 1),
        enabled_vars(_enabled_vars), disabled_vars(_disabled_vars)
    {
    }

    virtual ~MvaVariables() {}
    virtual void SetValue(const std::string& name, double value) = 0;
    virtual void AddEventVariables(bool istraining, int mass, std::string tree, double weight) = 0;

    bool IsEnabled(const std::string& name) const
    {
        return (enabled_vars.size() && enabled_vars.count(name)) || !disabled_vars.count(name);
    }

    void AddEvent(std::string tree, const ntuple::Event& event, int mass = 1, double sample_weight = 1.)
    {
        auto bb= event.jets_p4[0] + event.jets_p4[1];
        auto leptons= event.p4_1 + event.p4_2;
        auto leptonsMET= event.p4_1 + event.p4_2 + event.pfMET_p4;

        VAR("pt_b1", event.jets_p4[0].pt());
        VAR("pt_b2", event.jets_p4[1].pt());
        VAR("pt_hbb", bb.pt());
        VAR("pt_l1", event.p4_1.pt());
        VAR("pt_l2", event.p4_2.pt());
        VAR("pt_l1l2", leptons.pt());
        VAR("pt_l1l2met", leptonsMET.pt());
        VAR("pt_htautau", event.SVfit_p4.pt());
        VAR("pt_met", event.pfMET_p4.pt());
        VAR("HT_otherjets", Calculate_HT(event.jets_p4.begin()+2, event.jets_p4.end()));
        VAR("p_zeta", Calculate_Pzeta(event.p4_1, event.p4_2, event.pfMET_p4));
        VAR("p_zetavisible", Calculate_visiblePzeta(event.p4_1,event.p4_2));
        VAR("abs_dphi_l1l2", std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.p4_2)));
        VAR("dphi_l1l2", ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.p4_2));
        VAR("abs_dphi_b1b2", std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.jets_p4[0], event.jets_p4[1])));
        VAR("dphi_b1b2", ROOT::Math::VectorUtil::DeltaPhi(event.jets_p4[0], event.jets_p4[1]));
        VAR("abs_dphi_l1met", std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.pfMET_p4)));
        VAR("dphi_l1met", ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.pfMET_p4));
        VAR("abs_dphi_l2met", std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_2, event.pfMET_p4)));
        VAR("dphi_l2met", ROOT::Math::VectorUtil::DeltaPhi(event.p4_2, event.pfMET_p4));
        VAR("abs_dphi_l1l2met", std::abs(ROOT::Math::VectorUtil::DeltaPhi(leptons, event.pfMET_p4)));
        VAR("dphi_l1l2met", ROOT::Math::VectorUtil::DeltaPhi(leptons, event.pfMET_p4));
        VAR("abs_dphi_htautaumet", std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4)));
        VAR("dphi_htautaumet", ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4));
        VAR("abs_dphi_hbbmet", std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.pfMET_p4)));
        VAR("dphi_hbbmet", ROOT::Math::VectorUtil::DeltaPhi(bb, event.pfMET_p4));
        VAR("abs_dphi_hbbhatutau", std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.SVfit_p4)));
        VAR("dphi_hbbhtautau", ROOT::Math::VectorUtil::DeltaPhi(bb, event.SVfit_p4));
        VAR("abs_deta_l1l2", std::abs(event.p4_1.eta() - event.p4_2.eta()));
        VAR("deta_l1l2", event.p4_1.eta() - event.p4_2.eta());
        VAR("abs_deta_b1b2", std::abs(event.jets_p4[0].eta() - event.jets_p4[1].eta()));
        VAR("deta_b1b2", event.jets_p4[0].eta() - event.jets_p4[1].eta());
        VAR("abs_deta_l1met", std::abs(event.p4_1.eta()-event.pfMET_p4.eta()));
        VAR("deta_l1met", event.p4_1.eta()-event.pfMET_p4.eta());
        VAR("abs_deta_l2met", std::abs(event.p4_2.eta()-event.pfMET_p4.eta()));
        VAR("deta_l2met", event.p4_2.eta()-event.pfMET_p4.eta());
        VAR("abs_deta_l1l2met", std::abs(leptons.eta()-event.pfMET_p4.eta()));
        VAR("deta_l1l2met", leptons.eta()-event.pfMET_p4.eta());
        VAR("abs_deta_htautaumet", std::abs(event.SVfit_p4.eta()-event.pfMET_p4.eta()));
        VAR("deta_hatutaumet", event.SVfit_p4.eta()-event.pfMET_p4.eta());
        VAR("abs_deta_hbbmet", std::abs(bb.eta()-event.pfMET_p4.eta()));
        VAR("deta_hbbmet", bb.eta()-event.pfMET_p4.eta());
        VAR("abs_deta_hbbhtautau", std::abs(bb.eta()-event.SVfit_p4.eta()));
        VAR("deta_hbbhtautau", bb.eta()-event.SVfit_p4.eta());
        VAR("dR_l1l2", ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2));
        VAR("dR_b1b2", ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]));
        VAR("dR_l1met", ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.pfMET_p4));
        VAR("dR_l2met", ROOT::Math::VectorUtil::DeltaR(event.p4_2, event.pfMET_p4));
        VAR("dR_l1l2met", ROOT::Math::VectorUtil::DeltaR(leptons, event.pfMET_p4));
        VAR("dR_htautaumet", ROOT::Math::VectorUtil::DeltaR(event.SVfit_p4, event.pfMET_p4));
        VAR("dR_hbbmet", ROOT::Math::VectorUtil::DeltaR(bb, event.pfMET_p4));
        VAR("dR_hbbhtautau", ROOT::Math::VectorUtil::DeltaR(bb, event.SVfit_p4));
        VAR("dR_b1b2Pt_hbb", (ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]))*bb.Pt());
        VAR("dR_l1l2Pt_htautau", ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2)*event.SVfit_p4.Pt());
        VAR("mass_l1l2met", ROOT::Math::VectorUtil::InvariantMass(leptons,event.pfMET_p4));
        VAR("mass_htautau", event.SVfit_p4.M());
        VAR("mass_l1l2", (event.p4_1+event.p4_2).M());
        VAR("mass_hbb", bb.M());
        VAR("MT_l1", Calculate_MT(event.p4_1,event.pfMET_p4));
        VAR("MT_l2", Calculate_MT(event.p4_2,event.pfMET_p4));
        VAR("MT_htautau", Calculate_MT(event.SVfit_p4, event.pfMET_p4));
        VAR("MT_l1l2", Calculate_MT(leptons, event.pfMET_p4));
        VAR("MT_tot", Calculate_TotalMT(event.p4_1,event.p4_2,event.pfMET_p4)); //Total transverse mass
        VAR("MT2", std::min(Calculate_MT2(event.p4_1,event.p4_2,event.jets_p4[0], event.jets_p4[1], event.pfMET_p4),Calculate_MT2(event.p4_1, event.jets_p4[1], event.p4_2,event.jets_p4[0], event.pfMET_p4))); //Stransverse mass
//        VAR("mass_H", ROOT::Math::VectorUtil::InvariantMass(bb,event.SVfit_p4));
        VAR("mass_top1", four_bodies::Calculate_topPairMasses(event.p4_1, event.p4_2, event.jets_p4[0], event.jets_p4[1], event.pfMET_p4).first);
        VAR("mass_top2", four_bodies::Calculate_topPairMasses(event.p4_1, event.p4_2, event.jets_p4[0], event.jets_p4[1], event.pfMET_p4).second);
//        VAR("MX", four_bodies::Calculate_MX(event.p4_1, event.p4_2, event.jets_p4[0], event.jets_p4[1], event.pfMET_p4));
        VAR("dR_l1l2_boosted", four_bodies::Calculate_dR_boosted(event.p4_1, event.p4_2, event.SVfit_p4));
        VAR("dR_b1b2_boosted", four_bodies::Calculate_dR_boosted(event.jets_p4[0], event.jets_p4[1], bb));
        VAR("phi", four_bodies::Calculate_phi(event.p4_1,event.p4_2,event.jets_p4[0], event.jets_p4[1], event.SVfit_p4, bb));
        VAR("costheta_star_leptons", four_bodies::Calculate_cosThetaStar(event.SVfit_p4, bb));
        VAR("costheta_star_bjets", four_bodies::Calculate_cosThetaStar(bb, event.SVfit_p4));
        VAR("phi_1", four_bodies::Calculate_phi1(event.p4_1, event.p4_2, event.SVfit_p4, bb));
        VAR("phi_2", four_bodies::Calculate_phi1(event.jets_p4[0], event.jets_p4[1], bb, event.SVfit_p4));
        VAR("costheta_l1htautau", four_bodies::Calculate_cosTheta_2bodies(event.p4_1, event.SVfit_p4));
        VAR("costheta_b1hbb", four_bodies::Calculate_cosTheta_2bodies(event.jets_p4[0], bb));
        VAR("costheta_hbbhh", four_bodies::Calculate_cosTheta_h_hh(bb, bb+event.SVfit_p4));
        VAR("costheta_hbbhhmet", four_bodies::Calculate_cosTheta_h_hh(bb, bb+leptonsMET));
        VAR("costheta_htautauhh", four_bodies::Calculate_cosTheta_h_hh(event.SVfit_p4, bb+event.SVfit_p4));
        VAR("costheta_l1l2methh", four_bodies::Calculate_cosTheta_h_hh(leptonsMET, bb+event.SVfit_p4));
        VAR("costheta_l1l2methhmet", four_bodies::Calculate_cosTheta_h_hh(leptonsMET, bb+leptonsMET));


        const bool is_training = split_training_testing ? test_vs_training(gen) == 1 : true;
        AddEventVariables(is_training, mass, tree, sample_weight); // event.weight * sample_weight
    }

private:
    bool split_training_testing;
    std::mt19937 gen;
    std::uniform_int_distribution<> test_vs_training;
    VarNameSet enabled_vars, disabled_vars;
};
#undef VAR

}
