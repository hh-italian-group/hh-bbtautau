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
    virtual void AddEventVariables(bool issignal, bool istraining, int mass, double weight) = 0;

    bool IsEnabled(const std::string& name) const
    {
        return (enabled_vars.size() && enabled_vars.count(name)) || !disabled_vars.count(name);
    }

    void AddEvent(const ntuple::Event& event, bool is_signal, int mass = 1, double sample_weight = 1.)
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
        VAR("pt_otherjets", Calculate_sumPt_otherJet(event.jets_p4));
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
        VAR("deta_hbbhatutau", bb.eta()-event.SVfit_p4.eta());
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
        VAR("MT_hatautau", Calculate_MT(event.SVfit_p4, event.pfMET_p4));
        VAR("MT_l1l2", Calculate_MT(leptons, event.pfMET_p4));
        VAR("MT_tot", Calculate_TotalMT(event.p4_1,event.p4_2,event.pfMET_p4)); //Total transverse mass
        VAR("MT2", std::min(Calculate_MT2(event.p4_1,event.p4_2,event.jets_p4[0], event.jets_p4[1], event.pfMET_p4),Calculate_MT2(event.p4_1, event.jets_p4[1], event.p4_2,event.jets_p4[0], event.pfMET_p4))); //Stransverse mass
        VAR("mass_H", ROOT::Math::VectorUtil::InvariantMass(bb,event.SVfit_p4));
        VAR("mass_top1", Calculate_topPairMasses(event.p4_1, event.p4_2, event.jets_p4[0], event.jets_p4[1], event.pfMET_p4).first);
        VAR("mass_top2", Calculate_topPairMasses(event.p4_1, event.p4_2, event.jets_p4[0], event.jets_p4[1], event.pfMET_p4).second);
        VAR("MX", Calculate_MX(event.p4_1, event.p4_2, event.jets_p4[0], event.jets_p4[1], event.pfMET_p4));
        VAR("dR_l1l2_boosted", Calculate_dR_boosted(event.p4_1, event.p4_2, event.SVfit_p4));
        VAR("dR_b1b2_boosted", Calculate_dR_boosted(event.jets_p4[0], event.jets_p4[1], bb));
        VAR("phi", Calculate_phi_4bodies(event.p4_1,event.p4_2,event.jets_p4[0], event.jets_p4[1], event.SVfit_p4, bb));
        VAR("theta_star_leptons", Calculate_cosThetaStar(event.SVfit_p4, bb));
        VAR("theta_star_bjets", Calculate_cosThetaStar(bb, event.SVfit_p4));
        VAR("phi_1", Calculate_phi1(event.p4_1, event.p4_2, event.SVfit_p4, bb));
        VAR("phi_2", Calculate_phi1(event.jets_p4[0], event.jets_p4[1], bb, event.SVfit_p4));
        VAR("theta_1", Calculate_theta_2bodies(event.p4_1, event.SVfit_p4));
        VAR("theta_2", Calculate_theta_2bodies(event.jets_p4[0], bb));

        const bool is_training = split_training_testing ? test_vs_training(gen) == 1 : true;
        AddEventVariables(is_signal, is_training, mass, sample_weight); // event.weight * sample_weight
    }

private:
    bool split_training_testing;
    std::mt19937 gen;
    std::uniform_int_distribution<> test_vs_training;
    VarNameSet enabled_vars, disabled_vars;
};
#undef VAR

}
