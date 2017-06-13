/*! Definition of MvaVariables.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <iostream>
#include <random>
#include "h-tautau/Analysis/include/EventTuple.h"
#include "hh-bbtautau/Analysis/include/MT2.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/TextIO.h"

namespace analysis {
namespace mva_study{

//ENUM_ISTREAM_OPERATORS()
ENUM_OSTREAM_OPERATORS()

enum class SampleType { Sgn_Res = 1, Sgn_NonRes = 0, Bkg_TTbar = -1 };

struct SampleId {
    SampleType sampleType;
    int mass;

    SampleId() : sampleType(SampleType::Sgn_Res), mass(0) {}
    SampleId(SampleType _sampleType, int _mass = 0) : sampleType(_sampleType), mass(_mass) {}

    bool operator<(const SampleId& x) const
    {
        if (sampleType != x.sampleType) return static_cast<int>(sampleType) < static_cast<int>(x.sampleType);
        return mass < x.mass;
    }

    bool IsSignal() const { return sampleType == SampleType::Sgn_Res || sampleType == SampleType::Sgn_NonRes; }
    bool IsBackground() const { return !IsSignal(); }
    bool IsSM() const { return sampleType == SampleType::Sgn_NonRes; }

};

//static const SampleId Bkg{SampleType::Bkg_TTbar, -1};

ENUM_NAMES(SampleType) = {
    {SampleType::Sgn_Res, "Sgn_Res"},
    {SampleType::Sgn_NonRes, "SM"},
    {SampleType::Bkg_TTbar, "TT"}
};


inline std::ostream& operator<<(std::ostream& os, const SampleId& id)
{
    if(id.sampleType == SampleType::Sgn_NonRes || id.sampleType == SampleType::Bkg_TTbar)
        os << id.sampleType;
    else
        os << "M" << id.mass;
    return os;
}

inline std::istream& operator>>(std::istream& is, SampleId& id)
{
    std::string type;
    is >> type;
    id.mass = 0;
    if(!TryParse(type, id.sampleType)) {
        if(!type.size() || type.at(0) != 'M')
            throw exception("Bad sample id");
        id.sampleType = SampleType::Sgn_Res;
        id.mass = Parse<int>(type.substr(1));
    }/* else {
        if (type.at(0) == 'S'){
            id.sampleType = SampleType::Sgn_NonRes;
            id.mass = 0;
        }
        else{
            if (type.at(0) == 'T'){
                id.sampleType = SampleType::Bkg_TTbar;
                id.mass = -1;
            }
        }
    }*/
    return is;
}

#define VAR(name, formula) if(IsEnabled(name)) SetValue(name, formula)
class MvaVariables {
public:
    using VarNameSet = std::unordered_set<std::string>;

    MvaVariables(size_t _number_set = 1, uint_fast32_t seed = std::numeric_limits<uint_fast32_t>::max(),
                 const  VarNameSet& _enabled_vars = {}) :
        gen(seed), which_set(0, _number_set-1), enabled_vars(_enabled_vars)
    {
    }

    virtual ~MvaVariables() {}
    virtual void SetValue(const std::string& name, double value) = 0;
    virtual void AddEventVariables(size_t which_set, const SampleId& mass, double weight) = 0;

    bool IsEnabled(const std::string& name) const
    {
        return !enabled_vars.size() || enabled_vars.count(name);
    }

    void AddEvent(const ntuple::Event& event, const SampleId& mass , double sample_weight = 1.)
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
        VAR("pt_l1l2MET", leptonsMET.pt());
        VAR("pt_htautau", event.SVfit_p4.pt());
        VAR("pt_MET", event.pfMET_p4.pt());
        VAR("HT_otherjets", Calculate_HT(event.jets_p4.begin()+2, event.jets_p4.end()));
        VAR("p_zeta", Calculate_Pzeta(event.p4_1, event.p4_2, event.pfMET_p4));
        VAR("p_zetavisible", Calculate_visiblePzeta(event.p4_1,event.p4_2));
        VAR("abs_dphi_l1l2", std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.p4_2)));
        VAR("dphi_l1l2", ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.p4_2));
        VAR("abs_dphi_b1b2", std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.jets_p4[0], event.jets_p4[1])));
        VAR("dphi_b1b2", ROOT::Math::VectorUtil::DeltaPhi(event.jets_p4[0], event.jets_p4[1]));
        VAR("abs_dphi_l1MET", std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.pfMET_p4)));
        VAR("dphi_l1MET", ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.pfMET_p4));
        VAR("abs_dphi_l2MET", std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_2, event.pfMET_p4)));
        VAR("dphi_l2MET", ROOT::Math::VectorUtil::DeltaPhi(event.p4_2, event.pfMET_p4));
        VAR("abs_dphi_l1l2MET", std::abs(ROOT::Math::VectorUtil::DeltaPhi(leptons, event.pfMET_p4)));
        VAR("dphi_l1l2MET", ROOT::Math::VectorUtil::DeltaPhi(leptons, event.pfMET_p4));
        VAR("abs_dphi_htautauMET", std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4)));
        VAR("dphi_htautauMET", ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4));
        VAR("abs_dphi_hbbMET", std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.pfMET_p4)));
        VAR("dphi_hbbMET", ROOT::Math::VectorUtil::DeltaPhi(bb, event.pfMET_p4));
        VAR("abs_dphi_hbbhatutau", std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.SVfit_p4)));
        VAR("dphi_hbbhtautau", ROOT::Math::VectorUtil::DeltaPhi(bb, event.SVfit_p4));
        VAR("abs_deta_l1l2", std::abs(event.p4_1.eta() - event.p4_2.eta()));
        VAR("deta_l1l2", event.p4_1.eta() - event.p4_2.eta());
        VAR("abs_deta_b1b2", std::abs(event.jets_p4[0].eta() - event.jets_p4[1].eta()));
        VAR("deta_b1b2", event.jets_p4[0].eta() - event.jets_p4[1].eta());
        VAR("abs_deta_l1MET", std::abs(event.p4_1.eta()-event.pfMET_p4.eta()));
        VAR("deta_l1MET", event.p4_1.eta()-event.pfMET_p4.eta());
        VAR("abs_deta_l2MET", std::abs(event.p4_2.eta()-event.pfMET_p4.eta()));
        VAR("deta_l2MET", event.p4_2.eta()-event.pfMET_p4.eta());
        VAR("abs_deta_l1l2MET", std::abs(leptons.eta()-event.pfMET_p4.eta()));
        VAR("deta_l1l2MET", leptons.eta()-event.pfMET_p4.eta());
        VAR("abs_deta_htautauMET", std::abs(event.SVfit_p4.eta()-event.pfMET_p4.eta()));
        VAR("deta_htautauMET", event.SVfit_p4.eta()-event.pfMET_p4.eta());
        VAR("abs_deta_hbbMET", std::abs(bb.eta()-event.pfMET_p4.eta()));
        VAR("deta_hbbMET", bb.eta()-event.pfMET_p4.eta());
        VAR("abs_deta_hbbhtautau", std::abs(bb.eta()-event.SVfit_p4.eta()));
        VAR("deta_hbbhtautau", bb.eta()-event.SVfit_p4.eta());
        VAR("dR_l1l2", ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2));
        VAR("dR_b1b2", ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]));
        VAR("dR_l1MET", ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.pfMET_p4));
        VAR("dR_l2MET", ROOT::Math::VectorUtil::DeltaR(event.p4_2, event.pfMET_p4));
        VAR("dR_l1l2MET", ROOT::Math::VectorUtil::DeltaR(leptons, event.pfMET_p4));
        VAR("dR_htautauMET", ROOT::Math::VectorUtil::DeltaR(event.SVfit_p4, event.pfMET_p4));
        VAR("dR_hbbMET", ROOT::Math::VectorUtil::DeltaR(bb, event.pfMET_p4));
        VAR("dR_hbbhtautau", ROOT::Math::VectorUtil::DeltaR(bb, event.SVfit_p4));
        VAR("dR_b1b2Pt_hbb", (ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]))*bb.Pt());
        VAR("dR_l1l2Pt_htautau", ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2)*event.SVfit_p4.Pt());
        VAR("mass_l1l2MET", ROOT::Math::VectorUtil::InvariantMass(leptons,event.pfMET_p4));
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
        VAR("phi_1", four_bodies::Calculate_phi1(event.p4_1, event.p4_2, event.SVfit_p4, bb));
        VAR("phi_2", four_bodies::Calculate_phi1(event.jets_p4[0], event.jets_p4[1], bb, event.SVfit_p4));
        VAR("costheta_l1htautau", four_bodies::Calculate_cosTheta_2bodies(event.p4_1, event.SVfit_p4));
        VAR("costheta_l2htautau", four_bodies::Calculate_cosTheta_2bodies(event.p4_2, event.SVfit_p4));
        VAR("costheta_METhtautau", four_bodies::Calculate_cosTheta_2bodies(event.pfMET_p4, event.SVfit_p4));
        VAR("costheta_METhbb", four_bodies::Calculate_cosTheta_2bodies(event.pfMET_p4, bb));
        VAR("costheta_b1hbb", four_bodies::Calculate_cosTheta_2bodies(event.jets_p4[0], bb));
        VAR("costheta_hbbhh", four_bodies::Calculate_cosTheta_2bodies(bb, bb+event.SVfit_p4));
        VAR("costheta_hbbhhMET", four_bodies::Calculate_cosTheta_2bodies(bb, bb+leptonsMET));
        VAR("costheta_htautauhh", four_bodies::Calculate_cosTheta_2bodies(event.SVfit_p4, bb+event.SVfit_p4));
        VAR("costheta_l1l2METhh", four_bodies::Calculate_cosTheta_2bodies(leptonsMET, bb+event.SVfit_p4));
        VAR("costheta_l1l2METhhMET", four_bodies::Calculate_cosTheta_2bodies(leptonsMET, bb+leptonsMET));

        AddEventVariables(which_set(gen), mass, sample_weight); // event.weight * sample_weight
    }

private:
    std::mt19937 gen;
    std::uniform_int_distribution<size_t> which_set;
    VarNameSet enabled_vars;
};
#undef VAR
}
}
