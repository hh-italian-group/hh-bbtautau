/*! Definition of MvaVariables.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <iostream>
#include <random>
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"
#include "hh-bbtautau/Analysis/include/MT2.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "TMVA/Reader.h"

namespace analysis {
namespace mva_study{

using ::analysis::operator<<;
using ::analysis::operator>>;

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

    bool operator ==(const SampleId& x) const
    {
        return x.mass==mass && x.sampleType == sampleType;
    }

    bool operator !=(const SampleId& x) const {return !(x==*this);}

    bool IsSignal() const { return sampleType == SampleType::Sgn_Res || sampleType == SampleType::Sgn_NonRes; }
    bool IsBackground() const { return !IsSignal(); }
    bool IsSM() const { return sampleType == SampleType::Sgn_NonRes; }

    static const SampleId& MassTot()
    {
        static const SampleId mass_tot(SampleType::Sgn_Res, std::numeric_limits<int>::max());
        return mass_tot;
    }
    static const SampleId& Bkg()
    {
        static const SampleId bkg(SampleType::Bkg_TTbar, 0);
        return bkg;
    }
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
    else {
        if(id == SampleId::MassTot())
            os << "Mtot";
        else
            os << "M" << id.mass;
    }
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const std::pair<SampleId, SampleId>& id_pair)
{
    os << id_pair.first << ":" << id_pair.second;
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
    }
    return is;
}


class MvaVariablesBase {
public:
    virtual ~MvaVariablesBase() {}
    virtual void AddEvent(analysis::EventInfoBase& eventbase, const SampleId& mass , int spin, double sample_weight = 1., int which_test = -1) = 0;
    virtual double Evaluate() { throw exception("Not supported."); }
    virtual std::shared_ptr<TMVA::Reader> GetReader() = 0;
};

#define VAR(name, formula) if(IsEnabled(name)) SetValue(name, formula)
#define VAR_INT(name, formula) if(IsEnabled(name)) SetValue(name, formula, 'I')
class MvaVariables : public MvaVariablesBase {
public:
    using VarNameSet = std::unordered_set<std::string>;
    MvaVariables(size_t _number_set = 1, uint_fast32_t seed = std::numeric_limits<uint_fast32_t>::max(),
                 const VarNameSet& _enabled_vars = {}, const VarNameSet& _disabled_vars = {}) :
        gen(seed), which_set(0, _number_set-1), enabled_vars(_enabled_vars), disabled_vars(_disabled_vars)
    {
    }

    virtual ~MvaVariables() {}
    virtual void SetValue(const std::string& name, double value, char type = 'F') = 0;
    virtual void AddEventVariables(size_t which_set, const SampleId& mass, double weight, double sampleweight, int spin, std::string channel) = 0;
    bool IsEnabled(const std::string& name) const
    {
        return (!enabled_vars.size() && !disabled_vars.count(name)) || enabled_vars.count(name);
    }

    virtual void AddEvent(analysis::EventInfoBase& eventbase, const SampleId& mass , int spin, double sample_weight = 1., int which_test = -1) override
    {
        ntuple::Event event = *eventbase;

        const auto& Htt = eventbase.GetHiggsTTMomentum(false);
        const auto& Htt_sv = eventbase.GetHiggsTTMomentum(true);
        const auto& t1 = eventbase.GetLeg(1);
        const auto& t2 = eventbase.GetLeg(2);

        const auto& Hbb = eventbase.GetHiggsBB();
        const auto& b1 = Hbb.GetFirstDaughter();
        const auto& b2 = Hbb.GetSecondDaughter();

        const auto& met = eventbase.GetMET();

        VAR("pt_b1", b1.GetMomentum().Pt());
        VAR("pt_b2", b2.GetMomentum().Pt());
        VAR("pt_hbb", Hbb.GetMomentum().Pt());
        VAR("pt_l1", t1.GetMomentum().Pt());
        VAR("pt_l2", t2.GetMomentum().Pt());
        VAR("pt_l1l2", (t1.GetMomentum()+t2.GetMomentum()).Pt());
        VAR("pt_htautau", Htt.Pt());
        VAR("pt_htautau_sv", Htt_sv.Pt());
        VAR("pt_MET", met.GetMomentum().Pt());
        VAR("HT_otherjets", event.ht_other_jets);
        VAR("p_zeta", Calculate_Pzeta(t1.GetMomentum(), t2.GetMomentum(),  met.GetMomentum()));
        VAR("p_zetavisible", Calculate_visiblePzeta(t1.GetMomentum(), t2.GetMomentum()));

        VAR("abs_dphi_l1l2", std::abs(ROOT::Math::VectorUtil::DeltaPhi(t1.GetMomentum(), t2.GetMomentum())));
        VAR("dphi_l1l2", ROOT::Math::VectorUtil::DeltaPhi(t1.GetMomentum(), t2.GetMomentum()));
        VAR("abs_dphi_b1b2", std::abs(ROOT::Math::VectorUtil::DeltaPhi(b1.GetMomentum(), b2.GetMomentum())));
        VAR("dphi_b1b2", ROOT::Math::VectorUtil::DeltaPhi(b1.GetMomentum(), b2.GetMomentum()));
        VAR("abs_dphi_l1MET", std::abs(ROOT::Math::VectorUtil::DeltaPhi(t1.GetMomentum(), met.GetMomentum())));
        VAR("dphi_l1MET", ROOT::Math::VectorUtil::DeltaPhi(t1.GetMomentum(), met.GetMomentum()));
        VAR("abs_dphi_l2MET", std::abs(ROOT::Math::VectorUtil::DeltaPhi(t2.GetMomentum(), met.GetMomentum())));
        VAR("dphi_l2MET", ROOT::Math::VectorUtil::DeltaPhi(t2.GetMomentum(), met.GetMomentum()));
        VAR("abs_dphi_l1l2MET", std::abs(ROOT::Math::VectorUtil::DeltaPhi(t1.GetMomentum()+t2.GetMomentum(), met.GetMomentum())));
        VAR("dphi_l1l2MET", ROOT::Math::VectorUtil::DeltaPhi(t1.GetMomentum()+t2.GetMomentum(), met.GetMomentum()));
        VAR("abs_dphi_METhtautau", std::abs(ROOT::Math::VectorUtil::DeltaPhi(Htt, met.GetMomentum())));
        VAR("dphi_METhtautau", ROOT::Math::VectorUtil::DeltaPhi(Htt, met.GetMomentum()));
        VAR("abs_dphi_METhtautau_sv", std::abs(ROOT::Math::VectorUtil::DeltaPhi(Htt_sv, met.GetMomentum())));
        VAR("dphi_METhtautau_sv", ROOT::Math::VectorUtil::DeltaPhi(Htt_sv, met.GetMomentum()));
        VAR("abs_dphi_hbbMET", std::abs(ROOT::Math::VectorUtil::DeltaPhi(Hbb.GetMomentum(), met.GetMomentum())));
        VAR("dphi_hbbMET", ROOT::Math::VectorUtil::DeltaPhi(Hbb.GetMomentum(), met.GetMomentum()));
        VAR("abs_dphi_hbbhatutau", std::abs(ROOT::Math::VectorUtil::DeltaPhi(Hbb.GetMomentum(), Htt)));
        VAR("dphi_hbbhtautau", ROOT::Math::VectorUtil::DeltaPhi(Hbb.GetMomentum(), Htt));
        VAR("abs_dphi_hbbhatutau_sv", std::abs(ROOT::Math::VectorUtil::DeltaPhi(Hbb.GetMomentum(), Htt_sv)));
        VAR("dphi_hbbhtautau_sv", ROOT::Math::VectorUtil::DeltaPhi(Hbb.GetMomentum(), Htt_sv));

        VAR("abs_deta_l1l2", std::abs(t1.GetMomentum().eta() - t2.GetMomentum().eta()));
        VAR("deta_l1l2", t1.GetMomentum().eta() - t2.GetMomentum().eta());
        VAR("abs_deta_b1b2", std::abs(b1.GetMomentum().eta() - b2.GetMomentum().eta()));
        VAR("deta_b1b2", b1.GetMomentum().eta() - b2.GetMomentum().eta());
        VAR("abs_deta_l1MET", std::abs(t1.GetMomentum().eta()-met.GetMomentum().eta()));
        VAR("deta_l1MET", t1.GetMomentum().eta()-met.GetMomentum().eta());
        VAR("abs_deta_l2MET", std::abs(t2.GetMomentum().eta()-met.GetMomentum().eta()));
        VAR("deta_l2MET", t2.GetMomentum().eta()-met.GetMomentum().eta());
        VAR("abs_deta_l1l2MET", std::abs((t1.GetMomentum()+t2.GetMomentum()).eta()-met.GetMomentum().eta()));
        VAR("deta_l1l2MET", (t1.GetMomentum()+t2.GetMomentum()).eta()-met.GetMomentum().eta());
        VAR("abs_deta_METhtautau", std::abs(Htt.eta()-met.GetMomentum().eta()));
        VAR("deta_METhtautau", Htt.eta()-met.GetMomentum().eta());
        VAR("abs_deta_METhtautau_sv", std::abs(Htt_sv.eta()-met.GetMomentum().eta()));
        VAR("deta_METhtautau_sv", Htt_sv.eta()-met.GetMomentum().eta());
        VAR("abs_deta_hbbMET", std::abs(Hbb.GetMomentum().eta()-met.GetMomentum().eta()));
        VAR("deta_hbbMET", Hbb.GetMomentum().eta()-met.GetMomentum().eta());
        VAR("abs_deta_hbbhtautau", std::abs(Hbb.GetMomentum().eta()-Htt.eta()));
        VAR("deta_hbbhtautau", Hbb.GetMomentum().eta()-Htt.eta());
        VAR("abs_deta_hbbhtautau_sv", std::abs(Hbb.GetMomentum().eta()-Htt_sv.eta()));
        VAR("deta_hbbhtautau_sv", Hbb.GetMomentum().eta()-Htt_sv.eta());

        VAR("dR_l1l2", ROOT::Math::VectorUtil::DeltaR(t1.GetMomentum(), t2.GetMomentum()));
        VAR("dR_b1b2", ROOT::Math::VectorUtil::DeltaR(b1.GetMomentum(), b2.GetMomentum()));
        VAR("dR_l1MET", ROOT::Math::VectorUtil::DeltaR(t1.GetMomentum(), met.GetMomentum()));
        VAR("dR_l2MET", ROOT::Math::VectorUtil::DeltaR(t2.GetMomentum(), met.GetMomentum()));
        VAR("dR_l1l2MET", ROOT::Math::VectorUtil::DeltaR(t1.GetMomentum()+t2.GetMomentum(), met.GetMomentum()));
        VAR("dR_METhtautau", ROOT::Math::VectorUtil::DeltaR(Htt, met.GetMomentum()));
        VAR("dR_METhtautau_sv", ROOT::Math::VectorUtil::DeltaR(Htt_sv, met.GetMomentum()));
        VAR("dR_hbbMET", ROOT::Math::VectorUtil::DeltaR(Hbb.GetMomentum(), met.GetMomentum()));
        VAR("dR_hbbhtautau", ROOT::Math::VectorUtil::DeltaR(Hbb.GetMomentum(), Htt));
        VAR("dR_hbbhtautau_sv", ROOT::Math::VectorUtil::DeltaR(Hbb.GetMomentum(), Htt_sv));
        VAR("dR_b1b2Pt_hbb", (ROOT::Math::VectorUtil::DeltaR(b1.GetMomentum(), b2.GetMomentum()))*Hbb.GetMomentum().Pt());
        VAR("dR_l1l2Pt_htautau", ROOT::Math::VectorUtil::DeltaR(t1.GetMomentum(), t2.GetMomentum())*Htt.Pt());
        VAR("dR_l1l2Pt_htautau_sv", ROOT::Math::VectorUtil::DeltaR(t1.GetMomentum(), t2.GetMomentum())*Htt_sv.Pt());
        VAR("dR_b1b2_boosted", four_bodies::Calculate_dR_boosted(b1.GetMomentum(), b2.GetMomentum(), Hbb.GetMomentum()));
        VAR("dR_l1l2_boosted", four_bodies::Calculate_dR_boosted(t1.GetMomentum(), t2.GetMomentum(), Htt));
        VAR("dR_l1l2_boosted_sv", four_bodies::Calculate_dR_boosted(t1.GetMomentum(), t2.GetMomentum(), Htt_sv));

        VAR("mass_l1l2MET", ROOT::Math::VectorUtil::InvariantMass(t1.GetMomentum()+t2.GetMomentum(), met.GetMomentum()));
        VAR("mass_htautau", Htt.M());
        VAR("mass_htautau_sv", Htt_sv.M());
        VAR("mass_l1l2", (t1.GetMomentum()+t2.GetMomentum()).M());
        VAR("mass_hbb", Hbb.GetMomentum().M());
        VAR("MT_l1", Calculate_MT(t1.GetMomentum(),met.GetMomentum()));
        VAR("MT_l2", Calculate_MT(t2.GetMomentum(),met.GetMomentum()));
        VAR("MT_htautau", Calculate_MT(Htt, met.GetMomentum()));
        VAR("MT_htautau_sv", Calculate_MT(Htt_sv, met.GetMomentum()));
        VAR("MT_l1l2", Calculate_MT(t1.GetMomentum()+t2.GetMomentum(), met.GetMomentum()));
        VAR("MT_tot", Calculate_TotalMT(t1.GetMomentum(), t2.GetMomentum(),met.GetMomentum())); //Total transverse mass
        VAR("MT2", eventbase.GetMT2()); //Stransverse mass
//        VAR("MT2", std::min(Calculate_MT2_old(event.p4_1, event.p4_2, event.jets_p4[0], event.jets_p4[1], event.pfMET_p4), Calculate_MT2_old(event.p4_1, event.p4_2, event.jets_p4[1], event.jets_p4[0], event.pfMET_p4))); //Stransverse mass
        VAR("mass_top1", four_bodies::Calculate_topPairMasses(t1.GetMomentum(), t2.GetMomentum(), b1.GetMomentum(), b2.GetMomentum(), met.GetMomentum()).first);
        VAR("mass_top2", four_bodies::Calculate_topPairMasses(t1.GetMomentum(), t2.GetMomentum(), b1.GetMomentum(), b2.GetMomentum(), met.GetMomentum()).second);
        VAR("mass_X", four_bodies::Calculate_MX(t1.GetMomentum(), t2.GetMomentum(), b1.GetMomentum(), b2.GetMomentum(), met.GetMomentum()));
        VAR("mass_H", ROOT::Math::VectorUtil::InvariantMass(Hbb.GetMomentum(), Htt));
        VAR("mass_H_sv", ROOT::Math::VectorUtil::InvariantMass(Hbb.GetMomentum(), Htt_sv));
        VAR("mass_H_vis", ROOT::Math::VectorUtil::InvariantMass(Hbb.GetMomentum(), t1.GetMomentum()+t2.GetMomentum()));
        VAR("mass_H_kinfit", eventbase.GetKinFitResults().mass);
        VAR("mass_H_kinfit_chi2", eventbase.GetKinFitResults().chi2);

        VAR("phi", four_bodies::Calculate_phi(t1.GetMomentum(), t2.GetMomentum(), b1.GetMomentum(), b2.GetMomentum(), Htt, Hbb.GetMomentum()));
        VAR("phi_sv", four_bodies::Calculate_phi(t1.GetMomentum(), t2.GetMomentum(), b1.GetMomentum(), b2.GetMomentum(), Htt_sv, Hbb.GetMomentum()));
        VAR("phi_1", four_bodies::Calculate_phi1(t1.GetMomentum(), t2.GetMomentum(), Htt, Hbb.GetMomentum()));
        VAR("phi_1_sv", four_bodies::Calculate_phi1(t1.GetMomentum(), t2.GetMomentum(), Htt_sv, Hbb.GetMomentum()));
        VAR("phi_2", four_bodies::Calculate_phi1(b1.GetMomentum(), b2.GetMomentum(), Htt, Hbb.GetMomentum()));
        VAR("phi_2_sv", four_bodies::Calculate_phi1(b1.GetMomentum(), b2.GetMomentum(), Htt_sv, Hbb.GetMomentum()));
        VAR("costheta_star_leptons", four_bodies::Calculate_cosThetaStar(Htt, Hbb.GetMomentum()));
        VAR("costheta_star_leptons_sv", four_bodies::Calculate_cosThetaStar(Htt_sv, Hbb.GetMomentum()));
        VAR("costheta_l1htautau", four_bodies::Calculate_cosTheta_2bodies(t1.GetMomentum(), Htt));
        VAR("costheta_l1htautau_sv", four_bodies::Calculate_cosTheta_2bodies(t1.GetMomentum(), Htt_sv));
        VAR("costheta_l2htautau", four_bodies::Calculate_cosTheta_2bodies(t2.GetMomentum(), Htt));
        VAR("costheta_l2htautau_sv", four_bodies::Calculate_cosTheta_2bodies(t2.GetMomentum(), Htt_sv));
        VAR("costheta_METhtautau", four_bodies::Calculate_cosTheta_2bodies(met.GetMomentum(), Htt));
        VAR("costheta_METhtautau_sv", four_bodies::Calculate_cosTheta_2bodies(met.GetMomentum(), Htt_sv));
        VAR("costheta_METhbb", four_bodies::Calculate_cosTheta_2bodies(met.GetMomentum(), Hbb.GetMomentum()));
        VAR("costheta_b1hbb", four_bodies::Calculate_cosTheta_2bodies(b1.GetMomentum(), Hbb.GetMomentum()));

        VAR("costheta_hbbhh", four_bodies::Calculate_cosTheta_2bodies(Hbb.GetMomentum(), eventbase.GetResonanceMomentum(false,false)));
        VAR("costheta_hbbhh_sv", four_bodies::Calculate_cosTheta_2bodies(Hbb.GetMomentum(), eventbase.GetResonanceMomentum(true,false)));
        VAR("costheta_hbbhhMET", four_bodies::Calculate_cosTheta_2bodies(Hbb.GetMomentum(), eventbase.GetResonanceMomentum(false,true)));

        VAR("costheta_htautauhh", four_bodies::Calculate_cosTheta_2bodies(Htt, eventbase.GetResonanceMomentum(false,false)));
        VAR("costheta_htautauhh_sv", four_bodies::Calculate_cosTheta_2bodies(Htt, eventbase.GetResonanceMomentum(true,false)));
        VAR("costheta_htautauhhMET", four_bodies::Calculate_cosTheta_2bodies(Htt, eventbase.GetResonanceMomentum(false,true)));
        VAR("costheta_htautau_svhh", four_bodies::Calculate_cosTheta_2bodies(Htt_sv, eventbase.GetResonanceMomentum(false,false)));
        VAR("costheta_htautau_svhh_sv", four_bodies::Calculate_cosTheta_2bodies(Htt_sv, eventbase.GetResonanceMomentum(true,false)));
        VAR("costheta_htautau_svhhMET", four_bodies::Calculate_cosTheta_2bodies(Htt_sv, eventbase.GetResonanceMomentum(false,true)));

        VAR("costheta_l1l2METhh", four_bodies::Calculate_cosTheta_2bodies(t1.GetMomentum()+t2.GetMomentum()+met.GetMomentum(), eventbase.GetResonanceMomentum(false,false)));
        VAR("costheta_l1l2METhh_sv", four_bodies::Calculate_cosTheta_2bodies(t1.GetMomentum()+t2.GetMomentum()+met.GetMomentum(), eventbase.GetResonanceMomentum(true,false)));
        VAR("costheta_l1l2METhhMET", four_bodies::Calculate_cosTheta_2bodies(t1.GetMomentum()+t2.GetMomentum()+met.GetMomentum(), eventbase.GetResonanceMomentum(false,true)));

        VAR("mass", mass.mass);
        VAR_INT("channel", event.channelId);
        VAR_INT("spin", spin);

        size_t test = which_test ==-1 ? which_set(gen) : static_cast<size_t>(which_test);
        AddEventVariables(test, mass, event.weight_total, sample_weight, spin, ToString(event.channelId));
    }

private:
    std::mt19937_64 gen;
    std::uniform_int_distribution<size_t> which_set;
    VarNameSet enabled_vars, disabled_vars;
};
#undef VAR
}
}
