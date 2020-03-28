/*! Definition of MvaVariables.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/MvaVariables.h"

namespace analysis {
namespace mva_study{



SampleId::SampleId() : sampleType(SampleType::Sgn_Res), mass(0) {}
SampleId::SampleId(SampleType _sampleType, int _mass) : sampleType(_sampleType), mass(_mass) {}

bool SampleId::operator<(const SampleId& x) const
{
    if (sampleType != x.sampleType) return static_cast<int>(sampleType) < static_cast<int>(x.sampleType);
    return mass < x.mass;
}

bool SampleId::operator ==(const SampleId& x) const
{
    return x.mass==mass && x.sampleType == sampleType;
}

bool SampleId::operator !=(const SampleId& x) const {return !(x==*this);}

bool SampleId::IsSignal() const { return sampleType == SampleType::Sgn_Res || sampleType == SampleType::Sgn_NonRes; }
bool SampleId::IsBackground() const { return !IsSignal(); }
bool SampleId::IsSM() const { return sampleType == SampleType::Sgn_NonRes; }

const SampleId& SampleId::MassTot()
{
    static const SampleId mass_tot(SampleType::Sgn_Res, std::numeric_limits<int>::max());
    return mass_tot;
}

const SampleId& SampleId::Bkg()
{
    static const SampleId bkg(SampleType::Bkg_TTbar, 0);
    return bkg;
}

const SampleId& SampleId::SM()
{
    static const SampleId sm(SampleType::Sgn_NonRes, 0);
    return sm;
}

std::ostream& operator<<(std::ostream& os, const SampleId& id)
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

std::ostream& operator<<(std::ostream& os, const std::pair<SampleId, SampleId>& id_pair)
{
    os << id_pair.first << ":" << id_pair.second;
    return os;
}

std::istream& operator>>(std::istream& is, SampleId& id)
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

double MvaVariablesBase::Evaluate() { throw exception("Not supported."); }
double MvaVariablesBase::AddAndEvaluate(EventInfo& eventbase, const SampleId& mass, int spin,
                                        double sample_weight, int which_test)
{
    Lock lock(mutex);
    AddEvent(eventbase, mass, spin, sample_weight, which_test);
    return Evaluate();
}


#define VAR(name, formula) if(IsEnabled(name)) SetValue(name, formula)
#define VAR_INT(name, formula) if(IsEnabled(name)) SetValue(name, formula, 'I')

MvaVariables::MvaVariables(size_t _number_set, uint_fast32_t seed, const VarNameSet& _enabled_vars,
                           const VarNameSet& _disabled_vars) :
    gen(seed), which_set(0, _number_set-1), enabled_vars(_enabled_vars), disabled_vars(_disabled_vars)
{
}

const std::unordered_set<std::string>& MvaVariables::GetDisabledVars() const
{
    return disabled_vars;
}

bool MvaVariables::IsEnabled(const std::string& name) const
{
    return (!enabled_vars.size() && !disabled_vars.count(name)) || enabled_vars.count(name);
}

void MvaVariables::AddEvent(analysis::EventInfo& eventbase, const SampleId& mass, int spin, double sample_weight,
                            int which_test)
{
    using namespace ROOT::Math::VectorUtil;
    static constexpr double default_value = -999.;

    const auto& Htt = eventbase.GetHiggsTTMomentum(false);
    const auto& Htt_sv = eventbase.GetHiggsTTMomentum(true);
    const auto& t1 = eventbase.GetLeg(1).GetMomentum();
    const auto& t2 = eventbase.GetLeg(2).GetMomentum();

    const auto& Hbb = eventbase.GetHiggsBB().GetMomentum();
    const auto& b1 = eventbase.GetHiggsBB().GetFirstDaughter().GetMomentum();
    const auto& b2 = eventbase.GetHiggsBB().GetSecondDaughter().GetMomentum();

    const auto& met = eventbase.GetMET().GetMomentum();
    const bool SVfit_is_valid = eventbase.GetSVFitResults().has_valid_momentum;
    const auto& kinfit_results = eventbase.GetKinFitResults();
    const bool kinfit_is_valid = kinfit_results.HasValidMass();

    VAR_INT("decayMode_1", eventbase.GetLeg(1)->decayMode());
    VAR_INT("decayMode_2", eventbase.GetLeg(2)->decayMode());
    VAR("iso_1", eventbase.GetLeg(1)->iso());
    VAR("iso_2", eventbase.GetLeg(2)->iso());
    VAR("csv_1", eventbase.GetHiggsBB().GetFirstDaughter()->csv());
    VAR("csv_2", eventbase.GetHiggsBB().GetSecondDaughter()->csv());

    VAR("pt_b1", b1.Pt());
    VAR("pt_b2", b2.Pt());
    VAR("pt_hbb", Hbb.Pt());
    VAR("pt_l1", t1.Pt());
    VAR("pt_l2", t2.Pt());
    VAR("pt_l1l2", (t1+t2).Pt());
    VAR("pt_htautau_vis", Htt.Pt());
    VAR("pt_htautau", (Htt+met).Pt());
    VAR("pt_htautau_sv", SVfit_is_valid ? Htt_sv.Pt() : default_value);
    VAR("pt_MET", met.Pt());
    VAR("HT_otherjets", eventbase.GetHT(false,true));
    VAR("p_zeta", Calculate_Pzeta(t1, t2,  met));
    VAR("p_zetavisible", Calculate_visiblePzeta(t1, t2));

    VAR("abs_dphi_l1l2", std::abs(ROOT::Math::VectorUtil::DeltaPhi(t1, t2)));
    VAR("dphi_l1l2", ROOT::Math::VectorUtil::DeltaPhi(t1, t2));
    VAR("abs_dphi_b1b2", std::abs(ROOT::Math::VectorUtil::DeltaPhi(b1, b2)));
    VAR("dphi_b1b2", ROOT::Math::VectorUtil::DeltaPhi(b1, b2));
    VAR("abs_dphi_l1MET", std::abs(ROOT::Math::VectorUtil::DeltaPhi(t1, met)));
    VAR("dphi_l1MET", ROOT::Math::VectorUtil::DeltaPhi(t1, met));
    VAR("abs_dphi_l2MET", std::abs(ROOT::Math::VectorUtil::DeltaPhi(t2, met)));
    VAR("dphi_l2MET", ROOT::Math::VectorUtil::DeltaPhi(t2, met));
    VAR("abs_dphi_l1l2MET", std::abs(ROOT::Math::VectorUtil::DeltaPhi(t1+t2, met)));
    VAR("dphi_l1l2MET", ROOT::Math::VectorUtil::DeltaPhi(t1+t2, met));
    VAR("abs_dphi_METhtautau", std::abs(ROOT::Math::VectorUtil::DeltaPhi(Htt+met, met)));
    VAR("dphi_METhtautau", ROOT::Math::VectorUtil::DeltaPhi(Htt+met, met));
    VAR("abs_dphi_METhtautau_sv", SVfit_is_valid ? std::abs(ROOT::Math::VectorUtil::DeltaPhi(Htt_sv, met)) : default_value);
    VAR("dphi_METhtautau_sv", SVfit_is_valid ? ROOT::Math::VectorUtil::DeltaPhi(Htt_sv, met) : default_value);
    VAR("abs_dphi_hbbMET", std::abs(ROOT::Math::VectorUtil::DeltaPhi(Hbb, met)));
    VAR("dphi_hbbMET", ROOT::Math::VectorUtil::DeltaPhi(Hbb, met));
    VAR("abs_dphi_hbbhatutau_MET", std::abs(ROOT::Math::VectorUtil::DeltaPhi(Hbb, Htt+met)));
    VAR("dphi_hbbhtautau", ROOT::Math::VectorUtil::DeltaPhi(Hbb, Htt+met));
    VAR("abs_dphi_hbbhatutau_sv", SVfit_is_valid ? std::abs(ROOT::Math::VectorUtil::DeltaPhi(Hbb, Htt_sv)) : default_value);
    VAR("dphi_hbbhtautau_sv", SVfit_is_valid ? ROOT::Math::VectorUtil::DeltaPhi(Hbb, Htt_sv) : default_value);

    VAR("abs_deta_l1l2", std::abs(t1.eta() - t2.eta()));
    VAR("deta_l1l2", t1.eta() - t2.eta());
    VAR("abs_deta_b1b2", std::abs(b1.eta() - b2.eta()));
    VAR("deta_b1b2", b1.eta() - b2.eta());
    VAR("abs_deta_l1MET", std::abs(t1.eta()-met.eta()));
    VAR("deta_l1MET", t1.eta()-met.eta());
    VAR("abs_deta_l2MET", std::abs(t2.eta()-met.eta()));
    VAR("deta_l2MET", t2.eta()-met.eta());
    VAR("abs_deta_l1l2MET", std::abs((t1+t2).eta()-met.eta()));
    VAR("deta_l1l2MET", (t1+t2).eta()-met.eta());
    VAR("abs_deta_METhtautau", std::abs((Htt+met).eta()-met.eta()));
    VAR("deta_METhtautau", (Htt+met).eta()-met.eta());
    VAR("abs_deta_METhtautau_sv", SVfit_is_valid ? std::abs(Htt_sv.eta()-met.eta()) : default_value);
    VAR("deta_METhtautau_sv", SVfit_is_valid ? Htt_sv.eta()-met.eta() : default_value);
    VAR("abs_deta_hbbMET", std::abs(Hbb.eta()-met.eta()));
    VAR("deta_hbbMET", Hbb.eta()-met.eta());
    VAR("abs_deta_hbbhtautau", std::abs(Hbb.eta()-(Htt+met).eta()));
    VAR("deta_hbbhtautau", Hbb.eta()-(Htt+met).eta());
    VAR("abs_deta_hbbhtautau_sv", SVfit_is_valid ? std::abs(Hbb.eta()-Htt_sv.eta()) : default_value);
    VAR("deta_hbbhtautau_sv", SVfit_is_valid ? Hbb.eta()-Htt_sv.eta() : default_value);

    VAR("dR_l1l2", ROOT::Math::VectorUtil::DeltaR(t1, t2));
    VAR("dR_b1b2", ROOT::Math::VectorUtil::DeltaR(b1, b2));
    VAR("dR_l1MET", ROOT::Math::VectorUtil::DeltaR(t1, met));
    VAR("dR_l2MET", ROOT::Math::VectorUtil::DeltaR(t2, met));
    VAR("dR_l1l2MET", ROOT::Math::VectorUtil::DeltaR(t1+t2, met));
    VAR("dR_METhtautau", ROOT::Math::VectorUtil::DeltaR(Htt+met, met));
    VAR("dR_METhtautau_sv", SVfit_is_valid ? ROOT::Math::VectorUtil::DeltaR(Htt_sv, met) : default_value);
    VAR("dR_hbbMET", ROOT::Math::VectorUtil::DeltaR(Hbb, met));
    VAR("dR_hbbhtautau", ROOT::Math::VectorUtil::DeltaR(Hbb, Htt+met));
    VAR("dR_hbbhtautau_sv", SVfit_is_valid ? ROOT::Math::VectorUtil::DeltaR(Hbb, Htt_sv) : default_value);
    VAR("dR_b1b2Pt_hbb", (ROOT::Math::VectorUtil::DeltaR(b1, b2))*Hbb.Pt());
    VAR("dR_l1l2Pt_htautau", ROOT::Math::VectorUtil::DeltaR(t1, t2)*(Htt+met).Pt());
    VAR("dR_l1l2Pt_htautau_sv", SVfit_is_valid ? ROOT::Math::VectorUtil::DeltaR(t1, t2)*Htt_sv.Pt() : default_value);
    VAR("dR_b1b2_boosted", four_bodies::Calculate_dR_boosted(b1, b2, Hbb));
    VAR("dR_l1l2_boosted", four_bodies::Calculate_dR_boosted(t1, t2, Htt+met));
    VAR("dR_l1l2_boosted_sv", SVfit_is_valid ? four_bodies::Calculate_dR_boosted(t1, t2, Htt_sv) : default_value);

    VAR("mass_l1l2MET", InvariantMass(t1+t2, met));
    VAR("mass_htautau_vis", Htt.M());
    VAR("mass_htautau", (Htt+met).M());
    VAR("mass_htautau_sv", SVfit_is_valid ? Htt_sv.M() : default_value);
    VAR("mass_l1l2", (t1+t2).M());
    VAR("mass_hbb", Hbb.M());
    VAR("MT_l1", Calculate_MT(t1,met));
    VAR("MT_l2", Calculate_MT(t2,met));
    VAR("MT_htautau", Calculate_MT(Htt+met, met));
    VAR("MT_htautau_sv", SVfit_is_valid ? Calculate_MT(Htt_sv, met) : default_value);
    VAR("MT_l1l2", Calculate_MT(t1+t2, met));
    VAR("MT_tot", Calculate_TotalMT(t1, t2,met)); //Total transverse mass
    VAR("MT2", eventbase.GetMT2()); //Stransverse mass
    VAR("mass_top1", four_bodies::Calculate_topPairMasses(t1, t2, b1, b2, met).first);
    VAR("mass_top2", four_bodies::Calculate_topPairMasses(t1, t2, b1, b2, met).second);
    VAR("mass_X", four_bodies::Calculate_MX(t1, t2, b1, b2, met));
    VAR("mass_H", InvariantMass(Hbb, Htt+met));
    VAR("mass_H_sv", SVfit_is_valid ? InvariantMass(Hbb, Htt_sv) : default_value);
    VAR("mass_H_vis", InvariantMass(Hbb, t1+t2));
    VAR("mass_H_kinfit", kinfit_is_valid ? kinfit_results.mass : default_value);
    VAR("mass_H_kinfit_chi2", kinfit_is_valid ? kinfit_results.chi2 : default_value);

    VAR("phi", four_bodies::Calculate_phi(t1, t2, b1, b2, Htt+met, Hbb));
    VAR("phi_sv", SVfit_is_valid ? four_bodies::Calculate_phi(t1, t2, b1, b2, Htt_sv, Hbb) : default_value);
    VAR("phi_1", four_bodies::Calculate_phi1(t1, t2, Htt+met, Hbb));
    VAR("phi_1_sv", SVfit_is_valid ? four_bodies::Calculate_phi1(t1, t2, Htt_sv, Hbb) : default_value);
    VAR("phi_2", four_bodies::Calculate_phi1(b1, b2, Htt+met, Hbb));
    VAR("phi_2_sv", SVfit_is_valid ? four_bodies::Calculate_phi1(b1, b2, Htt_sv, Hbb) : default_value);
    VAR("costheta_star_leptons", four_bodies::Calculate_cosThetaStar(Htt+met, Hbb));
    VAR("costheta_star_leptons_sv", SVfit_is_valid ? four_bodies::Calculate_cosThetaStar(Htt_sv, Hbb)
            : default_value);
    VAR("costheta_l1htautau", four_bodies::Calculate_cosTheta_2bodies(t1, Htt+met));
    VAR("costheta_l1htautau_sv", SVfit_is_valid ? four_bodies::Calculate_cosTheta_2bodies(t1, Htt_sv)
            : default_value);
    VAR("costheta_l2htautau", four_bodies::Calculate_cosTheta_2bodies(t2, Htt+met));
    VAR("costheta_l2htautau_sv", SVfit_is_valid ? four_bodies::Calculate_cosTheta_2bodies(t2, Htt_sv)
            : default_value);
    VAR("costheta_METhtautau", four_bodies::Calculate_cosTheta_2bodies(met, Htt+met));
    VAR("costheta_METhtautau_sv", SVfit_is_valid ? four_bodies::Calculate_cosTheta_2bodies(met, Htt_sv)
            : default_value);
    VAR("costheta_METhbb", four_bodies::Calculate_cosTheta_2bodies(met, Hbb));
    VAR("costheta_b1hbb", four_bodies::Calculate_cosTheta_2bodies(b1, Hbb));

    VAR("costheta_hbbhh",
            four_bodies::Calculate_cosTheta_2bodies(Hbb, eventbase.GetResonanceMomentum(false,false)));
    VAR("costheta_hbbhh_sv", SVfit_is_valid
            ? four_bodies::Calculate_cosTheta_2bodies(Hbb, eventbase.GetResonanceMomentum(true,false))
            : default_value);
    VAR("costheta_hbbhhMET",
            four_bodies::Calculate_cosTheta_2bodies(Hbb, eventbase.GetResonanceMomentum(false,true)));

    VAR("costheta_htautauhh",
            four_bodies::Calculate_cosTheta_2bodies(Htt+met, eventbase.GetResonanceMomentum(false,false)));
    VAR("costheta_htautauhh_sv", SVfit_is_valid
            ? four_bodies::Calculate_cosTheta_2bodies(Htt+met, eventbase.GetResonanceMomentum(true,false))
            : default_value);
    VAR("costheta_htautauhhMET",
            four_bodies::Calculate_cosTheta_2bodies(Htt+met, eventbase.GetResonanceMomentum(false,true)));
    VAR("costheta_htautau_svhh", SVfit_is_valid
            ? four_bodies::Calculate_cosTheta_2bodies(Htt_sv, eventbase.GetResonanceMomentum(false,false))
            : default_value);
    VAR("costheta_htautau_svhh_sv", SVfit_is_valid
            ? four_bodies::Calculate_cosTheta_2bodies(Htt_sv, eventbase.GetResonanceMomentum(true,false))
            : default_value);
    VAR("costheta_htautau_svhhMET", SVfit_is_valid
            ? four_bodies::Calculate_cosTheta_2bodies(Htt_sv, eventbase.GetResonanceMomentum(false,true))
            : default_value);

    VAR("costheta_l1l2METhh",
            four_bodies::Calculate_cosTheta_2bodies(t1+t2+met, eventbase.GetResonanceMomentum(false,false)));
    VAR("costheta_l1l2METhh_sv", SVfit_is_valid
            ? four_bodies::Calculate_cosTheta_2bodies(t1+t2+met, eventbase.GetResonanceMomentum(true,false))
            : default_value);
    VAR("costheta_l1l2METhhMET",
            four_bodies::Calculate_cosTheta_2bodies(t1+t2+met, eventbase.GetResonanceMomentum(false,true)));

    VAR("mass", mass.mass);
    VAR_INT("channel", eventbase->channelId);
    VAR_INT("spin", spin);
    VAR("kl", spin);

    size_t test = which_test ==-1 ? which_set(gen) : static_cast<size_t>(which_test);
    AddEventVariables(test, mass, eventbase->weight_total, sample_weight, spin,
                      ToString(static_cast<Channel>(eventbase->channelId)));
}

#undef VAR
#undef VAR_INT

}
}
