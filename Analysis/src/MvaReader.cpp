/*! Definition of MvaReader.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/MvaReader.h"

namespace analysis {
namespace mva_study{

MvaVariablesEvaluation::MvaVariablesEvaluation(const std::string& _method_name, const std::string& _bdt_weights,
                                               const VarNameSet& _enabled_vars) :
    MvaVariables(1, 0, _enabled_vars), variable_float(max_n_vars), method_name(_method_name), bdt_weights(_bdt_weights),
    reader(new TMVA::Reader), n_vars(0), is_initialized(false)
{
}

void MvaVariablesEvaluation::SetValue(const std::string& name, double value, char /*type*/)
{
    if (!name_indices.count(name)){
        ++n_vars;
        name_indices[name] = n_vars - 1;
        reader->AddVariable(name, &variable_float.at(n_vars - 1));
    }
    variable_float.at(name_indices.at(name)) = static_cast<float>(value);
}

double MvaVariablesEvaluation::Evaluate()
{
    if(!is_initialized) {
        reader->BookMVA(method_name, bdt_weights);
        is_initialized = true;
    }
    return reader->EvaluateMVA(method_name);
}

std::shared_ptr<TMVA::Reader> MvaVariablesEvaluation::GetReader() { return reader; }

LegacyMvaVariables::LegacyMvaVariables(const std::string& _method_name, const std::string& bdt_weights, bool _isLow)
    : method_name(_method_name), isLow(_isLow), reader(new TMVA::Reader)
{
    reader->AddVariable("dphi_mumet", &dphi_mumet);
    reader->AddVariable("dphi_metsv", &dphi_metsv);
    if (isLow){
        reader->AddVariable("dR_bb*bb_h.Pt()", &dR_bbbb);
        reader->AddVariable("dR_taumu*SVfit_p4.Pt()", &dR_taumusvfit);
    }
    else {
        reader->AddVariable("dR_bb", &dR_bb);
        reader->AddVariable("dR_taumu", &dR_taumu);
    }
    reader->AddVariable("mT1", &mT1);
    reader->AddVariable("mT2", &mT2);
    reader->AddVariable("dphi_bbmet", &dphi_bbmet);
    reader->AddVariable("dphi_bbsv", &dphi_bbsv);

    reader->BookMVA(method_name, bdt_weights);
}

void LegacyMvaVariables::AddEvent(analysis::EventInfo& eventbase, const SampleId& /*mass*/ , int /* spin*/,
                                  double /*sample_weight*/, int /*which_test*/)
{
    const auto& Htt = eventbase.GetHiggsTTMomentum(false);
    const auto& Htt_sv = eventbase.GetHiggsTTMomentum(true);
    const auto& t1 = eventbase.GetLeg(1);
    const auto& t2 = eventbase.GetLeg(2);

    const auto& Hbb = eventbase.GetHiggsBB();
    const auto& b1 = Hbb.GetFirstDaughter();
    const auto& b2 = Hbb.GetSecondDaughter();

    const auto& met = eventbase.GetMET();

    dphi_mumet = static_cast<float>(std::abs(ROOT::Math::VectorUtil::DeltaPhi(t1.GetMomentum(), met.GetMomentum())));
    dphi_metsv = static_cast<float>(std::abs(ROOT::Math::VectorUtil::DeltaPhi(Htt_sv, met.GetMomentum())));
    dR_bb = static_cast<float>(std::abs(ROOT::Math::VectorUtil::DeltaR(b1.GetMomentum(), b2.GetMomentum())));
    dR_bbbb = static_cast<float>(ROOT::Math::VectorUtil::DeltaR(b1.GetMomentum(), b2.GetMomentum())
                                 * Hbb.GetMomentum().Pt());
    dR_taumu = static_cast<float>(std::abs(ROOT::Math::VectorUtil::DeltaR(t1.GetMomentum(), t2.GetMomentum())));
    dR_taumusvfit = static_cast<float>(ROOT::Math::VectorUtil::DeltaR(t1.GetMomentum(), t2.GetMomentum())*Htt.Pt());
    mT1 = static_cast<float>(Calculate_MT(t1.GetMomentum(), met.GetMomentum()));
    mT2 = static_cast<float>(Calculate_MT(t2.GetMomentum(), met.GetMomentum()));
    dphi_bbmet = static_cast<float>(std::abs(ROOT::Math::VectorUtil::DeltaPhi(Hbb.GetMomentum(), met.GetMomentum())));
    dphi_bbsv = static_cast<float>(std::abs(ROOT::Math::VectorUtil::DeltaPhi(Hbb.GetMomentum(), Htt_sv)));
}

double LegacyMvaVariables::Evaluate() { return reader->EvaluateMVA(method_name); }
std::shared_ptr<TMVA::Reader> LegacyMvaVariables::GetReader() { return reader; }


bool MvaReader::MvaKey::operator<(const MvaKey& other) const
{
    if(method_name != other.method_name) return method_name < other.method_name;
    if(mass != other.mass) return mass < other.mass;
    return spin < other.spin;
}

MvaReader::VarsPtr MvaReader::Add(const MvaKey& key, const std::string& bdt_weights,
                                  const std::unordered_set<std::string>& enabled_vars, bool is_legacy, bool is_Low)
{
    if(methods.count(key))
        throw exception("MVA method with (name, mass, spin) = '%1%, %2%, %3%' already exists.")
            % key.method_name % key.mass % key.spin;

    return methods[key] = CreateMvaVariables(key.method_name, bdt_weights, enabled_vars, is_legacy, is_Low);
}

double MvaReader::Evaluate(const MvaKey& key, EventInfo* event)
{
    auto iter = methods.find(key);
    if(iter == methods.end())
        throw exception("Method '%1%' not found.") % key.method_name;
    return iter->second->AddAndEvaluate(*event, SampleId(SampleType::Sgn_Res, key.mass), key.spin);
}

MvaReader::VarsPtr MvaReader::CreateMvaVariables(const std::string& method_name, const std::string&  bdt_weights,
                                                 const std::unordered_set<std::string>& enabled_vars, bool is_legacy,
                                                 bool isLow)
{
    if (is_legacy) return std::make_shared<LegacyMvaVariables>(method_name, bdt_weights, isLow);
    return std::make_shared<MvaVariablesEvaluation>(method_name, bdt_weights, enabled_vars);
}

}
}
