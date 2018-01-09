/*! Definition of MvaReader.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "TMVA/Reader.h"
#include "MvaVariables.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"



namespace analysis {
namespace mva_study{

class MvaVariablesEvaluation : public MvaVariables {
private:
    using DataVectorF = std::vector<float>;
    static constexpr size_t max_n_vars = 1000;
    DataVectorF variable_float;
    std::map<std::string, size_t> name_indices;
    std::string method_name, bdt_weights;
    std::shared_ptr<TMVA::Reader> reader;
    size_t n_vars;
    bool is_initialized;

public:
    MvaVariablesEvaluation(const std::string& _method_name, const std::string& _bdt_weights,
                           const VarNameSet& _enabled_vars) :
        MvaVariables(1, 0, _enabled_vars), variable_float(max_n_vars), method_name(_method_name), bdt_weights(_bdt_weights),
        reader(new TMVA::Reader), n_vars(0), is_initialized(false)
    {}

    virtual void SetValue(const std::string& name, double value, char /*type*/) override
    {
        if (!name_indices.count(name)){
            ++n_vars;
            name_indices[name] = n_vars - 1;
            reader->AddVariable(name, &variable_float.at(n_vars - 1));
        }
        variable_float.at(name_indices.at(name)) = static_cast<float>(value);
    }

    virtual void AddEventVariables(size_t /*istraining*/, const SampleId& /*mass*/, double /*weight*/, double /*sampleweight*/, int /*spin*/, std::string /*channel*/) override {}

    virtual double Evaluate() override
    {
        if(!is_initialized) {
            reader->BookMVA(method_name, bdt_weights);
            is_initialized = true;
        }
        return reader->EvaluateMVA(method_name);
    }

    virtual std::shared_ptr<TMVA::Reader> GetReader() override { return reader;}
};

class LegacyMvaVariables : public MvaVariablesBase {
private:
    float dphi_mumet, dphi_metsv, dR_bb, dR_bbbb, dR_taumu, dR_taumusvfit, mT1, mT2, dphi_bbmet, dphi_bbsv;
    std::string method_name;
    bool isLow;
    std::shared_ptr<TMVA::Reader> reader;

public:
    LegacyMvaVariables(const std::string& _method_name, const std::string& bdt_weights, bool _isLow)
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

    virtual void AddEvent(analysis::EventInfoBase& eventbase, const SampleId& /*mass*/ , int /* spin*/, double /*sample_weight*/, int /*which_test*/, double /*weight_bkg*/) override
    {
        const auto& Htt = eventbase.GetHiggsTTMomentum(false);
        const auto& Htt_sv = eventbase.GetHiggsTTMomentum(true);
        const auto& t1 = eventbase.GetLeg(1);
        const auto& t2 = eventbase.GetLeg(2);

        const auto& Hbb = eventbase.GetHiggsBB();
        const auto& b1 = Hbb.GetFirstDaughter();
        const auto& b2 = Hbb.GetSecondDaughter();

        const auto& met = eventbase.GetMET();

        dphi_mumet = std::abs(ROOT::Math::VectorUtil::DeltaPhi(t1.GetMomentum(), met.GetMomentum()));
        dphi_metsv = std::abs(ROOT::Math::VectorUtil::DeltaPhi(Htt_sv, met.GetMomentum()));
        dR_bb = std::abs(ROOT::Math::VectorUtil::DeltaR(b1.GetMomentum(), b2.GetMomentum()));
        dR_bbbb = ROOT::Math::VectorUtil::DeltaR(b1.GetMomentum(), b2.GetMomentum())*Hbb.GetMomentum().Pt();
        dR_taumu = std::abs(ROOT::Math::VectorUtil::DeltaR(t1.GetMomentum(), t2.GetMomentum()));
        dR_taumusvfit = ROOT::Math::VectorUtil::DeltaR(t1.GetMomentum(), t2.GetMomentum())*Htt.Pt();
        mT1 = static_cast<float>(Calculate_MT(t1.GetMomentum(), met.GetMomentum()));
        mT2 = static_cast<float>(Calculate_MT(t2.GetMomentum(), met.GetMomentum()));
        dphi_bbmet = std::abs(ROOT::Math::VectorUtil::DeltaPhi(Hbb.GetMomentum(), met.GetMomentum()));
        dphi_bbsv = std::abs(ROOT::Math::VectorUtil::DeltaPhi(Hbb.GetMomentum(), Htt_sv));
    }
    virtual double Evaluate() override { return reader->EvaluateMVA(method_name); }

    virtual std::shared_ptr<TMVA::Reader> GetReader() override { return reader;}

};

class MvaReader {
public:
    using Vars = MvaVariablesBase;
    using VarsPtr = std::shared_ptr<Vars>;
    using Range = ::analysis::Range<int>;

    struct VarRange { Range range; VarsPtr vars; };
    using VarRangeMap = std::map<int, VarRange>;
    using MethodRanges = std::map<std::string, VarRangeMap>;

    VarsPtr AddRange(const Range& mass_range, const std::string& method_name, const std::string& bdt_weights,
                  const std::unordered_set<std::string>& enabled_vars, bool is_legacy = false, bool is_Low = true)
    {
        if(method_vars.count(method_name))
            throw exception("MVA method with name '%1%' already exists.") % method_name;

        auto& vars = method_vars[method_name];
        auto bound = vars.lower_bound(mass_range.min());
        for(size_t n = 0; n < 2 && bound != vars.end() && bound->second.range != mass_range; ++n, ++bound) {
            if(bound->second.range.Overlaps(mass_range, false))
                throw exception("Overlapping ranges.");
        }

        vars[mass_range.max()].range = mass_range;
        return vars[mass_range.max()].vars = CreateMvaVariables(method_name, bdt_weights, enabled_vars, is_legacy,
                                                                is_Low);
    }

    double Evaluate(EventInfoBase* event, int mass, const std::string& method_name, int spin)
    {
        auto& vars = FindMvaVariables(mass, method_name);
        return vars.AddAndEvaluate(*event, SampleId(SampleType::Sgn_Res, mass), spin);
    }

    Vars& FindMvaVariables(int mass, const std::string& method_name = "")
    {
        if(!method_vars.count(method_name))
            throw exception("Method '%1%' not found.") % method_name;
        const auto& vars = method_vars.at(method_name);
        auto low_bound = vars.lower_bound(mass);
        if(low_bound == vars.end() || !low_bound->second.range.Contains(mass))
            throw exception("Range not found for method = '%1%' for mass = %2%.") % method_name % mass;
        return *low_bound->second.vars;
    }

private:
    VarsPtr CreateMvaVariables(const std::string& method_name, const std::string&  bdt_weights,
                               const std::unordered_set<std::string>& enabled_vars, bool is_legacy, bool isLow = true)
    {
        if (is_legacy) return std::make_shared<LegacyMvaVariables>(method_name, bdt_weights, isLow);
        return std::make_shared<MvaVariablesEvaluation>(method_name, bdt_weights, enabled_vars);
    }

private:
    MethodRanges method_vars;
};

}
}
