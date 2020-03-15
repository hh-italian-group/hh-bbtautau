/*! Definition of MvaReader.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "TMVA/Reader.h"
#include "MvaVariables.h"

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
                           const VarNameSet& _enabled_vars);

    virtual void SetValue(const std::string& name, double value, char /*type*/) override;
    virtual void AddEventVariables(size_t /*istraining*/, const SampleId& /*mass*/, double /*weight*/,
                                   double /*sampleweight*/, int /*spin*/, std::string /*channel*/) override {}

    virtual double Evaluate() override;
    virtual std::shared_ptr<TMVA::Reader> GetReader() override;
};

class LegacyMvaVariables : public MvaVariablesBase {
private:
    float dphi_mumet, dphi_metsv, dR_bb, dR_bbbb, dR_taumu, dR_taumusvfit, mT1, mT2, dphi_bbmet, dphi_bbsv;
    std::string method_name;
    bool isLow;
    std::shared_ptr<TMVA::Reader> reader;

public:
    LegacyMvaVariables(const std::string& _method_name, const std::string& bdt_weights, bool _isLow);
    virtual void AddEvent(analysis::EventInfo& eventbase, const SampleId& /*mass*/ , int /* spin*/,
                          double /*sample_weight*/, int /*which_test*/) override;
    virtual double Evaluate() override;
    virtual std::shared_ptr<TMVA::Reader> GetReader() override;
};

class MvaReader {
public:
    using Vars = MvaVariablesBase;
    using VarsPtr = std::shared_ptr<Vars>;

    struct MvaKey {
        std::string method_name;
        int mass, spin;

        bool operator<(const MvaKey& other) const;
    };

    using MethodMap = std::map<MvaKey, VarsPtr>;

    VarsPtr Add(const MvaKey& key, const std::string& bdt_weights, const std::unordered_set<std::string>& enabled_vars,
                bool is_legacy = false, bool is_Low = true);
    double Evaluate(const MvaKey& key, EventInfo* event);

private:
    VarsPtr CreateMvaVariables(const std::string& method_name, const std::string&  bdt_weights,
                               const std::unordered_set<std::string>& enabled_vars, bool is_legacy, bool isLow = true);

private:
    MethodMap methods;
};

}
}
