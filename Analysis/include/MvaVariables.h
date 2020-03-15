/*! Definition of MvaVariables.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <random>
#include <TMVA/Reader.h>
#include "h-tautau/Analysis/include/EventInfo.h"

namespace analysis {
namespace mva_study{

using ::analysis::operator<<;
using ::analysis::operator>>;

enum class SampleType { Sgn_Res = 1, Sgn_NonRes = 0, Bkg_TTbar = -1 };

struct SampleId {
    SampleType sampleType;
    int mass;

    SampleId();
    SampleId(SampleType _sampleType, int _mass = 0);

    bool operator<(const SampleId& x) const;
    bool operator ==(const SampleId& x) const;
    bool operator !=(const SampleId& x) const;
    bool IsSignal() const;
    bool IsBackground() const;
    bool IsSM() const;

    static const SampleId& MassTot();
    static const SampleId& Bkg();
    static const SampleId& SM();
};

//static const SampleId Bkg{SampleType::Bkg_TTbar, -1};

ENUM_NAMES(SampleType) = {
    {SampleType::Sgn_Res, "Sgn_Res"},
    {SampleType::Sgn_NonRes, "SM"},
    {SampleType::Bkg_TTbar, "TT"}
};


std::ostream& operator<<(std::ostream& os, const SampleId& id);
std::ostream& operator<<(std::ostream& os, const std::pair<SampleId, SampleId>& id_pair);
std::istream& operator>>(std::istream& is, SampleId& id);

class MvaVariablesBase {
public:
    using Mutex = std::mutex;
    using Lock = std::lock_guard<Mutex>;

    virtual ~MvaVariablesBase() {}
    virtual void AddEvent(EventInfo& eventbase, const SampleId& mass, int  spin, double sample_weight = 1.,
                          int which_test = -1) = 0;
    virtual double Evaluate();
    virtual std::shared_ptr<TMVA::Reader> GetReader() = 0;

    double AddAndEvaluate(EventInfo& eventbase, const SampleId& mass, int spin,
                          double sample_weight = 1., int which_test = -1);

private:
    Mutex mutex;
};

class MvaVariables : public MvaVariablesBase {
public:
    using VarNameSet = std::unordered_set<std::string>;
    MvaVariables(size_t _number_set = 1, uint_fast32_t seed = std::numeric_limits<uint_fast32_t>::max(),
                 const VarNameSet& _enabled_vars = {}, const VarNameSet& _disabled_vars = {});

    virtual ~MvaVariables() override {}
    virtual void SetValue(const std::string& name, double value, char type = 'F') = 0;
    virtual void AddEventVariables(size_t which_set, const SampleId& mass, double weight, double sampleweight,
                                   int spin, std::string channel) = 0;
    const std::unordered_set<std::string>& GetDisabledVars() const;
    bool IsEnabled(const std::string& name) const;

    virtual void AddEvent(analysis::EventInfo& eventbase, const SampleId& mass, int spin,
                          double sample_weight = 1., int which_test = -1) override;

private:
    std::mt19937_64 gen;
    std::uniform_int_distribution<size_t> which_set;
    VarNameSet enabled_vars, disabled_vars;
};

}
}
