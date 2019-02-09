/*! Definition of MvaVariablesStudy class, the main class for Mva studies.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "hh-bbtautau/Analysis/include/MvaVariables.h"

namespace analysis {
namespace mva_study{

using DataVector = std::vector<double>;
using VarData = std::map<std::string, DataVector>;
using SampleIdVarData = std::map<SampleId, VarData>;

struct ChannelSpin{
    std::string channel;
    int spin;

    ChannelSpin(std::string _channel, int _spin);
    bool operator<(const ChannelSpin& x) const;
};

class MvaVariablesStudy : public MvaVariables {
private:
    std::map<std::string, double> variables;
    using SplittedSampleIdVarData = std::map<size_t, SampleIdVarData>;
    std::map<ChannelSpin,SplittedSampleIdVarData> all_variables;

public:
    using MvaVariables::MvaVariables;

    virtual void SetValue(const std::string& name, double value, char /*type = 'F'*/) override;
    virtual void AddEventVariables(size_t set, const SampleId& mass,  double /*weight*/, double /*sampleweight*/,
                                   int spin, std::string channel) override;
    const SampleIdVarData& GetSampleVariables(std::string channel, int spin, size_t set = 0) const;
    virtual std::shared_ptr<TMVA::Reader> GetReader() override;
};

struct Name_ND {
private:
    std::set<std::string> names;

public:
    using const_iterator = std::set<std::string>::const_iterator;

    Name_ND() {}
    Name_ND(const std::string& name);
    Name_ND(std::initializer_list<std::string> _names);

    bool operator<(const Name_ND& x) const;
    const_iterator begin() const;
    const_iterator end() const;
    size_t size() const;
    size_t count(const std::string& name) const;
    void insert(const std::string& name);
    const std::string& at(size_t n) const;

    template<typename Set>
    bool IsSubset(const Set& set) const
    {
        for(const auto& name : names)
            if(!set.count(name)) return false;
        return true;
    }

    bool operator==(const Name_ND& other) const;
    bool operator !=(const Name_ND& other) const;

    std::string ToString(const std::string& sep = "_") const;
};

std::ostream& operator<<(std::ostream& os, const Name_ND& name);

}
}
