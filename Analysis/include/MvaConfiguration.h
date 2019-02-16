/*! Definition of MvaVariablesStudy class, the main class for Mva studies.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "hh-bbtautau/Analysis/include/MvaVariables.h"

namespace analysis {
namespace mva_study{

using ::analysis::operator<<;
using ::analysis::operator>>;

struct MvaOptionBase {

    MvaOptionBase(const std::string& _name, bool _isSignificant,
                  const std::set<std::string>& _excludedBoostTypes);

    virtual ~MvaOptionBase() {}

    bool IsEnabled(const std::string& boostType) const;
    bool IsSignificant() const;
    const std::string& Name() const;
    bool IsRange() const;

    virtual std::string GetNameSuffix(size_t n) const = 0;
    virtual std::string GetConfigString(size_t n) const = 0;
    virtual size_t GetNumberOfRangeEntries() const = 0;
    virtual double GetNumericValue(size_t n) const;

private:
    std::string name;
    bool isSignificant;
    std::set<std::string> excludedBoostTypes;
};

template<typename T>
struct MvaNumericOption : MvaOptionBase {
    using Value = T;

    MvaNumericOption(const std::string& _name, const Value& _value, bool _isSignificant,
              const std::set<std::string>& _excludedBoostTypes)
        : MvaOptionBase(_name, _isSignificant, _excludedBoostTypes), isRange(false),
          value(_value)
    { }

    MvaNumericOption(const std::string& _name, const RangeWithStep<Value>& _range, bool _isSignificant,
              const std::set<std::string>& _excludedBoostTypes)
        : MvaOptionBase(_name, _isSignificant, _excludedBoostTypes), isRange(true),
          range(_range)
    { }

    virtual std::string GetNameSuffix(size_t n) const override
    {
        return ToString(n, "_");
    }

    virtual std::string GetConfigString(size_t n) const override
    {
        return ToString(n, "=");
    }

    virtual size_t GetNumberOfRangeEntries() const override
    {
        return isRange ? range.n_grid_points() : 1;
    }

    virtual double GetNumericValue(size_t n) const override { return GetPoint(n); }

protected:

    Value GetPoint(size_t n) const
    {
        if(n >= GetNumberOfRangeEntries())
            throw exception("Option value index is out of range.");
        return isRange ? range.grid_point_value(n) : value;
    }

    std::string ToString(size_t n, const std::string& sep) const
    {
        std::ostringstream ss;
        ss << std::boolalpha << Name() << sep << GetPoint(n);
        return ss.str();
    }

private:
    bool isRange;
    Value value;
    RangeWithStep<Value> range;
};

template<typename T>
struct MvaOption : MvaNumericOption<T> {
    using MvaNumericOption<T>::MvaNumericOption;
};

template<>
struct MvaOption<bool> : MvaNumericOption<bool> {
    using MvaNumericOption<bool>::MvaNumericOption;

    virtual std::string GetConfigString(size_t n) const override;
};

template<>
struct MvaOption<std::string> : MvaOptionBase {
    using Value = std::string;

    MvaOption(const std::string& _name, const Value& _value, bool _isSignificant,
              const std::set<std::string>& _excludedBoostTypes);

    MvaOption(const std::string& _name, const std::vector<Value>& _values, bool _isSignificant,
              const std::set<std::string>& _excludedBoostTypes);

    virtual std::string GetNameSuffix(size_t n) const override;
    virtual std::string GetConfigString(size_t n) const override;
    virtual size_t GetNumberOfRangeEntries() const override;

protected:
    const Value GetPoint(size_t n) const;
    std::string ToString(size_t n, const std::string& sep) const;

private:
    std::vector<Value> values;
};

struct MvaOptionCollection {
    void AddOption(std::shared_ptr<MvaOptionBase> opt);
    std::string GetName(const std::vector<size_t>& pos) const;
    std::string GetConfigString(const std::vector<size_t>& pos) const;
    std::vector<size_t> GetPositionLimits();
    std::shared_ptr<const MvaOptionBase> at(const std::string& name) const;
    std::shared_ptr<const MvaOptionBase> AtIndex(size_t n) const;
    size_t GetIndex(const std::string& name) const;
    size_t count(const std::string& name) const;

    double GetNumericValue(const std::vector<size_t>& point, const std::string& name) const;
    std::map<std::string, size_t> GetOptionNames() const;

private:
    std::vector<std::shared_ptr<MvaOptionBase>> options;
    std::map<std::string, size_t> option_names;
};

struct SampleEntry{
    std::string filename;
    double weight{-1};
    SampleId id;
    int spin;
};

std::ostream& operator<<(std::ostream& os, const SampleEntry& entry);
std::istream& operator>>(std::istream& is, SampleEntry& entry);

using SampleEntryCollection = std::vector<SampleEntry>;

struct SampleEntryList {
    std::string name;
    SampleEntryCollection files;
};

struct MvaSetup {
    std::string name;
    std::vector<Channel> channels;
    Range<int> mass_range;
    std::map<std::string, std::string> param_list;
    std::map<std::string, std::string> param_range;
    std::map<std::string, std::string> params;
    std::map<std::string, std::string> disabled_params;
    std::vector<std::string> significant_params;
    std::vector<std::string> variables;
    bool use_mass_var;

    MvaSetup();

    std::vector<std::string> GetParamList(const std::string& name) const;
    MvaOptionCollection CreateOptionCollection(bool run_on_grid) const;

private:
    std::map<std::string, std::set<std::string>> InvertDisabledParams() const;
};

using MvaSetupCollection = std::unordered_map<std::string, MvaSetup>;
using SampleEntryListCollection = std::unordered_map<std::string, SampleEntryList>;

static const Range<int> low_mass(250, 320), medium_mass(340, 400), high_mass(450,900);
static const std::vector<Range<int>> ranges{low_mass, medium_mass, high_mass};

//static const Range<int> all_mass(250, 900);
//static const std::vector<Range<int>> ranges{all_mass};

}
}
