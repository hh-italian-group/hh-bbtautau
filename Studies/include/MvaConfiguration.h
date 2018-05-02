/*! Definition of MvaVariablesStudy class, the main class for Mva studies. 
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "hh-bbtautau/Analysis/include/MvaVariables.h"
#include <iterator>

namespace  analysis {
namespace mva_study{

using ::analysis::operator<<;
using ::analysis::operator>>;

struct MvaOptionBase {

    MvaOptionBase(const std::string& _name, bool _isSignificant,
                  const std::set<std::string>& _excludedBoostTypes) :
        name(_name), isSignificant(_isSignificant), excludedBoostTypes(_excludedBoostTypes)
    { }

    virtual ~MvaOptionBase() {}

    bool IsEnabled(const std::string& boostType) const
    {
        return !excludedBoostTypes.count(boostType);
    }
    bool IsSignificant() const { return isSignificant; }
    const std::string& Name() const { return name; }
    bool IsRange() const { return GetNumberOfRangeEntries() > 1; }


    virtual std::string GetNameSuffix(size_t n) const = 0;
    virtual std::string GetConfigString(size_t n) const = 0;
    virtual size_t GetNumberOfRangeEntries() const = 0;
    virtual double GetNumericValue(size_t n) const { return n; }

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

    virtual std::string GetConfigString(size_t n) const override
    {
        std::ostringstream ss;
        bool value = GetPoint(n);
        if(!value)
            ss << "!";
        ss << Name();
        return ss.str();
    }
};

template<>
struct MvaOption<std::string> : MvaOptionBase {
    using Value = std::string;

    MvaOption(const std::string& _name, const Value& _value, bool _isSignificant,
              const std::set<std::string>& _excludedBoostTypes)
        : MvaOptionBase(_name, _isSignificant, _excludedBoostTypes)
    {
        values.push_back(_value);
    }

    MvaOption(const std::string& _name, const std::vector<Value>& _values, bool _isSignificant,
              const std::set<std::string>& _excludedBoostTypes)
        : MvaOptionBase(_name, _isSignificant, _excludedBoostTypes), values(_values)
    { }

    virtual std::string GetNameSuffix(size_t n) const override
    {
        return GetPoint(n);
    }

    virtual std::string GetConfigString(size_t n) const override
    {
        return ToString(n, "=");
    }

    virtual size_t GetNumberOfRangeEntries() const override
    {
        return values.size();
    }

protected:

    const Value GetPoint(size_t n) const
    {
        if(n >= GetNumberOfRangeEntries())
            throw exception("Option value index is out of range.");
        return  values.at(n);
    }

    std::string ToString(size_t n, const std::string& sep) const
    {
        if(n >= values.size())
            throw exception("Option value index is out of range.");

        std::ostringstream ss;
        ss << std::boolalpha << Name() << sep << GetPoint(n);
        return ss.str();

    }

private:
    std::vector<Value> values;
};

struct MvaOptionCollection {

    void AddOption(std::shared_ptr<MvaOptionBase> opt)
    {
        options.push_back(opt);
        option_names[opt->Name()] = options.size() - 1;
    }

    std::string GetName(const std::vector<size_t>& pos) const
    {
        if(pos.size() != options.size())
            throw exception("Invalid position");
        std::ostringstream ss;
        bool first = true;
        for(size_t n = 0; n < pos.size(); ++n) {
            const auto& opt = options.at(n);
            if(opt->IsSignificant() || opt->IsRange() ) {
                if(!first)
                    ss << "_";
                ss << opt->GetNameSuffix(pos.at(n));
                first = false;
            }
        }
        return ss.str();
    }

    std::string GetConfigString(const std::vector<size_t>& pos) const
    {
        if(pos.size() != options.size())
            throw exception("Invalid position");

        std::ostringstream ss;
        bool first = true;
        for(size_t n = 0; n < pos.size(); ++n) {
            const auto& opt = options.at(n);
            if(!first)
                ss << ":";
            ss << opt->GetConfigString(pos.at(n));
            first = false;
        }
        return ss.str();
    }

    std::vector<size_t> GetPositionLimits()
    {
        std::vector<size_t> pos(options.size());
        for(size_t n = 0; n < pos.size(); ++n) {
            pos.at(n) = options.at(n)->GetNumberOfRangeEntries();
        }
        return pos;
    }

    std::shared_ptr<const MvaOptionBase> at(const std::string& name) const
    {
        return options.at(option_names.at(name));
    }

    std::shared_ptr<const MvaOptionBase> AtIndex(size_t n) const
    {
        return options.at(n);
    }

    size_t GetIndex(const std::string& name) const
    {
        return option_names.at(name);
    }

    size_t count(const std::string& name) const { return option_names.count(name); }

    double GetNumericValue(const std::vector<size_t>& point, const std::string& name) const
    {
        if(!this->count(name)) return -1;
        const auto& option = this->at(name);
        const size_t pos = point.at(this->GetIndex(name));
        return option->GetNumericValue(pos);
    }

    std::map<std::string, size_t> GetOptionNames() const {return option_names;}

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

inline std::ostream& operator<<(std::ostream& os, const SampleEntry& entry)
{
    os << entry.filename << " " << entry.id << " " << entry.weight << " " <<  entry.spin;
    return os;
}

inline std::istream& operator>>(std::istream& is, SampleEntry& entry)
{
    std::string str;
    std::getline(is, str);
    const auto columns = SplitValueList(str, true);
    if(columns.size() != 4)
        throw exception("Invalid sample entry");
    entry.filename = columns.at(0);
    entry.id = Parse<SampleId>(columns.at(1));
    entry.spin = Parse<int>(columns.at(2));
    entry.weight = NumericalExpression(columns.at(3));
    return is;
}

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

    MvaSetup() : use_mass_var(false) {}

    std::vector<std::string> GetParamList(const std::string& name) const
    {
        if(!param_list.count(name))
            throw exception("Param list with name '%1%' not found.") % name;
        return SplitValueList(param_list.at(name), false);
    }

    MvaOptionCollection CreateOptionCollection(bool run_on_grid) const
    {
        MvaOptionCollection mvaOptions;
        auto excludedBoostTypes = InvertDisabledParams();
        std::set<std::string> significant_params_set(significant_params.begin(), significant_params.end());

        for(const auto& param : params) {
            if(run_on_grid && (param_list.count(param.first) || param_range.count(param.first))) continue;
            bool b_value;
            int i_value;
            double d_value;

            std::shared_ptr<MvaOptionBase> opt;
            if(param.second.find('.') != std::string::npos && TryParse(param.second, d_value)) {
                opt = std::make_shared<MvaOption<double>>(param.first, d_value,
                                                          significant_params_set.count(param.first),
                                                          excludedBoostTypes[param.first]);
            } else if(TryParse(param.second, i_value)) {
                opt = std::make_shared<MvaOption<int>>(param.first, i_value,
                                                       significant_params_set.count(param.first),
                                                       excludedBoostTypes[param.first]);
            } else if(TryParse(param.second, b_value)) {
                opt = std::make_shared<MvaOption<bool>>(param.first, b_value,
                                                        significant_params_set.count(param.first),
                                                        excludedBoostTypes[param.first]);
            } else {
                opt = std::make_shared<MvaOption<std::string>>(param.first, param.second,
                                                               significant_params_set.count(param.first),
                                                               excludedBoostTypes[param.first]);
            }
            mvaOptions.AddOption(opt);
        }

        for(const auto& param : param_list) {
            if(mvaOptions.count(param.first)) continue;
            const auto list = GetParamList(param.first);
            auto opt = std::make_shared<MvaOption<std::string>>(
                    param.first, list, significant_params_set.count(param.first),
                    excludedBoostTypes[param.first]);
            mvaOptions.AddOption(opt);
        }

        for(const auto& param : param_range) {
            if(mvaOptions.count(param.first)) continue;
            std::shared_ptr<MvaOptionBase> opt;
            const bool is_int = param.second.find('.') == std::string::npos;
            if(is_int) {
                auto range = Parse<RangeWithStep<int>>(param.second);
                opt = std::make_shared<MvaOption<int>>(param.first, range,
                                                       significant_params_set.count(param.first),
                                                       excludedBoostTypes[param.first]);
            } else {
                auto range = Parse<RangeWithStep<double>>(param.second);
                opt = std::make_shared<MvaOption<double>>(param.first, range,
                                                          significant_params_set.count(param.first),
                                                          excludedBoostTypes[param.first]);
            }
            mvaOptions.AddOption(opt);
        }

        return mvaOptions;
    }
private:
    std::map<std::string, std::set<std::string>> InvertDisabledParams() const
    {
        std::map<std::string, std::set<std::string>> excludedBoosType;
        for(const auto& param : disabled_params){
            std::vector<std::string> option = SplitValueList(param.second, false);
            for (const auto& opt : option){
                excludedBoosType[opt].insert(param.first);
            }
        }
        return excludedBoosType;
    }
};
using MvaSetupCollection = std::unordered_map<std::string, MvaSetup>;
using SampleEntryListCollection = std::unordered_map<std::string, SampleEntryList>;

static const Range<int> low_mass(250, 320), medium_mass(340, 400), high_mass(450,900);
static const std::vector<Range<int>> ranges{low_mass, medium_mass, high_mass};

//static const Range<int> all_mass(250, 900);
//static const std::vector<Range<int>> ranges{all_mass};

}
}
