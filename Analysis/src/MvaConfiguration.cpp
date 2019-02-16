/*! Definition of MvaVariablesStudy class, the main class for Mva studies.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/MvaConfiguration.h"

namespace  analysis {
namespace mva_study{

MvaOptionBase::MvaOptionBase(const std::string& _name, bool _isSignificant,
              const std::set<std::string>& _excludedBoostTypes) :
    name(_name), isSignificant(_isSignificant), excludedBoostTypes(_excludedBoostTypes)
{
}

bool MvaOptionBase::IsEnabled(const std::string& boostType) const
{
    return !excludedBoostTypes.count(boostType);
}

bool MvaOptionBase::IsSignificant() const { return isSignificant; }
const std::string& MvaOptionBase::Name() const { return name; }
bool MvaOptionBase::IsRange() const { return GetNumberOfRangeEntries() > 1; }

double MvaOptionBase::GetNumericValue(size_t n) const { return n; }



std::string MvaOption<bool>::GetConfigString(size_t n) const
{
    std::ostringstream ss;
    bool value = GetPoint(n);
    if(!value)
        ss << "!";
    ss << Name();
    return ss.str();
}

MvaOption<std::string>::MvaOption(const std::string& _name, const Value& _value, bool _isSignificant,
                                  const std::set<std::string>& _excludedBoostTypes)
    : MvaOptionBase(_name, _isSignificant, _excludedBoostTypes)
{
    values.push_back(_value);
}

MvaOption<std::string>::MvaOption(const std::string& _name, const std::vector<Value>& _values, bool _isSignificant,
          const std::set<std::string>& _excludedBoostTypes)
    : MvaOptionBase(_name, _isSignificant, _excludedBoostTypes), values(_values)
{ }

std::string MvaOption<std::string>::GetNameSuffix(size_t n) const
{
    return GetPoint(n);
}

std::string MvaOption<std::string>::GetConfigString(size_t n) const
{
    return ToString(n, "=");
}

size_t MvaOption<std::string>::GetNumberOfRangeEntries() const
{
    return values.size();
}

const MvaOption<std::string>::Value MvaOption<std::string>::GetPoint(size_t n) const
{
    if(n >= GetNumberOfRangeEntries())
        throw exception("Option value index is out of range.");
    return  values.at(n);
}

std::string MvaOption<std::string>::ToString(size_t n, const std::string& sep) const
{
    if(n >= values.size())
        throw exception("Option value index is out of range.");

    std::ostringstream ss;
    ss << std::boolalpha << Name() << sep << GetPoint(n);
    return ss.str();

}

void MvaOptionCollection::AddOption(std::shared_ptr<MvaOptionBase> opt)
{
    options.push_back(opt);
    option_names[opt->Name()] = options.size() - 1;
}

std::string MvaOptionCollection::GetName(const std::vector<size_t>& pos) const
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

std::string MvaOptionCollection::GetConfigString(const std::vector<size_t>& pos) const
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

std::vector<size_t> MvaOptionCollection::GetPositionLimits()
{
    std::vector<size_t> pos(options.size());
    for(size_t n = 0; n < pos.size(); ++n) {
        pos.at(n) = options.at(n)->GetNumberOfRangeEntries();
    }
    return pos;
}

std::shared_ptr<const MvaOptionBase> MvaOptionCollection::at(const std::string& name) const
{
    return options.at(option_names.at(name));
}

std::shared_ptr<const MvaOptionBase> MvaOptionCollection::AtIndex(size_t n) const
{
    return options.at(n);
}

size_t MvaOptionCollection::GetIndex(const std::string& name) const
{
    return option_names.at(name);
}

size_t MvaOptionCollection::count(const std::string& name) const { return option_names.count(name); }

double MvaOptionCollection::GetNumericValue(const std::vector<size_t>& point, const std::string& name) const
{
    if(!this->count(name)) return -1;
    const auto& option = this->at(name);
    const size_t pos = point.at(this->GetIndex(name));
    return option->GetNumericValue(pos);
}

std::map<std::string, size_t> MvaOptionCollection::GetOptionNames() const {return option_names;}

std::ostream& operator<<(std::ostream& os, const SampleEntry& entry)
{
    os << entry.filename << " " << entry.id << " " << entry.weight << " " <<  entry.spin;
    return os;
}

std::istream& operator>>(std::istream& is, SampleEntry& entry)
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


MvaSetup::MvaSetup() : use_mass_var(false) {}

std::vector<std::string> MvaSetup::GetParamList(const std::string& name) const
{
    if(!param_list.count(name))
        throw exception("Param list with name '%1%' not found.") % name;
    return SplitValueList(param_list.at(name), false);
}

MvaOptionCollection MvaSetup::CreateOptionCollection(bool run_on_grid) const
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
std::map<std::string, std::set<std::string>> MvaSetup::InvertDisabledParams() const
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

}
}
