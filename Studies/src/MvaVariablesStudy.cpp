/*! Definition of MvaVariablesStudy class, the main class for Mva studies.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Studies/include/MvaVariablesStudy.h"

namespace analysis {
namespace mva_study {

ChannelSpin::ChannelSpin(std::string _channel, int _spin) : channel(_channel), spin(_spin) {}

bool ChannelSpin::operator<(const ChannelSpin& x) const
{
    if (channel != x.channel) return channel < x.channel;
    return spin < x.spin;
}

void MvaVariablesStudy::SetValue(const std::string& name, double value, char /*type = 'F'*/)
{
    variables[name] = value;
}

void MvaVariablesStudy::AddEventVariables(size_t set, const SampleId& mass,  double /*weight*/,
                                          double /*sampleweight*/, int spin, std::string channel)
{
    ChannelSpin chsp(channel,spin);
    VarData& sample_vars = all_variables[chsp][set][mass];
    for(const auto& name_value : variables) {
        const std::string& name = name_value.first;
        const double value = name_value.second;
        sample_vars[name].push_back(value);
    }
}

const SampleIdVarData& MvaVariablesStudy::GetSampleVariables(std::string channel, int spin, size_t set) const
{
    if(set >= all_variables.size())
        throw exception("Sample part is out of range.");
    ChannelSpin chsp(channel,spin);
    std::cout<<chsp.channel<<"  "<<chsp.spin<<" "<<set<<std::endl;
    return all_variables.at(chsp).at(set);
}

std::shared_ptr<TMVA::Reader> MvaVariablesStudy::GetReader()
{
    throw exception ("GetReader not supported.");
}

Name_ND::Name_ND(const std::string& name) { names.insert(name); }
Name_ND::Name_ND(std::initializer_list<std::string> _names) : names(_names.begin(), _names.end()){}

bool Name_ND::operator<(const Name_ND& x) const
{
    if (names.size() != x.names.size()) return names.size()<x.names.size();
    if (names.size() == 0) return false;
    auto x_iter = x.names.begin();
    for(auto iter = names.begin(); iter != names.end(); ++iter, ++x_iter){
        if(*iter != *x_iter) return *iter < *x_iter;
    }
    return false;
}

Name_ND::const_iterator Name_ND::begin() const { return names.begin(); }
Name_ND::const_iterator Name_ND::end() const { return names.end(); }
size_t Name_ND::size() const { return names.size(); }
size_t Name_ND::count(const std::string& name) const { return names.count(name); }
void Name_ND::insert(const std::string& name) { names.insert(name); }

const std::string& Name_ND::at(size_t n) const
{
    if(n >= size())
        throw exception("Number of name dimensions = %1%, which is less than %2%") % size() % n;
    return *std::next(begin(), static_cast<std::iterator_traits<const_iterator>::difference_type>(n));
}

bool Name_ND::operator==(const Name_ND& other) const
{
    if (names.size() != other.size()) return false;
    auto iter = names.begin(), other_iter = other.names.begin();
    for(; iter != names.end(); ++iter, ++other_iter) {
        if(*iter != *other_iter) return false;
    }
    return true;
}

bool Name_ND::operator !=(const Name_ND& other) const { return !(*this == other); }

std::string Name_ND::ToString(const std::string& sep) const
{
    std::ostringstream ss;
    auto iter = names.begin();
    if(iter != names.end())
        ss << *iter;
    for(iter = std::next(iter); iter != names.end(); ++iter)
        ss << sep << *iter;
    return ss.str();
}

std::ostream& operator<<(std::ostream& os, const Name_ND& name)
{
    if(name.size() > 0){
        os << name.ToString(",");
    }
    return os;
}

}
}
