/*! Definition of classes that describe TupleSkimmer configuration.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Instruments/include/SkimmerConfig.h"

namespace analysis {
namespace tuple_skimmer {

bool Setup::ApplyTauIdCuts() const { return !tau_id_cuts.empty(); }

FileDescriptor::FileDescriptor(const std::vector<std::string>& _inputs, const std::string& _output,
                               const std::string& _cross_section) :
    output(_output), cross_section(_cross_section), first_input_is_ref(false)
{
    for(const auto& entry: _inputs){
        std::size_t pos = entry.find(":");
        std::string partial = entry.substr(0,pos);
        if (partial == "part") input_is_partial.emplace_back(true);
        else input_is_partial.emplace_back(false);
        std::string str = pos == std::string::npos ? entry : entry.substr(pos+1);
        inputs.emplace_back(str);
    }
}

bool FileDescriptor::HasCrossSection() const { return !cross_section.empty(); }

double FileDescriptor::GetCrossSection(const CrossSectionProvider& xs_provider) const
{
    if(!HasCrossSection()) throw exception("Not provided cross-section");
    return xs_provider.GetCrossSection(cross_section);
}

bool SkimJob::ProduceMergedOutput() const { return merged_output.size(); }

double SkimJob::GetCrossSection(const CrossSectionProvider& xs_provider) const
{
    if(!cross_section.empty()) return xs_provider.GetCrossSection(cross_section);
    if(merged_output.empty()) return -1;
    double xs = 0;
    for(auto file : files)
        xs += file.GetCrossSection(xs_provider);

    return xs;
}

CrossSectionProvider::CrossSectionProvider(const std::string& cross_section_cfg)
{
    xs_items.Parse(cross_section_cfg);
}

double CrossSectionProvider::GetCrossSection(const std::string& process) const
{
    const auto& items = xs_items.GetItems();
    if(!items.count(process)) throw exception("Not found process %1% in the cross-section list") %process;
    const auto& item = items.at(process);
    if(!item.Has("xs")) throw exception("Not found cross-section for process %1%") %process ;
    return item.Get<double>("xs");

}
} //tuple_skimmer
} // namespace analysis
