/*! Definition of classes that describe TupleSkimmer configuration.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Instruments/include/SkimmerConfig.h"

namespace analysis {
namespace tuple_skimmer {

bool Setup::ApplyTauIdCuts() const { return !tau_id_cuts.empty(); }

void Setup::UpdateTauIdIndices()
{
    tau_id_cut_indices.clear();
    const auto& bit_refs = TauIdResults::GetBitRefsByName();
    for(const auto& cut_name : tau_id_cuts) {
        if(!bit_refs.count(cut_name))
            throw exception("Unknown tau ID '%1%'") % cut_name;
        tau_id_cut_indices.insert(bit_refs.at(cut_name));
    }
}

FileDescriptor::FileDescriptor(const std::vector<std::string>& _inputs) :
    inputs(_inputs), output(inputs.size() ? inputs.front() : ""),
    cross_section(std::numeric_limits<double>::quiet_NaN()), first_input_is_ref(true) {}

FileDescriptor::FileDescriptor(const std::vector<std::string>& _inputs, const std::string& _output) :
    inputs(_inputs), output(_output), cross_section(std::numeric_limits<double>::quiet_NaN()),
    first_input_is_ref(false) {}

FileDescriptor::FileDescriptor(const std::vector<std::string>& _inputs, double _cross_section) :
    output(inputs.size() ? inputs.front() : ""), cross_section(_cross_section),
    first_input_is_ref(false)
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

bool FileDescriptor::HasCrossSection() const { return !std::isnan(cross_section); }
double FileDescriptor::GetCrossSectionWeight() const { return HasCrossSection() ? cross_section : 1; }

bool SkimJob::ProduceMergedOutput() const { return merged_output.size(); }

} // namespace tuple_skimmer
} // namespace analysis
