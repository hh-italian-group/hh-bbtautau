/*! Definition of classes that describe TupleSkimmer configuration.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <limits>
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/McCorrections/include/WeightingMode.h"

namespace analysis {
namespace tuple_skimmer {

struct Setup {
    std::string name;
    std::set<EventEnergyScale> energy_scales;
    std::set<Channel> channels;
    std::set<std::string> tau_ids;
    Period period;
    DiscriminatorWP btag_wp;
    mc_corrections::WeightingMode common_weights;

    std::set<uint32_t> tau_id_hashes;

    void UpdateTauIdHashes()
    {
        tau_id_hashes.clear();
        for(const auto& id : tau_ids)
            tau_id_hashes.insert(tools::hash(id));
    }
};

using SetupCollection = std::unordered_map<std::string, Setup>;

struct FileDescriptor {
    std::vector<std::string> inputs;
    std::string output;
    double cross_section;
    bool first_input_is_ref;

    FileDescriptor(const std::vector<std::string>& _inputs = {}) :
        inputs(_inputs), output(inputs.size() ? inputs.front() : ""),
        cross_section(std::numeric_limits<double>::quiet_NaN()), first_input_is_ref(true) {}

    FileDescriptor(const std::vector<std::string>& _inputs, const std::string& _output) :
        inputs(_inputs), output(_output), cross_section(std::numeric_limits<double>::quiet_NaN()),
        first_input_is_ref(false) {}

    FileDescriptor(const std::vector<std::string>& _inputs, double _cross_section) :
        inputs(_inputs), output(inputs.size() ? inputs.front() : ""), cross_section(_cross_section),
        first_input_is_ref(false) {}

    bool HasCrossSection() const { return !std::isnan(cross_section); }
    double GetCrossSectionWeight() const { return HasCrossSection() ? cross_section : 1; }
};

struct SkimJob {
    std::string name;
    std::string merged_output;
    std::vector<FileDescriptor> files;
    bool apply_common_weights{true};
    mc_corrections::WeightingMode weights;
    bool ignore_trigger_info{false};

    bool ProduceMergedOutput() const { return merged_output.size(); }
};

using SkimJobCollection = std::unordered_map<std::string, SkimJob>;

} // namespace tuple_skimmer
} // namespace analysis
