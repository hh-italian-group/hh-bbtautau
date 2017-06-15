/*! Definition of classes that describe TupleSkimmer configuration.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <limits>
#include "h-tautau/Analysis/include/AnalysisTypes.h"

namespace analysis {
namespace tuple_skimmer {

struct Setup {
    std::string name;
    std::set<EventEnergyScale> energy_scales;
    std::set<Channel> channels;
    std::set<std::string> tau_ids;
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
};

struct SkimJob {
    std::string name;
    std::string dy_weight, sm_weight, tt_weight, wjets_weight;
    std::string merged_output;
    std::vector<FileDescriptor> files;

    bool ProduceMergedOutput() const { return merged_output.size(); }
};

using SkimJobCollection = std::unordered_map<std::string, SkimJob>;

} // namespace tuple_skimmer
} // namespace analysis
