/*! Definition of classes that describe TupleSkimmer configuration.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <limits>
#include "AnalysisTools/Core/include/Tools.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/McCorrections/include/WeightingMode.h"
#include "hh-bbtautau/Analysis/include/AnalysisCategories.h"

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
    unsigned n_splits{0};
    unsigned split_seed{0};

    std::set<uint32_t> tau_id_hashes;

    //light setup
    bool apply_mass_cut{false};
    bool apply_charge_cut{false};
    bool keep_genJets{false}, keep_genParticles{false}, keep_MET_cov{false};
    std::map<std::string, double> tau_id_cut;
    std::map<uint32_t, double> tau_id_cut_hashes;

    std::map<SelectionCut,analysis::EllipseParameters> massWindowParams;

    void UpdateTauIdHashes()
    {
        tau_id_hashes.clear();
        for(const auto& id : tau_ids)
            tau_id_hashes.insert(tools::hash(id));
        tau_id_cut_hashes.clear();
        for(const auto& id_cut: tau_id_cut)
            tau_id_cut_hashes[tools::hash(id_cut.first)] = id_cut.second;
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
    bool apply_dm_fix{true};
    mc_corrections::WeightingMode weights;

    bool ProduceMergedOutput() const { return merged_output.size(); }
};

using SkimJobCollection = std::unordered_map<std::string, SkimJob>;

} // namespace tuple_skimmer
} // namespace analysis
