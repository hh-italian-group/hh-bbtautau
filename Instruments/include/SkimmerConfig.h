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
    Period period;
    DiscriminatorWP btag_wp;
    mc_corrections::WeightingMode common_weights;
    unsigned n_splits{0};
    unsigned split_seed{0};

    //light setup
    bool apply_mass_cut{false}, apply_charge_cut{false}, apply_bb_cut{true};
    bool keep_genJets{false}, keep_genParticles{false}, keep_MET_cov{false};
    std::set<std::string> tau_id_cuts;
    std::set<size_t> tau_id_cut_indices;

    std::map<SelectionCut,analysis::EllipseParameters> massWindowParams;

    bool ApplyTauIdCuts() const { return !tau_id_cuts.empty(); }

    void UpdateTauIdIndices()
    {
        tau_id_cut_indices.clear();
        const auto& bit_refs = TauIdResults::GetBitRefsByName();
        for(const auto& cut_name : tau_id_cuts) {
            if(!bit_refs.count(cut_name))
                throw exception("Unknown tau ID '%1%'") % cut_name;
            tau_id_cut_indices.insert(bit_refs.at(cut_name));
        }
    }

};

using SetupCollection = std::unordered_map<std::string, Setup>;

struct FileDescriptor {
    std::vector<std::string> inputs;
    std::string output;
    double cross_section;
    bool first_input_is_ref;
    std::vector<bool> input_is_partial;

    FileDescriptor(const std::vector<std::string>& _inputs = {}) :
        inputs(_inputs), output(inputs.size() ? inputs.front() : ""),
        cross_section(std::numeric_limits<double>::quiet_NaN()), first_input_is_ref(true) {}

    FileDescriptor(const std::vector<std::string>& _inputs, const std::string& _output) :
        inputs(_inputs), output(_output), cross_section(std::numeric_limits<double>::quiet_NaN()),
        first_input_is_ref(false) {}

    FileDescriptor(const std::vector<std::string>& _inputs, double _cross_section) :
        output(inputs.size() ? inputs.front() : ""), cross_section(_cross_section),
        first_input_is_ref(false) {
            unsigned n = 0;
            for(const auto& entry: _inputs){
                std::size_t pos = entry.find(":");
                std::string partial = entry.substr(0,pos);
                if (partial == "part:") input_is_partial[n] = true;
                std::string str = entry.substr(pos+1);
                inputs.emplace_back(str);
                n++;
            }
    }

    bool HasCrossSection() const { return !std::isnan(cross_section); }
    double GetCrossSectionWeight() const { return HasCrossSection() ? cross_section : 1; }
};

struct SkimJob {
    std::string name;
    std::string merged_output;
    std::vector<FileDescriptor> files;
    bool apply_common_weights{true};
    mc_corrections::WeightingMode weights;

    bool ProduceMergedOutput() const { return merged_output.size(); }
};

using SkimJobCollection = std::unordered_map<std::string, SkimJob>;

} // namespace tuple_skimmer
} // namespace analysis
