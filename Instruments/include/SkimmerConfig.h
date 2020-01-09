/*! Definition of classes that describe TupleSkimmer configuration.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/JetTools/include/BTagger.h"
#include "h-tautau/McCorrections/include/WeightingMode.h"
#include "hh-bbtautau/Analysis/include/AnalysisCategories.h"
#include "h-tautau/Analysis/include/SignalObjectSelector.h"

namespace analysis {
namespace tuple_skimmer {

struct Setup {
    std::string name;
    std::set<UncertaintySource> unc_sources;
    std::set<Channel> channels;
    Period period;
    SignalMode mode;
    DiscriminatorWP btag_wp;
    std::vector<std::string> cachePaths;
    mc_corrections::WeightingMode common_weights;
    JetOrdering jet_ordering;
    unsigned n_splits{0};
    unsigned split_seed{0};

    //light setup
    bool apply_mass_cut{false}, apply_charge_cut{false}, apply_bb_cut{true};
    bool keep_genJets{false}, keep_genParticles{false}, keep_MET_cov{true};
    bool apply_kinfit{true};
    bool applyTauId{true};
    std::set<std::string> tau_id_cuts;
    std::set<size_t> tau_id_cut_indices;

    std::map<SelectionCut,analysis::EllipseParameters> massWindowParams;

    bool ApplyTauIdCuts() const;

};

using SetupCollection = std::unordered_map<std::string, Setup>;

struct FileDescriptor {
    std::vector<std::string> inputs;
    std::string output;
    double cross_section;
    bool first_input_is_ref;
    std::vector<bool> input_is_partial;

    FileDescriptor(const std::vector<std::string>& _inputs = {});
    FileDescriptor(const std::vector<std::string>& _inputs, const std::string& _output);
    FileDescriptor(const std::vector<std::string>& _inputs, double _cross_section);

    bool HasCrossSection() const;
    double GetCrossSectionWeight() const;
};

struct SkimJob {
    std::string name;
    std::string merged_output;
    std::vector<FileDescriptor> files;
    bool apply_common_weights{true};
    bool isData{false};
    mc_corrections::WeightingMode weights;

    bool ProduceMergedOutput() const;
};

using SkimJobCollection = std::unordered_map<std::string, SkimJob>;

} // namespace tuple_skimmer
} // namespace analysis
