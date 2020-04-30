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
    std::set<SignalMode> mode;
    DiscriminatorWP btag_wp;
    bool use_cache;
    mc_corrections::WeightingMode common_weights;
    BTaggerKind jet_ordering;
    unsigned n_splits{0};
    unsigned split_seed{0};
    std::string xs_cfg;

    //light setup
    bool apply_mass_cut{false}, apply_charge_cut{false}, apply_bb_cut{true}, apply_tau_iso{false};
    bool keep_genJets{false}, keep_genParticles{false}, keep_MET_cov{true};
    bool apply_kinfit{true};
    bool applyTauId{true};
    std::set<std::string> tau_id_cuts;
    std::set<size_t> tau_id_cut_indices;

    std::map<SelectionCut,analysis::EllipseParameters> massWindowParams;

    bool ApplyTauIdCuts() const;

};

using SetupCollection = std::unordered_map<std::string, Setup>;

class CrossSectionProvider{
public:
    CrossSectionProvider(std::string cross_section_cfg);
    double GetCrossSection(const std::string& process) const;

private:
    PropertyConfigReader xs_items;
};

struct FileDescriptor {
    std::vector<std::string> inputs;
    std::string output;
    std::string cross_section;
    bool first_input_is_ref;
    std::vector<bool> input_is_partial;

    FileDescriptor(const std::vector<std::string>& _inputs, const std::string& _output, const std::string& _cross_section);

    bool HasCrossSection() const;
    double GetCrossSection(const CrossSectionProvider& xs_provider) const;

};

struct SkimJob {
    std::string name;
    std::string merged_output;
    std::vector<FileDescriptor> files;
    bool apply_common_weights{true};
    bool isData{false};
    mc_corrections::WeightingMode weights;
    std::string cross_section;

    bool ProduceMergedOutput() const;
    double GetCrossSection(const CrossSectionProvider& xs_provider) const;
};

using SkimJobCollection = std::unordered_map<std::string, SkimJob>;


} // namespace tuple_skimmer
} // namespace analysis
