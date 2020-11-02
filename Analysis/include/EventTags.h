/*! Definition of event tags.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "hh-bbtautau/Analysis/include/EventAnalyzerDataId.h"

namespace analysis {
namespace bbtautau {

struct EventTags {
    using DataId = EventAnalyzerDataId;

    std::vector<DataId> dataIds;
    std::vector<double> weights;
};

class EventTagCreator {
public:
    using DataId = EventTags::DataId;

    EventTagCreator(const EventCategorySet& _categories, const EventSubCategorySet& _subCategories,
                    const std::map<SelectionCut, analysis::EllipseParameters>& _massWindowParams,
                    const std::set<UncertaintySource>& _event_unc_sources,
                    const std::set<UncertaintySource>& _norm_unc_sources,
                    bool _use_kinFit, bool _use_svFit);

    int CreateVBFTag(const LorentzVectorM& VBF1_p4, const LorentzVectorM& VBF2_p4, bool is_VBF,
                     bool pass_vbf_trigger) const;

    EventTags CreateEventTags(const std::vector<DataId>& dataIds_base, const std::vector<double>& weights,
                              int num_central_jets, bool has_b_pair, int num_btag_loose, int num_btag_medium,
                              int num_btag_tight, bool is_vbf, bool is_boosted, int vbf_tag_raw,
                              const LorentzVectorM& SVfit_p4, const LorentzVectorM& MET_p4,
                              double m_bb, double m_tt_vis, int kinFit_convergence) const;

private:
    const EventCategorySet& categories;
    const EventSubCategorySet& subCategories;
    const std::map<SelectionCut,analysis::EllipseParameters>& massWindowParams;
    const std::set<UncertaintySource>& event_unc_sources, norm_unc_sources;
    const bool use_kinFit, use_svFit;
};

} // namespace bbtautau
} // namespace analysis
