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
                    const std::set<UncertaintySource>& _norm_unc_sources);

    static const std::pair<float, VBF_Category>& FindVBFCategory(float dnn_score_TT_dl, float dnn_score_DY, float dnn_score_TT_lep,
                                                     float dnn_score_qqHH_sm, float dnn_score_ggHH, float dnn_score_TT_FH,
                                                     float dnn_score_ttH, float dnn_score_ttH_tautau, float dnn_score_tth_bb,
                                                     float dnn_score_qqHH, float dnn_score_qqHH_vbf_c2v, float dnn_score_TT_sl, std::string mdnn_version) const;
    //int CreateVBFTag(const LorentzVectorM& VBF1_p4, const LorentzVectorM& VBF2_p4, bool is_VBF,
    //                 bool pass_vbf_trigger) const;

    EventTags CreateEventTags(const std::vector<DataId>& dataIds_base, const std::vector<double>& weights,
                              const std::vector<float>& btag_weight_Loose,
                              const std::vector<float>& btag_weight_Medium,
                              const std::vector<float>& btag_weight_Tight,
                              int num_central_jets, bool has_b_pair, int num_btag_loose, int num_btag_medium,
                              int num_btag_tight, bool is_vbf, bool is_boosted, std::pair<float, VBF_Category> vbf_cat,
                              const LorentzVectorM& SVfit_p4, const LorentzVectorM& MET_p4,
                              double m_bb, double m_tt_vis, int kinFit_convergence) const;

private:
    const EventCategorySet& categories;
    const EventSubCategorySet& subCategories;
    const std::map<SelectionCut,analysis::EllipseParameters>& massWindowParams;
    const std::set<UncertaintySource>& event_unc_sources, norm_unc_sources;
    //const bool use_kinFit, use_svFit;
};

} // namespace bbtautau
} // namespace analysis
