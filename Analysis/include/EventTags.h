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
    std::vector<float> weights;
};

struct btag_corrections{
    UncertaintySource unc_source;
    UncertaintyScale unc_scale;
    bool tune;
    constexpr bool operator<(const btag_corrections & other) const
    {
        if(unc_source != other.unc_source) return unc_source < other.unc_source;
        if(unc_scale != other.unc_scale) return unc_scale < other.unc_scale;
        return tune < other.tune;
    }
};

class EventTagCreator {
public:
    using DataId = EventTags::DataId;
    using UncMap = std::map<std::pair<UncertaintySource, UncertaintyScale>, float>;

    EventTagCreator(const EventCategorySet& _categories, const EventSubCategorySet& _subCategories,
                    const std::set<UncertaintySource>& _event_unc_sources,
                    const std::set<UncertaintySource>& _norm_unc_sources, bool _use_IterativeFit,
                    std::string _json_file);

    static std::pair<float, VBF_Category> FindVBFCategory(float dnn_score_TT_dl, float dnn_score_TT_sl, float dnn_score_TT_lep,
                                                            float dnn_score_TT_FH, float dnn_score_DY, float dnn_score_ggHH,
                                                            float dnn_score_ttH, float dnn_score_ttH_tautau, float dnn_score_tth_bb,
                                                            float dnn_score_qqHH, float dnn_score_qqHH_vbf_c2v, float dnn_score_qqHH_sm);

    EventTags CreateEventTags(const DataId& dataId_base, float weight, bool is_data, float weight_btag_Loose,
            float weight_btag_Medium, float weight_btag_Tight, float weight_btag_IterativeFit, int num_central_jets,
            bool has_b_pair, int num_btag_loose, int num_btag_medium, int num_btag_tight, bool is_vbf, bool is_boosted, bool tune,
            const std::pair<float,VBF_Category>& vbf_cat, const LorentzVectorM& SVfit_p4, const UncMap& unc_map, float m_bb, float m_tt_vis,
            int kinFit_convergence, int SVfit_valid) const;

private:
    const EventCategorySet& categories;
    const EventSubCategorySet& subCategories;
    const std::set<UncertaintySource>& event_unc_sources, norm_unc_sources;
    const bool use_IterativeFit;
    std::string json_file;
    //float iterativeFit_correction;
    std::map<btag_corrections, float> iterativeFit_corrections;
};

} // namespace bbtautau
} // namespace analysis
