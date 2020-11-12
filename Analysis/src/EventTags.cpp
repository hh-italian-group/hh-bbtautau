/*! Definition of event tags.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/EventTags.h"
#include "h-tautau/Cuts/include/hh_bbtautau_Run2.h"

namespace analysis {
namespace bbtautau {

EventTagCreator::EventTagCreator(const EventCategorySet& _categories, const EventSubCategorySet& _subCategories,
                                  const std::map<SelectionCut, analysis::EllipseParameters>& _massWindowParams,
                                  const std::set<UncertaintySource>& _event_unc_sources,
                                  const std::set<UncertaintySource>& _norm_unc_sources,
                                  bool _use_IterativeFit, const Channel _channel, const Period _period ):
        categories(_categories), subCategories(_subCategories), massWindowParams(_massWindowParams),
        event_unc_sources(_event_unc_sources), norm_unc_sources(_norm_unc_sources), use_IterativeFit(_use_IterativeFit),
        channel(_channel),period(_period)
{
}

std::pair<float, VBF_Category> EventTagCreator::FindVBFCategory(float dnn_score_TT_dl, float dnn_score_TT_sl, float dnn_score_TT_lep,
                                                                 float dnn_score_TT_FH, float dnn_score_DY, float dnn_score_ggHH,
                                                                 float dnn_score_ttH, float dnn_score_ttH_tautau, float dnn_score_tth_bb,
                                                                 float dnn_score_qqHH, float dnn_score_qqHH_vbf_c2v, float dnn_score_qqHH_sm)
{
    std::vector<std::pair<float, VBF_Category>> category_and_dnn = {
        { dnn_score_qqHH_sm + dnn_score_qqHH + dnn_score_qqHH_vbf_c2v, VBF_Category::qqHH },
        { dnn_score_ggHH , VBF_Category::ggHH},
        { dnn_score_TT_dl + dnn_score_TT_lep + dnn_score_TT_sl , VBF_Category::TT_L },
        { dnn_score_TT_FH, VBF_Category::TT_FH},
        { dnn_score_ttH + dnn_score_ttH_tautau + dnn_score_tth_bb, VBF_Category::ttH},
        { dnn_score_DY, VBF_Category::DY},
        };
    return *std::max_element(category_and_dnn.begin(), category_and_dnn.end());
}

EventTags EventTagCreator::CreateEventTags(const std::vector<DataId>& dataIds_base,
                                            const std::vector<float>& weights,
                                            const std::vector<float>& btag_weight_Loose,
                                            const std::vector<float>& btag_weight_Medium,
                                            const std::vector<float>& btag_weight_Tight,
                                            const std::vector<float>& btag_weight_IterativeFit,
                                            int num_central_jets, bool has_b_pair,
                                            int num_btag_loose, int num_btag_medium, int num_btag_tight,
                                            bool is_vbf, bool is_boosted, const std::pair<float,VBF_Category>& vbf_cat,
                                            const LorentzVectorM& SVfit_p4, const LorentzVectorM& MET_p4,
                                            float m_bb, float m_tt_vis, int kinFit_convergence, int SVfit_valid) const
{
    static const EventCategory evtCategory_2j = EventCategory::Parse("2j");
    const std::map<std::pair<Channel, Period>, float> iterativeFit_correction = {
    {{Channel::ETau, Period::Run2016}, 1.0120}, // 1.01203715762
    {{Channel::MuTau, Period::Run2016}, 1.013}, // 1.01296331792
    {{Channel::TauTau, Period::Run2016}, 1.0101}, // 1.01010364517
    {{Channel::ETau, Period::Run2017}, 0.9949}, // 0.994930642536
    {{Channel::MuTau, Period::Run2017}, 0.9993}, // 0.999263715405
    {{Channel::TauTau, Period::Run2017}, 0.9547}, // 0.954725518543
    {{Channel::ETau, Period::Run2018}, 1.0004}, // 1.0003768284
    {{Channel::MuTau, Period::Run2018}, 1.0039}, // 1.00388888346
    {{Channel::TauTau, Period::Run2018}, 0.9795}, // 0.97949255088
    };
    const std::map<DiscriminatorWP, size_t> bjet_counts = {
        { DiscriminatorWP::Loose, num_btag_loose },
        { DiscriminatorWP::Medium, num_btag_medium },
        { DiscriminatorWP::Tight, num_btag_tight}
    };

    float iterativeFit_corr = iterativeFit_correction.at({channel, period});

    std::map<DiscriminatorWP, float > btag_weights = {};
    if((btag_weight_Loose.size()!=0 || btag_weight_Medium.size()!=0 || btag_weight_Tight.size()!=0) && !use_IterativeFit){
        btag_weights = {
            { DiscriminatorWP::Loose, btag_weight_Loose.at(0) },
            { DiscriminatorWP::Medium, btag_weight_Medium.at(0) },
            { DiscriminatorWP::Tight, btag_weight_Tight.at(0)},
        };
    }
    else if(btag_weight_IterativeFit.size()!=0 && use_IterativeFit){
        btag_weights = {
            { DiscriminatorWP::Loose, btag_weight_IterativeFit.at(0)*iterativeFit_corr },
            { DiscriminatorWP::Medium, btag_weight_IterativeFit.at(0)*iterativeFit_corr },
            { DiscriminatorWP::Tight, btag_weight_IterativeFit.at(0)*iterativeFit_corr },
        };
    }

    EventTags evt_tags;

    for (size_t i = 0; i < dataIds_base.size(); ++i) {
        const auto& dataId_base = dataIds_base.at(i);
        float weight = weights.at(i);

        if(dataId_base.Get<EventCategory>() != evtCategory_2j
                || !event_unc_sources.count(dataId_base.Get<UncertaintySource>()))
            continue;

        for(const auto& category : categories) {
            if(!category.Contains(static_cast<size_t>(num_central_jets), bjet_counts, is_vbf, is_boosted, vbf_cat.second))
                continue;
            EventSubCategory evtSubCategory;
            float btag_sf = 1;
            if(category.HasBtagConstraint() && !btag_weights.empty()){
                  btag_sf = btag_weights.at(category.BtagWP());
            }
            weight *= btag_sf;


            if(has_b_pair) {
              const EllipseParameters& window = category.HasBoostConstraint() && category.IsBoosted() ? cuts::hh_bbtautau_Run2::hh_tag::boosted_window : cuts::hh_bbtautau_Run2::hh_tag::resolved_window;
              evtSubCategory.SetCutResult(SelectionCut::mh, SVfit_valid && window.IsInside(SVfit_p4.M(), m_bb));
              evtSubCategory.SetCutResult(SelectionCut::KinematicFitConverged, kinFit_convergence>0);
            }

            for(const auto& subCategory : subCategories) {
                if(!evtSubCategory.Implies(subCategory)) continue;
                const auto& dataId = dataId_base.Set(category).Set(subCategory);
                evt_tags.dataIds.push_back(dataId);
                evt_tags.weights.push_back(weight);
            }
        }
    }

    return evt_tags;
}
} // namespace bbtautau
} // namespace analysis
