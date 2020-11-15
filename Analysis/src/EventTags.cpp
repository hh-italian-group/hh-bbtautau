/*! Definition of event tags.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/EventTags.h"
#include "h-tautau/Cuts/include/hh_bbtautau_Run2.h"

namespace analysis {
namespace bbtautau {

EventTagCreator::EventTagCreator(const EventCategorySet& _categories, const EventSubCategorySet& _subCategories,
                                  const std::set<UncertaintySource>& _event_unc_sources,
                                  const std::set<UncertaintySource>& _norm_unc_sources,
                                  bool _use_IterativeFit, Channel _channel, Period _period) :
        categories(_categories), subCategories(_subCategories), event_unc_sources(_event_unc_sources),
        norm_unc_sources(_norm_unc_sources), use_IterativeFit(_use_IterativeFit), channel(_channel), period(_period)
{
    if(use_IterativeFit) {
        static const std::map<std::pair<Channel, Period>, float> iterativeFit_corrections = {
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
        auto iter = iterativeFit_corrections.find({channel, period});
        if(iter == iterativeFit_corrections.end())
            throw exception("Normalization correction for iterative fit is not available for %1% %2%")
                    % period % channel;
        iterativeFit_correction = iter->second;
    }
}

std::pair<float, VBF_Category> EventTagCreator::FindVBFCategory(float dnn_score_TT_dl, float dnn_score_TT_sl,
        float dnn_score_TT_lep, float dnn_score_TT_FH, float dnn_score_DY, float dnn_score_ggHH, float dnn_score_ttH,
        float dnn_score_ttH_tautau, float dnn_score_tth_bb, float dnn_score_qqHH, float dnn_score_qqHH_vbf_c2v,
        float dnn_score_qqHH_sm)
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

EventTags EventTagCreator::CreateEventTags(const DataId& dataId_base, float weight, bool is_data,
        float weight_btag_Loose, float weight_btag_Medium, float weight_btag_Tight, float weight_btag_IterativeFit,
        int num_central_jets, bool has_b_pair, int num_btag_loose, int num_btag_medium, int num_btag_tight,
        bool is_vbf, bool is_boosted, const std::pair<float,VBF_Category>& vbf_cat, const LorentzVectorM& SVfit_p4,
        float m_bb, float /*m_tt_vis*/, int kinFit_convergence, int SVfit_valid) const
{
    const std::map<DiscriminatorWP, size_t> bjet_counts = {
        { DiscriminatorWP::Loose, num_btag_loose },
        { DiscriminatorWP::Medium, num_btag_medium },
        { DiscriminatorWP::Tight, num_btag_tight}
    };

    const std::map<DiscriminatorWP, float > btag_weights = {
        { DiscriminatorWP::Loose, weight_btag_Loose },
        { DiscriminatorWP::Medium, weight_btag_Medium },
        { DiscriminatorWP::Tight, weight_btag_Tight },
    };
    const auto get_weight_btag = [&](DiscriminatorWP wp) {
        if(use_IterativeFit) return weight_btag_IterativeFit * iterativeFit_correction;
        return btag_weights.at(wp);
    };

    EventTags evt_tags;

    if(!event_unc_sources.count(dataId_base.Get<UncertaintySource>()))
        return evt_tags;

    for(const auto& category : categories) {
        if(!category.Contains(static_cast<size_t>(num_central_jets), bjet_counts, is_vbf, is_boosted, vbf_cat.second))
            continue;
        EventSubCategory evtSubCategory;
        const float btag_sf = category.HasBtagConstraint() && !is_data ? get_weight_btag(category.BtagWP()) : 1.f;
        const float cat_weight = weight * btag_sf;


        if(has_b_pair) {
            const EllipseParameters& window = category.HasBoostConstraint() && category.IsBoosted()
                    ? cuts::hh_bbtautau_Run2::hh_tag::boosted_window : cuts::hh_bbtautau_Run2::hh_tag::resolved_window;
            evtSubCategory.SetCutResult(SelectionCut::mh, SVfit_valid && window.IsInside(SVfit_p4.M(), m_bb));
            evtSubCategory.SetCutResult(SelectionCut::KinematicFitConverged, kinFit_convergence > 0);
        }

        for(const auto& subCategory : subCategories) {
            if(!evtSubCategory.Implies(subCategory)) continue;
            const auto& dataId = dataId_base.Set(category).Set(subCategory);
            evt_tags.dataIds.push_back(dataId);
            evt_tags.weights.push_back(cat_weight);
        }
    }

    return evt_tags;
}
} // namespace bbtautau
} // namespace analysis
