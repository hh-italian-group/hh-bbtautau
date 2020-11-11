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
                                  bool _use_IterativeFit):
        categories(_categories), subCategories(_subCategories), massWindowParams(_massWindowParams),
        event_unc_sources(_event_unc_sources), norm_unc_sources(_norm_unc_sources), use_IterativeFit(_use_IterativeFit)
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

    const std::map<DiscriminatorWP, size_t> bjet_counts = {
        { DiscriminatorWP::Loose, num_btag_loose },
        { DiscriminatorWP::Medium, num_btag_medium },
        { DiscriminatorWP::Tight, num_btag_tight}
    };

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
            { DiscriminatorWP::Loose, btag_weight_IterativeFit.at(0) },
            { DiscriminatorWP::Medium, btag_weight_IterativeFit.at(0) },
            { DiscriminatorWP::Tight, btag_weight_IterativeFit.at(0)},
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
            /*
            if(has_b_pair) {
                if(category.HasBoostConstraint() && category.IsBoosted()) {
                    const bool isInsideBoostedCut = cuts::hh_bbtautau_Run2::hh_tag::IsInsideBoostedMassWindow(SVfit_p4.M(), m_bb);
                    evtSubCategory.SetCutResult(SelectionCut::mh, isInsideBoostedCut);
                }
                else {
                    if( massWindowParams.count(SelectionCut::mh))
                        throw exception("Category mh inconsistent with the false requirement of SVfit.");
                    if(massWindowParams.count(SelectionCut::mh)) {
                        const bool cut_result = massWindowParams.at(SelectionCut::mh).IsInside(SVfit_p4.M(), m_bb);
                        //&& event.GetSVFitResults(ana_setup.allow_calc_svFit).has_valid_momentum
                        evtSubCategory.SetCutResult(SelectionCut::mh, );
                    }
                    if(massWindowParams.count(SelectionCut::mhVis))
                        evtSubCategory.SetCutResult(SelectionCut::mhVis,massWindowParams.at(SelectionCut::mhVis)
                                .IsInside(m_tt_vis, m_bb));

                    if(massWindowParams.count(SelectionCut::mhMET))
                        evtSubCategory.SetCutResult(SelectionCut::mhMET,massWindowParams.at(SelectionCut::mhMET)
                                .IsInside((SVfit_p4+MET_p4).M(), m_bb));

                }
                */
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
