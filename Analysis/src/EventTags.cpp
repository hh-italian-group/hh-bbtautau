/*! Definition of event tags.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/EventTags.h"
#include "h-tautau/Cuts/include/hh_bbtautau_Run2.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
/*

1. in the constructor you load json into std::map<std::pair<UncSource, UncScale>, float> ---> Fatto
2. when you calculate btag weight you know what is uncSource and uncScale -> multiply it by the value from the map
for btag related unc, you should multiply by both sf and unc up/down correction
*/

namespace analysis {
namespace bbtautau {
/*namespace{
std::string get_name(UncertaintySource unc_source, UncertaintyScale unc_scale){
    std::string complete_name;
    if(unc_scale!= UncertaintyScale::Central && unc_scale!= UncertaintyScale::Up && unc_scale!= UncertaintyScale::Down)
        throw exception("Incorrect uncertainty scale");

    if(unc_source==UncertaintySource::None){
        if(unc_scale==UncertaintyScale::Central)
            complete_name="Central";
        else
            throw exception("Central Value can't have different scales");
    }
    else{
        if(unc_scale!=UncertaintyScale::Central)
            complete_name=ToString(unc_source)+"_"+ToString(unc_scale);
        else{
            complete_name=ToString(unc_source);
        }
    }
    return complete_name;

}*/
}

EventTagCreator::EventTagCreator(const EventCategorySet& _categories, const EventSubCategorySet& _subCategories,
                                  const std::set<UncertaintySource>& _event_unc_sources,
                                  const std::set<UncertaintySource>& _norm_unc_sources,
                                  bool _use_IterativeFit, std::string _json_file)  :
        categories(_categories), subCategories(_subCategories), event_unc_sources(_event_unc_sources),
        norm_unc_sources(_norm_unc_sources), use_IterativeFit(_use_IterativeFit),
        json_file(_json_file)
{
    if(use_IterativeFit) {
        std::vector<UncertaintyScale> unc_scales= {UncertaintyScale::Central, UncertaintyScale::Up, UncertaintyScale::Down};
        boost::property_tree::ptree loadPtreeRoot;
        boost::property_tree::read_json(json_file, loadPtreeRoot);
        for(UncertaintySource unc_source:event_unc_sources){
            boost::property_tree::ptree::const_assoc_iterator it = loadPtreeRoot.find(ToString(unc_source));
            if( it == loadPtreeRoot.not_found() ) continue;
            auto unc_source_s = loadPtreeRoot.get_child(ToString(unc_source));
                for (UncertaintyScale unc_scale: unc_scales ){
                    boost::property_tree::ptree::const_assoc_iterator it2 = unc_source_s.find(ToString(unc_scale));
                    if( it2 == unc_source_s.not_found() ) continue;
                        std::string btag_corr_to_find= ToString(unc_source)+"."+ToString(unc_scale);
                        float btag_corr_value = loadPtreeRoot.get<float>((btag_corr_to_find).c_str());
                        //std::cout << typeid(unc_source).name() << std::endl;
                        //std::cout << typeid(unc_scale).name() << std::endl;
                        //std::cout << typeid(btag_corr_value).name() << std::endl;
                        std::pair<UncertaintySource, UncertaintyScale> pair = {unc_source,unc_scale};
                        iterativeFit_corrections.insert({pair,btag_corr_value});
                    }
        }

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
        { dnn_score_TT_dl + dnn_score_TT_lep + dnn_score_TT_sl+dnn_score_TT_FH, VBF_Category::TT},
        { dnn_score_ttH + dnn_score_ttH_tautau + dnn_score_tth_bb, VBF_Category::ttH},
        { dnn_score_DY, VBF_Category::DY},
    };
    return *std::max_element(category_and_dnn.begin(), category_and_dnn.end());
}

static std::set<UncertaintySource> btag_uncs = {
    UncertaintySource::btag_lf, UncertaintySource::btag_hf, UncertaintySource::btag_hfstats1,
    UncertaintySource::btag_hfstats2, UncertaintySource::btag_lfstats1,
    UncertaintySource::btag_lfstats2, UncertaintySource::btag_cferr1, UncertaintySource::btag_cferr2
};


EventTags EventTagCreator::CreateEventTags(const DataId& dataId_base, float weight, bool is_data,
        float weight_btag_Loose, float weight_btag_Medium, float weight_btag_Tight, float weight_btag_IterativeFit,
        int num_central_jets, bool has_b_pair, int num_btag_loose, int num_btag_medium, int num_btag_tight,
        bool is_vbf, bool is_boosted, const std::pair<float,VBF_Category>& vbf_cat, const LorentzVectorM& SVfit_p4,
        const UncMap& unc_map, float m_bb, float /*m_tt_vis*/, int kinFit_convergence, int SVfit_valid) const
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
    const auto get_weight_btag = [&](DiscriminatorWP wp, UncertaintySource unc_source, UncertaintyScale unc_scale) {
        if(use_IterativeFit) {
            if(iterativeFit_corrections.find((unc_source, unc_scale))){
                return weight_btag_IterativeFit * iterativeFit_corrections.at((unc_source, unc_scale));
            }
            else
                return weight_btag_IterativeFit ; //* iterativeFit_corrections.at((UncertaintySource::None, UncertaintyScale::Central))
        }
        return btag_weights.at(wp);
    };

    EventTags evt_tags;
    if(!event_unc_sources.count(dataId_base.Get<UncertaintySource>()))
        return evt_tags;


    for(const auto& category : categories) {
        if(!category.Contains(static_cast<size_t>(num_central_jets), bjet_counts, is_vbf, is_boosted, vbf_cat.second))
            continue;
        EventSubCategory evtSubCategory;
        const float btag_sf = category.HasBtagConstraint() && !is_data ? get_weight_btag(category.BtagWP(), dataId_base.Get<UncertaintySource>(), dataId_base.Get<UncertaintyScale>()) : 1.f;
        //const float cat_weight = weight * btag_sf;
        const float cat_weight = btag_uncs.count(dataId_base.Get<UncertaintySource>())? weight * btag_sf * unc_map.at(category.BtagWP()) : weight * btag_sf ;

        if(has_b_pair) {
            const EllipseParameters& window = category.HasBoostConstraint() && category.IsBoosted()
                    ? cuts::hh_bbtautau_Run2::hh_tag::boosted_window : cuts::hh_bbtautau_Run2::hh_tag::resolved_window;
            evtSubCategory.SetCutResult(SelectionCut::mh, SVfit_valid && window.IsInside(SVfit_p4.M(), m_bb));
            evtSubCategory.SetCutResult(SelectionCut::KinematicFitConverged, kinFit_convergence > 0);
        }
        for(const auto& subCategory : subCategories) {
            if(!evtSubCategory.Implies(subCategory)) continue;
            const auto dataId = dataId_base.Set(category).Set(subCategory);
            evt_tags.dataIds.push_back(dataId);
            evt_tags.weights.push_back(cat_weight);
            if(!is_data && dataId_base.Get<UncertaintySource>()==UncertaintySource::None) {
                for(const auto& [key, weight_shift] : unc_map) {
                    if(!category.HasBtagConstraint() && btag_uncs.count(key.first)) continue;
                    const auto dataId_scaled = dataId.Set(key.first).Set(key.second);
                    const float btag_sf_scaled = get_weight_btag(category.BtagWP(), dataId_scaled.Get<UncertaintySource>(), dataId_scaled.Get<UncertaintyScale>());
                    const float cat_weight_shifted = btag_uncs.count(dataId_scaled.Get<UncertaintySource>())? weight_shift * btag_sf * unc_map.at(category.BtagWP()) : weight_shift * btag_sf_scaled ;
                    evt_tags.dataIds.push_back(dataId_scaled);
                    evt_tags.weights.push_back(cat_weight *weight_shift);
                }
            }
        }


    }

    return evt_tags;
}
} // namespace bbtautau
} // namespace analysis
