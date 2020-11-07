/*! Various hh event weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "hh-bbtautau/McCorrections/include/EventWeights_HH.h"

#include "hh-bbtautau/McCorrections/include/HH_nonResonant_weight.h"
#include "hh-bbtautau/McCorrections/include/TTbar_weight.h"
#include "hh-bbtautau/McCorrections/include/NJets_HT_weight.h"

namespace analysis {
namespace mc_corrections {

EventWeights_HH::EventWeights_HH(Period period, const BTagger& bTagger, const WeightingMode& mode) :
        EventWeights(period, bTagger, mode)
{
    if (period == Period::Run2016) {
        if(mode.empty() || mode.count(WeightType::DY))
            CreateProvider<NJets_HT_weight>(WeightType::DY, "DY", Full_Cfg_Name("2016/dyjets_weights_2016.cfg"));
        if(mode.empty() || mode.count(WeightType::TTbar))
            CreateProvider<TTbar_weight>(WeightType::TTbar, Full_Cfg_Name("2016/ttbar_weights_full.cfg"));
        if(mode.empty() || mode.count(WeightType::Wjets))
            CreateProvider<NJets_HT_weight>(WeightType::Wjets, "Wjets", Full_Cfg_Name("2016/wjets_weights_2016.cfg"));
    }
    else if (period == Period::Run2017){
        if(mode.empty() || mode.count(WeightType::DY))
            CreateProvider<NJets_HT_weight>(WeightType::DY, "DY", Full_Cfg_Name("2017/dyjets_weights_2017.cfg"));
        if(mode.empty() || mode.count(WeightType::Wjets))
            CreateProvider<NJets_HT_weight>(WeightType::Wjets, "Wjets", Full_Cfg_Name("2017/wjets_weights_2017.cfg"));
    }
    else if (period == Period::Run2018){
        if(mode.empty() || mode.count(WeightType::DY))
            CreateProvider<NJets_HT_weight>(WeightType::DY, "DY", Full_Cfg_Name("2018/dyjets_weights_2018.cfg"));
        if(mode.empty() || mode.count(WeightType::Wjets))
            CreateProvider<NJets_HT_weight>(WeightType::Wjets, "Wjets", Full_Cfg_Name("2018/wjets_weights_2018.cfg"));
    }
    else
        throw exception("Period %1% is not supported (EventWeights_HH).") % period;

    if(mode.empty() || mode.count(WeightType::BSM_to_SM))
        CreateProvider<NonResHH_EFT::WeightProvider>(WeightType::BSM_to_SM,
                FullBSMtoSM_Name("coefficientsByBin_extended_3M_costHHSim_19-4.txt"));
}

ntuple::ProdSummary EventWeights_HH::GetSummaryWithWeights(const std::shared_ptr<TFile>& file,
                                                           const WeightingMode& weighting_mode,
                                                           const boost::optional<double>& max_gen_weight,
                                                           bool control_duplicates) const
{
    static const WeightingMode shape_weights(WeightType::PileUp, WeightType::BSM_to_SM, WeightType::DY,
                                             WeightType::TTbar, WeightType::Wjets, WeightType::GenEventWeight);
    static const WeightingMode shape_weights_withTopPt = shape_weights | WeightingMode(WeightType::TopPt);

    auto summary_tuple = ntuple::CreateSummaryTuple("summary", file.get(), true, ntuple::TreeState::Full);
    auto summary = ntuple::MergeSummaryTuple(*summary_tuple);
    summary.totalShapeWeight = 0;
    summary.totalShapeWeight_withTopPt = 0;
    summary.totalShapeWeight_withPileUp_Up = 0;
    summary.totalShapeWeight_withPileUp_Down = 0;


    const auto mode = shape_weights & weighting_mode;
    //TopPt
    const auto mode_withTopPt = shape_weights_withTopPt & weighting_mode;
    const bool calc_withTopPt = mode_withTopPt.count(WeightType::TopPt);
    //PileUp
    auto mode_withoutPileUp = mode;
    mode_withoutPileUp.erase(WeightType::PileUp);
    const bool calc_withPileUp = mode.count(WeightType::PileUp);

    if(mode.size() || mode_withTopPt.size()) {
        auto all_events = ntuple::CreateExpressTuple("all_events", file.get(), true, ntuple::TreeState::Full);

        using EventIdSet = std::set<EventIdentifier>;
        EventIdSet processed_events;
        auto pu_weight_provider = GetProviderT<mc_corrections::PileUpWeightEx>(mc_corrections::WeightType::PileUp);

        for(const auto& event : *all_events) {
            if(max_gen_weight && std::abs(event.genEventWeight) > *max_gen_weight) continue;
            if(control_duplicates) {
                const EventIdentifier Id(event.run, event.lumi, event.evt);
                if(processed_events.count(Id)) {
                    // std::cout << "WARNING: duplicated express event " << Id << std::endl;
                    continue;
                }
                processed_events.insert(Id);
            }

            summary.totalShapeWeight += GetTotalWeight(event, mode);
            if(calc_withTopPt)
                summary.totalShapeWeight_withTopPt += GetTotalWeight(event, mode_withTopPt);
            if(calc_withPileUp){
                summary.totalShapeWeight_withPileUp_Up += GetTotalWeight(event, mode_withoutPileUp)
                                                          * pu_weight_provider->Get(event, UncertaintyScale::Up);
                summary.totalShapeWeight_withPileUp_Down += GetTotalWeight(event, mode_withoutPileUp)
                                                          * pu_weight_provider->Get(event, UncertaintyScale::Down);
            }
        }
    }
    else{
        try{
            auto all_events = ntuple::CreateExpressTuple("all_events", file.get(), true,
                                                         ntuple::TreeState::Full);
            summary.totalShapeWeight = all_events->GetEntries();
            if(calc_withTopPt)
                summary.totalShapeWeight_withTopPt = all_events->GetEntries();
            if(calc_withPileUp){
                summary.totalShapeWeight_withPileUp_Up = all_events->GetEntries();
                summary.totalShapeWeight_withPileUp_Down = all_events->GetEntries();
            }
        } catch(std::exception& ) {}
    }
    return summary;
}

std::map<UncertaintyScale, std::vector<double>> EventWeights_HH::GetTotalShapeWeights(
        const std::shared_ptr<TFile>& file, const WeightingMode& weighting_mode,
        const std::vector<NonResHH_EFT::Point>& eft_points, bool orthogonal) const
{
    static const WeightingMode shape_weights(WeightType::PileUp, WeightType::BSM_to_SM, WeightType::DY,
                                             WeightType::TTbar, WeightType::Wjets, WeightType::GenEventWeight);
    size_t N = eft_points.size();
    std::map<UncertaintyScale, std::vector<double>> total_weights_scale;
    for(const auto& scale : GetAllUncertaintyScales())
        total_weights_scale[scale].resize(N, 0);

    const auto mode = shape_weights & weighting_mode;
    auto mode_withoutPileUp = mode;
    mode_withoutPileUp.erase(WeightType::PileUp);

    if(mode.size()) {
        auto eft_weights_provider = GetProviderT<NonResHH_EFT::WeightProvider>(WeightType::BSM_to_SM);
        auto all_events = ntuple::CreateExpressTuple("all_events", file.get(), true, ntuple::TreeState::Full);
        auto pu_weight_provider = GetProviderT<mc_corrections::PileUpWeightEx>(mc_corrections::WeightType::PileUp);

        for(const auto& event : *all_events) {
            for(const auto& scale : GetAllUncertaintyScales()){
                for(size_t n = 0; n < N; ++n) {
                    if(orthogonal && (event.evt % N) != n) continue;
                    eft_weights_provider->SetTargetPoint(eft_points.at(n));
                    const double pu_weight = weighting_mode.count(WeightType::PileUp) ?
                                             pu_weight_provider->Get(event, scale) : 1.;
                    total_weights_scale.at(scale).at(n) += GetTotalWeight(event, mode_withoutPileUp) * pu_weight;
                }
            }
        }
    }
    return total_weights_scale;
}

std::string EventWeights_HH::FullBSMtoSM_Name(const std::string& fileName)
{
    static const std::string path = "hh-bbtautau/McCorrections/data";
    return FullName(fileName, path);
}

std::string EventWeights_HH::Full_Cfg_Name(const std::string& fileName)
{
    static const std::string path = "hh-bbtautau/Analysis/config";
    return FullName(fileName, path);
}

} // namespace mc_corrections
} // namespace analysis
