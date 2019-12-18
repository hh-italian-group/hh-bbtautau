/*! Various hh event weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "hh-bbtautau/McCorrections/include/EventWeights_HH.h"

#include "hh-bbtautau/McCorrections/include/HH_nonResonant_weight.h"
#include "hh-bbtautau/McCorrections/include/TTbar_weight.h"
#include "hh-bbtautau/McCorrections/include/NJets_HT_weight.h"

namespace analysis {
namespace mc_corrections {

EventWeights_HH::EventWeights_HH(Period period, JetOrdering jet_ordering, DiscriminatorWP btag_wp,
                                 bool use_LLR_weights, bool applyTauId, WeightingMode mode) :
    EventWeights(period, jet_ordering, btag_wp, applyTauId, mode)
{
    if (period == Period::Run2016){
        std::string dy_weights = use_LLR_weights ? Full_Cfg_Name("2016/dyjets_weights_LLR.cfg")
                                                 : Full_Cfg_Name("2016/dyjets_weights.cfg");
        if(mode.empty() || mode.count(WeightType::DY))
            providers[WeightType::DY] = std::make_shared<NJets_HT_weight>("DY", dy_weights);
        if(mode.empty() || mode.count(WeightType::TTbar))
            providers[WeightType::TTbar] = std::make_shared<TTbar_weight>(Full_Cfg_Name("2016/ttbar_weights_full.cfg"));
        std::string wjet_weights =
                use_LLR_weights ? Full_Cfg_Name("2016/wjets_weights_LLR.cfg") : Full_Cfg_Name("2016/wjets_weights.cfg");
        if(mode.empty() || mode.count(WeightType::Wjets))
            providers[WeightType::Wjets] = std::make_shared<NJets_HT_weight>("Wjets", wjet_weights);
    }
    else if (period == Period::Run2017){
        std::string dy_weights = Full_Cfg_Name("2017/dyjets_weights_2017.cfg");
        if(mode.empty() || mode.count(WeightType::DY))
            providers[WeightType::DY] = std::make_shared<NJets_HT_weight>("DY", dy_weights);
        std::string wjet_weights = Full_Cfg_Name("2017/wjets_weights_2017.cfg");
        if(mode.empty() || mode.count(WeightType::Wjets))
            providers[WeightType::Wjets] = std::make_shared<NJets_HT_weight>("Wjets", wjet_weights);

    }
    else if (period == Period::Run2018){
        std::string dy_weights = Full_Cfg_Name("2018/dyjets_weights_2018.cfg");
        if(mode.empty() || mode.count(WeightType::DY))
            providers[WeightType::DY] = std::make_shared<NJets_HT_weight>("DY", dy_weights);
        std::string wjet_weights = Full_Cfg_Name("2018/wjets_weights_2018.cfg");
        if(mode.empty() || mode.count(WeightType::Wjets))
            providers[WeightType::Wjets] = std::make_shared<NJets_HT_weight>("Wjets", wjet_weights);

    }
    else
        throw exception("Period %1% is not supported (EventWeights_HH).") % period;

    if(mode.empty() || mode.count(WeightType::BSM_to_SM))
        providers[WeightType::BSM_to_SM] = std::make_shared<NonResHH_EFT::WeightProvider>(
                    FullBSMtoSM_Name("coefficientsByBin_extended_3M_costHHSim_19-4.txt"));
}

ntuple::ProdSummary EventWeights_HH::GetSummaryWithWeights(const std::shared_ptr<TFile>& file,
                                                           const WeightingMode& weighting_mode) const
{
    static const WeightingMode shape_weights = { WeightType::PileUp, WeightType::BSM_to_SM, WeightType::DY,
                                               WeightType::TTbar, WeightType::Wjets, WeightType::GenEventWeight };
    static const WeightingMode shape_weights_withTopPt = shape_weights | WeightingMode({WeightType::TopPt});

    auto summary_tuple = ntuple::CreateSummaryTuple("summary", file.get(), true, ntuple::TreeState::Full);
    auto summary = ntuple::MergeSummaryTuple(*summary_tuple);
    summary.totalShapeWeight = 0;
    summary.totalShapeWeight_withTopPt = 0;

    const auto mode = shape_weights & weighting_mode;
    const auto mode_withTopPt = shape_weights_withTopPt & weighting_mode;
    const bool calc_withTopPt = mode_withTopPt.count(WeightType::TopPt);

    if(mode.size() || mode_withTopPt.size()) {
        auto all_events = ntuple::CreateExpressTuple("all_events", file.get(), true,
                                                       ntuple::TreeState::Full);

        using EventIdSet = std::set<EventIdentifier>;
        EventIdSet processed_events;
        for(const auto& event : *all_events) {
            const EventIdentifier Id(event.run, event.lumi, event.evt);
            if(processed_events.count(Id)) {
//                    std::cout << "WARNING: duplicated express event " << Id << std::endl;
                continue;
            }
            processed_events.insert(Id);
            summary.totalShapeWeight += GetTotalWeight(event, mode);
            if(calc_withTopPt)
                summary.totalShapeWeight_withTopPt += GetTotalWeight(event, mode_withTopPt);
        }
    }
    else{
        try{
            auto all_events = ntuple::CreateExpressTuple("all_events", file.get(), true,
                                                         ntuple::TreeState::Full);
            summary.totalShapeWeight = all_events->GetEntries();
            if(calc_withTopPt)
                summary.totalShapeWeight_withTopPt = all_events->GetEntries();
        } catch(std::exception& ) {}
    }
    return summary;
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
