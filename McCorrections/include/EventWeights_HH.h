/*! Various hh event weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/McCorrections/include/EventWeights.h"
#include "HH_nonResonant_weight.h"
#include "TTbar_weight.h"
#include "NJets_HT_weight.h"

namespace analysis {
namespace mc_corrections {

class EventWeights_HH : public EventWeights {
public:
    EventWeights_HH(Period period, DiscriminatorWP btag_wp, bool use_LLR_weights, WeightingMode mode = {}) :
        EventWeights(period, btag_wp, mode)
    {
        if(period != Period::Run2016)
            throw exception("Period %1% is not supported.") % period;
        std::string dy_weights =
                use_LLR_weights ? Full_Cfg_Name("dyjets_weights_LLR.cfg") : Full_Cfg_Name("dyjets_weights.cfg");
        if(mode.empty() || mode.count(WeightType::DY))
            providers[WeightType::DY] = std::make_shared<NJets_HT_weight>("DY", dy_weights);
        if(mode.empty() || mode.count(WeightType::TTbar))
            providers[WeightType::TTbar] = std::make_shared<TTbar_weight>(Full_Cfg_Name("ttbar_weights_full.cfg"));
        if(mode.empty() || mode.count(WeightType::BSM_to_SM))
            providers[WeightType::BSM_to_SM] = std::make_shared<NonResHH_EFT::WeightProvider>(
                        FullBSMtoSM_Name("coefficientsByBin_extended_3M_costHHSim_19-4.txt"));
        std::string wjet_weights =
                use_LLR_weights ? Full_Cfg_Name("wjets_weights_LLR.cfg") : Full_Cfg_Name("wjets_weights.cfg");
        if(mode.empty() || mode.count(WeightType::Wjets))
            providers[WeightType::Wjets] = std::make_shared<NJets_HT_weight>("Wjets", wjet_weights);
    }

    ntuple::ProdSummary GetSummaryWithWeights(const std::shared_ptr<TFile>& file,
                                              const WeightingMode& weighting_mode) const
    {
        static const WeightingMode shape_weights = { WeightType::PileUp, WeightType::BSM_to_SM, WeightType::DY,
                                                   WeightType::TTbar, WeightType::Wjets};
        static const WeightingMode shape_weights_withTopPt = shape_weights | WeightingMode({WeightType::TopPt});

        auto summary_tuple = ntuple::CreateSummaryTuple("summary", file.get(), true, ntuple::TreeState::Full);
        auto summary = ntuple::MergeSummaryTuple(*summary_tuple);
        summary.totalShapeWeight = 0;
        summary.totalShapeWeight_withTopPt = 0;

        const auto mode = shape_weights & weighting_mode;
        const auto mode_withTopPt = shape_weights_withTopPt & weighting_mode;
        const bool calc_withTopPt = mode_withTopPt.count(WeightType::TopPt);
        if(mode.size() || mode_withTopPt.size()) {
            ntuple::ExpressTuple all_events("all_events", file.get(), true);
            for(const auto& event : all_events) {
                summary.totalShapeWeight += GetTotalWeight(event, mode);
                if(calc_withTopPt)
                    summary.totalShapeWeight_withTopPt += GetTotalWeight(event, mode_withTopPt);
            }
        }
        return summary;
    }

private:
    static std::string FullBSMtoSM_Name(const std::string& fileName)
    {
        static const std::string path = "hh-bbtautau/McCorrections/data";
        return FullName(fileName, path);
    }

    static std::string Full_Cfg_Name(const std::string& fileName)
    {
        static const std::string path = "hh-bbtautau/Analysis/config";
        return FullName(fileName, path);
    }
};

} // namespace mc_corrections
} // namespace analysis
