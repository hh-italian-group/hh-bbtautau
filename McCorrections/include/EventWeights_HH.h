/*! Various hh event weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/McCorrections/include/EventWeights.h"
#include "HH_BMStoSM_weight.h"
#include "TTbar_weight.h"
#include "NJets_HT_weight.h"

namespace analysis {
namespace mc_corrections {

class EventWeights_HH : public EventWeights {
public:
    EventWeights_HH(Period period, DiscriminatorWP btag_wp, bool use_LLR_weights) :
        EventWeights(period, btag_wp)
    {
        if(period != Period::Run2016)
            throw exception("Period %1% is not supported.") % period;
        std::string dy_weights =
                use_LLR_weights ? Full_Cfg_Name("dyjets_weights_LLR.cfg") : Full_Cfg_Name("dyjets_weights.cfg");
        providers[WeightType::DY] = std::make_shared<NJets_HT_weight>("DY", dy_weights);
        providers[WeightType::TTbar] = std::make_shared<TTbar_weight>(Full_Cfg_Name("ttbar_weights_full.cfg"));
        providers[WeightType::BSM_to_SM] = std::make_shared<HH_BMStoSM_weight>(
                    FullBSMtoSM_Name("weight_SM.root"), "weight");
        std::string wjet_weights =
                use_LLR_weights ? Full_Cfg_Name("wjets_weights_LLR.cfg") : Full_Cfg_Name("wjets_weights.cfg");
        providers[WeightType::Wjets] = std::make_shared<NJets_HT_weight>("Wjets", wjet_weights);
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
