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
    EventWeights_HH(Period period, DiscriminatorWP btag_wp) :
        EventWeights(period, btag_wp)
    {
        if(period != Period::Run2016)
            throw exception("Period %1% is not supported.") % period;
        providers[WeightType::DY] = std::make_shared<NJets_HT_weight>("DY", Full_Cfg_Name("dyjets_weights.cfg"));
        providers[WeightType::TTbar] = std::shared_ptr<TTbar_weight_multiChannel>(new TTbar_weight_multiChannel
                    { { Channel::ETau, Full_Cfg_Name("ttbar_weights_eTau.cfg") },
                      { Channel::MuTau, Full_Cfg_Name("ttbar_weights_muTau.cfg") },
                      { Channel::TauTau, Full_Cfg_Name("ttbar_weights_tauTau.cfg") } });
        providers[WeightType::BSM_to_SM] = std::make_shared<HH_BMStoSM_weight>(
                    FullBSMtoSM_Name("weight_SM.root"), "weight_node_BSM");
        providers[WeightType::Wjets] = std::make_shared<NJets_HT_weight>("Wjets", Full_Cfg_Name("wjets_weights.cfg"));
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
