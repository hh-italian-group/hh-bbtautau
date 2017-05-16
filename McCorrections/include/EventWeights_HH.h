/*! Various hh event weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/McCorrections/include/EventWeights.h"
#include "HH_BMStoSM_weight.h"
#include "DY_weight.h"
#include "TTbar_weight.h"

namespace analysis {
namespace mc_corrections {

class EventWeights_HH : public EventWeights {
public:
    using Event = ntuple::Event;
    using ExpressEvent = ntuple::ExpressEvent;
    using HH_BMStoSM_weightPtr = std::shared_ptr<HH_BMStoSM_weight>;
    using TTbar_weightPtr = std::shared_ptr<TTbar_weight>;
    using DY_weightPtr = std::shared_ptr<DY_weight>;

    EventWeights_HH(const Channel& channel, Period period, DiscriminatorWP btag_wp) :
        EventWeights(period, btag_wp)
    {
        dy_weight = DY_weightPtr(new class DY_weight(Full_Cfg_Name("dyjets_weights.cfg")));
        if(channel == Channel::ETau) ttbar_weight = TTbar_weightPtr(new class TTbar_weight(Full_Cfg_Name("ttbar_weights_eTau.cfg")));
        if(channel == Channel::MuTau) ttbar_weight = TTbar_weightPtr(new class TTbar_weight(Full_Cfg_Name("ttbar_weights_muTau.cfg")));
        if(channel == Channel::TauTau) ttbar_weight = TTbar_weightPtr(new class TTbar_weight(Full_Cfg_Name("ttbar_weights_tauTau.cfg")));
        sm_weight = HH_BMStoSM_weightPtr(new class
                                         HH_BMStoSM_weight(FullBSMtoSM_Name("weight_SM.root"),"weight_node_BSM"));

    }

    template<typename Event>
    double GetBSMtoSMweight(const Event& event) const {return sm_weight ? sm_weight->Get(event) : 1;}
    template<typename Event>
    double GetTTbar_weight(const Event& event) const {return ttbar_weight ? ttbar_weight->Get(event) : 1;}
    template<typename Event>
    double GetDY_weight(const Event& event) const {return dy_weight ? dy_weight->Get(event) : 1;}

    double GetTotalWeight(const Event& event,  bool apply_btag_weight = false, bool apply_bsm_to_sm_weight = false,
                          bool apply_ttbar_weight = false, bool apply_dy_weight = false)
    {
        double weight = GetPileUpWeight(event) * GetLeptonTotalWeight(event) * GetTopPtWeight(event);
        if(apply_btag_weight)
            weight *= GetBtagWeight(event);
        if(apply_bsm_to_sm_weight)
            weight *= GetBSMtoSMweight(event);
        if(apply_ttbar_weight)
            weight *= GetTTbar_weight(event);
        if(apply_dy_weight)
            weight *= GetDY_weight(event);
        return weight;
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

    static std::string FullName(const std::string& fileName, const std::string& path)
    {
        return path + "/" + fileName;
    }

private:
    DY_weightPtr dy_weight;
    TTbar_weightPtr ttbar_weight;
    HH_BMStoSM_weightPtr sm_weight;
};

} // namespace mc_corrections
} // namespace analysis
