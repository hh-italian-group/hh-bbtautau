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
    using Event = ntuple::Event;
    using ExpressEvent = ntuple::ExpressEvent;
    using HH_BMStoSM_weightPtr = std::shared_ptr<HH_BMStoSM_weight>;
    using TTbar_weightPtr = std::shared_ptr<TTbar_weight>;
    using NJets_HT_weightPtr = std::shared_ptr<NJets_HT_weight>;

    EventWeights_HH(const Channel& channel, Period period, DiscriminatorWP btag_wp) :
        EventWeights(period, btag_wp)
    {
        dy_weight = NJets_HT_weightPtr(new class DY_weight(Full_Cfg_Name("dyjets_weights.cfg")));
        if(channel == Channel::ETau) ttbar_weight = TTbar_weightPtr(new class TTbar_weight(Full_Cfg_Name("ttbar_weights_eTau.cfg")));
        if(channel == Channel::MuTau) ttbar_weight = TTbar_weightPtr(new class TTbar_weight(Full_Cfg_Name("ttbar_weights_muTau.cfg")));
        if(channel == Channel::TauTau) ttbar_weight = TTbar_weightPtr(new class TTbar_weight(Full_Cfg_Name("ttbar_weights_tauTau.cfg")));
        sm_weight = HH_BMStoSM_weightPtr(new class
                                         HH_BMStoSM_weight(FullBSMtoSM_Name("weight_SM.root"),"weight_node_BSM"));
        wjets_weight = NJets_HT_weightPtr(new class Wjets_weight(Full_Cfg_Name("wjets_weights.cfg")));

    }

    template<typename Event>
    double GetBSMtoSMweight(const Event& event) const {return sm_weight ? sm_weight->Get(event) : 1;}
    template<typename Event>
    double GetTTbar_weight(const Event& event) const {return ttbar_weight ? ttbar_weight->Get(event) : 1;}
    template<typename Event>
    double GetDY_weight(const Event& event) const {return dy_weight ? dy_weight->Get(event) : 1;}
    template<typename Event>
    double GetWjets_weight(const Event& event) const {return wjets_weight ? wjets_weight->Get(event) : 1;}

    double GetTotalWeight(const Event& event,  bool apply_btag_weight = false, bool apply_bsm_to_sm_weight = false,
                          bool apply_ttbar_weight = false, bool apply_dy_weight = false, bool apply_wjets_weight = false)
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
        if(apply_wjets_weight)
            weight *= GetWjets_weight(event);
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
    NJets_HT_weightPtr dy_weight;
    TTbar_weightPtr ttbar_weight;
    HH_BMStoSM_weightPtr sm_weight;
    NJets_HT_weightPtr wjets_weight;
};

} // namespace mc_corrections
} // namespace analysis
