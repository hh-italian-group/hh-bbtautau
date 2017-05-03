/*! Various hh event weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "../../../h-tautau/McCorrections/include/EventWeights.h"
#include "HH_BMStoSM_weight.h"
#include "DY_weight.h"
#include "TTbar_weight.h"

namespace analysis {
namespace mc_corrections {

class EventWeights_HH : public EventWeights {
public:
    using Event = ntuple::Event;
    using HH_BMStoSM_weightPtr = std::shared_ptr<HH_BMStoSM_weight>;

    EventWeights_HH()
    {
        sm_weight = HH_BMStoSM_weightPtr(new class
                                         HH_BMStoSM_weight(FullBSMtoSM_Name("weight_SM.root"),"weight_node_BSM"));

    }

    double GetBSMtoSMweight(const Event& event) const {return sm_weight ? sm_weight->Get(event) : 1;}

    virtual double GetTotalWeight(const Event& event, bool apply_btag_weight = false, bool apply_bsm_to_sm_weight = false) override
    {
        double weight = GetPileUpWeight(event) * GetLeptonTotalWeight(event);
        if(apply_btag_weight)
            weight *= GetBtagWeight(event);
        if(apply_bsm_to_sm_weight)
            weight *= GetBSMtoSMweight(event);
        return weight;
    }

private:
    static std::string FullBSMtoSM_Name(const std::string& fileName)
    {
        static const std::string path = "hh-bbtautau/McCorrections/data";
        return FullName(fileName, path);
    }

private:
    PileUpWeightPtr pileUp;
    LeptonWeightsPtr lepton;
    BTagWeightPtr bTag;
    HH_BMStoSM_weightPtr sm_weight;
};

} // namespace mc_corrections
} // namespace analysis
