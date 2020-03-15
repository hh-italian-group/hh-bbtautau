/*! Final analysis step for the tauTau channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseEventAnalyzer.h"

namespace analysis {

class bbtautauAnalyzer : public BaseEventAnalyzer {
public:

    bbtautauAnalyzer(const AnalyzerArguments& _args) : BaseEventAnalyzer(_args, Channel::TauTau) {}

protected:
    virtual EventRegion DetermineEventRegion(EventInfo& eventInfoBase, EventCategory /*eventCategory*/) override
    {
        const std::array<const LepCandidate*, 2> taus = { &eventInfoBase.GetLeg(1), &eventInfoBase.GetLeg(2) };
        std::array<EventRegion, 2> regions;
        const bool os = taus.at(0)->GetCharge() * taus.at(1)->GetCharge() == -1;

        for(size_t n = 0; n < taus.size(); ++n) {
            if(!SetRegionIsoRange(*taus.at(n), regions.at(n)))
                return EventRegion::Unknown();
            regions.at(n).SetCharge(os);
        }

        const auto signal_wp = signalObjectSelector.GetTauVSjetDiscriminator().second;
        if(regions.at(0).GetLowerIso() >= signal_wp) return regions.at(1);
        if(regions.at(1).GetLowerIso() >= signal_wp) return regions.at(0);
        return EventRegion::Unknown();
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbtautauAnalyzer, analysis::AnalyzerArguments)
