/*! Final analysis step for the muTau channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseEventAnalyzer.h"

namespace analysis {

class bbmutauAnalyzer : public BaseEventAnalyzer<MuonCandidate, TauCandidate> {
public:
    using Base = BaseEventAnalyzer<MuonCandidate, TauCandidate>;
    using Base::BaseEventAnalyzer;

protected:
    virtual EventRegion DetermineEventRegion(EventInfo& event, EventCategory /*eventCategory*/) override
    {
        const MuonCandidate& muon = event.GetFirstLeg();
        const TauCandidate& tau = event.GetSecondLeg();

        const bool os = !ana_setup.apply_os_cut || muon.GetCharge() * tau.GetCharge() == -1;
        const bool iso = !ana_setup.apply_iso_cut || tau->byIsolationMVA(DiscriminatorWP::Medium);
        return EventRegion(os, iso);
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbmutauAnalyzer, analysis::AnalyzerArguments)
