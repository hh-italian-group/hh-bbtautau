/*! Final analysis step for the eTau channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseEventAnalyzer.h"

namespace analysis {

class bbetauAnalyzer : public BaseEventAnalyzer<ElectronCandidate, TauCandidate> {
public:
    using Base = BaseEventAnalyzer<ElectronCandidate, TauCandidate>;
    using Base::BaseEventAnalyzer;

protected:
    virtual EventRegion DetermineEventRegion(EventInfo& event, EventCategory /*eventCategory*/) override
    {
        static const std::vector<std::string> trigger_patterns = {
            "HLT_Ele25_eta2p1_WPTight_Gsf_v"
        };

        const ElectronCandidate& electron = event.GetFirstLeg();
        const TauCandidate& tau = event.GetSecondLeg();

        if(!event.GetTriggerResults().AnyAcceptAndMatch(trigger_patterns)) return EventRegion::Unknown();

        const bool os = !ana_setup.apply_os_cut || electron.GetCharge() * tau.GetCharge() == -1;
        const bool iso = !ana_setup.apply_iso_cut || tau->byIsolationMVA(DiscriminatorWP::Medium);
        return EventRegion(os, iso);
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbetauAnalyzer, analysis::AnalyzerArguments)
