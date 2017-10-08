/*! Final analysis step for the tauTau channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseEventAnalyzer.h"

namespace analysis {

class bbtautauAnalyzer : public BaseEventAnalyzer<TauCandidate, TauCandidate> {
public:
    using Base = BaseEventAnalyzer<TauCandidate, TauCandidate>;
    using Base::BaseEventAnalyzer;

protected:
    virtual EventRegion DetermineEventRegion(EventInfo& event, EventCategory /*eventCategory*/) override
    {
        static const std::vector<std::string> trigger_patterns = {
            "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v", "HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v"
        };

        const TauCandidate& tau_1 = event.GetFirstLeg();
        const TauCandidate& tau_2 = event.GetSecondLeg();

        if(!event.GetTriggerResults().AnyAcceptAndMatch(trigger_patterns)) return EventRegion::Unknown();

        if(ana_setup.apply_iso_cut &&
                !(tau_1->byIsolationMVA(DiscriminatorWP::VLoose) && tau_2->byIsolationMVA(DiscriminatorWP::VLoose) &&
                (tau_1->byIsolationMVA(DiscriminatorWP::Medium) || tau_2->byIsolationMVA(DiscriminatorWP::Medium))))
            return EventRegion::Unknown();

        const bool os = !ana_setup.apply_os_cut || tau_1.GetCharge() * tau_2.GetCharge() == -1;
        const bool iso = !ana_setup.apply_iso_cut ||
                ( tau_1->byIsolationMVA(DiscriminatorWP::Medium) && tau_2->byIsolationMVA(DiscriminatorWP::Medium) );
        return EventRegion(os, iso);
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbtautauAnalyzer, analysis::AnalyzerArguments)
