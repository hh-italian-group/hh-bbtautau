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


        static const std::vector<std::string> trigger_patterns = ana_setup.trigger.at(ChannelId());

        static const std::vector<DiscriminatorWP> working_points = {
            DiscriminatorWP::VLoose, DiscriminatorWP::Loose, DiscriminatorWP::Medium
        };

        const ElectronCandidate& electron = event.GetFirstLeg();
        const TauCandidate& tau = event.GetSecondLeg();

        EventRegion region;

        if(!event.GetTriggerResults().AnyAcceptAndMatch(trigger_patterns)) return region;

        const bool os = !ana_setup.apply_os_cut || electron.GetCharge() * tau.GetCharge() == -1;
        region.SetCharge(os);

        for(auto wp = working_points.rbegin(); wp != working_points.rend(); ++wp) {
            if(tau->tauID(TauIdDiscriminator::byIsolationMVArun2v1DBoldDMwLT, *wp)) {
                region.SetLowerIso(*wp);
                if(wp != working_points.rbegin())
                    region.SetUpperIso(*(--wp));
                break;
            }
        }

        return region;
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbetauAnalyzer, analysis::AnalyzerArguments)
