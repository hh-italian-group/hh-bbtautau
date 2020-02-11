/*! Final analysis step for the tauTau channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseEventAnalyzer.h"

namespace analysis {

class bbtautauAnalyzer : public BaseEventAnalyzer {
public:

    bbtautauAnalyzer(const AnalyzerArguments& _args) : BaseEventAnalyzer(_args, Channel::TauTau) {}

protected:
    virtual EventRegion DetermineEventRegion(EventInfoBase& eventInfoBase, EventCategory /*eventCategory*/) override
    {

        static const std::vector<DiscriminatorWP> working_points = {
            DiscriminatorWP::VLoose, DiscriminatorWP::Loose, DiscriminatorWP::Medium
        };

        const LepCandidate& tau_1 = eventInfoBase.GetFirstLeg();
        const LepCandidate& tau_2 = eventInfoBase.GetSecondLeg();

        if (std::abs(tau_1.GetMomentum().eta())>2.1 || std::abs(tau_2.GetMomentum().eta())>2.1) return EventRegion::Unknown();
        EventRegion region_tau1, region_tau2;

        const bool os = tau_1.GetCharge() * tau_2.GetCharge() == -1;
        region_tau1.SetCharge(os);
        region_tau2.SetCharge(os);

        for(auto wp_1 = working_points.rbegin(); wp_1 != working_points.rend(); ++wp_1) {
            if(tau_1->Passed(signalObjectSelector.GetTauVSjetDiscriminator().first, *wp_1)) {
                region_tau1.SetLowerIso(*wp_1);
                if(wp_1 != working_points.rbegin())
                    region_tau1.SetUpperIso(*(--wp_1));
                break;
            }
        }

        for(auto wp_2 = working_points.rbegin(); wp_2 != working_points.rend(); ++wp_2) {
            if(tau_2->Passed(signalObjectSelector.GetTauVSjetDiscriminator().first, *wp_2)) {
                region_tau2.SetLowerIso(*wp_2);
                if(wp_2 != working_points.rbegin())
                    region_tau2.SetUpperIso(*(--wp_2));
                break;
            }
        }

        if(!region_tau1.HasLowerIso() || !region_tau2.HasLowerIso()) return EventRegion::Unknown();
        if(region_tau1.GetLowerIso() >= signalObjectSelector.GetTauVSjetDiscriminator().second) return region_tau2;
        if(region_tau2.GetLowerIso() >= signalObjectSelector.GetTauVSjetDiscriminator().second) return region_tau1;
        return EventRegion::Unknown();
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbtautauAnalyzer, analysis::AnalyzerArguments)
