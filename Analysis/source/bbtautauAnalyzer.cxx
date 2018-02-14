/*! Final analysis step for the tauTau channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseEventAnalyzer.h"
#include "h-tautau/McCorrections/include/TauIdWeight.h"

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

        static const std::vector<DiscriminatorWP> working_points = {
            DiscriminatorWP::VLoose, DiscriminatorWP::Loose, DiscriminatorWP::Medium
        };

        const TauCandidate& tau_1 = event.GetFirstLeg();
        const TauCandidate& tau_2 = event.GetSecondLeg();

//        EventRegion region;

//        if(!event.GetTriggerResults().AnyAcceptAndMatch(trigger_patterns)
//                || !tau_1->byIsolationMVA(DiscriminatorWP::Medium)) return EventRegion::Unknown();

//        const bool os = !ana_setup.apply_os_cut || tau_1.GetCharge() * tau_2.GetCharge() == -1;
//        region.SetCharge(os);

//        for(auto wp = working_points.rbegin(); wp != working_points.rend(); ++wp) {
//            if(tau_2->byIsolationMVA(*wp)) {
//                region.SetLowerIso(*wp);
//                if(wp != working_points.rbegin())
//                    region.SetUpperIso(*(--wp));
//                break;
//            }
//        }

//        return region;

        EventRegion region_tau1, region_tau2;

        if(!event.GetTriggerResults().AnyAcceptAndMatch(trigger_patterns)) return EventRegion::Unknown();

        const bool os = !ana_setup.apply_os_cut || tau_1.GetCharge() * tau_2.GetCharge() == -1;
        region_tau1.SetCharge(os);
        region_tau2.SetCharge(os);

        for(auto wp_1 = working_points.rbegin(); wp_1 != working_points.rend(); ++wp_1) {
            if(tau_1->byIsolationMVA(*wp_1)) {
                region_tau1.SetLowerIso(*wp_1);
                if(wp_1 != working_points.rbegin())
                    region_tau1.SetUpperIso(*(--wp_1));
                break;
            }
        }

        for(auto wp_2 = working_points.rbegin(); wp_2 != working_points.rend(); ++wp_2) {
            if(tau_2->byIsolationMVA(*wp_2)) {
                region_tau2.SetLowerIso(*wp_2);
                if(wp_2 != working_points.rbegin())
                    region_tau2.SetUpperIso(*(--wp_2));
                break;
            }
        }

        if(!region_tau1.HasLowerIso() || !region_tau2.HasLowerIso()) return EventRegion::Unknown();
        if(region_tau1.GetLowerIso() >= DiscriminatorWP::Medium) return region_tau2;
        if(region_tau2.GetLowerIso() >= DiscriminatorWP::Medium) return region_tau1;
        return EventRegion::Unknown();
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbtautauAnalyzer, analysis::AnalyzerArguments)
