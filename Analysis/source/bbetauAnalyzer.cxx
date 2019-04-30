/*! Final analysis step for the eTau channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseEventAnalyzer.h"

namespace analysis {

class bbetauAnalyzer : public BaseEventAnalyzer {
public:
    bbetauAnalyzer(const AnalyzerArguments& _args) : BaseEventAnalyzer(_args, Channel::ETau) {}

protected:
    virtual EventRegion DetermineEventRegion(EventInfoBase& eventInfoBase, EventCategory /*eventCategory*/) override
    {
        static const std::vector<DiscriminatorWP> working_points = {
            DiscriminatorWP::VLoose, DiscriminatorWP::Loose, DiscriminatorWP::Medium
        };

        const LepCandidate& electron = eventInfoBase.GetFirstLeg();
        const LepCandidate& tau = eventInfoBase.GetSecondLeg();

        EventRegion region;

        if(ana_setup.mode == SignalMode::HTT || ana_setup.mode == SignalMode::HTT_sync){
            double mt = analysis::Calculate_MT(electron.GetMomentum(),eventInfoBase.GetMET().GetMomentum());
            if(mt >= cuts::H_tautau_2016::mt) return EventRegion::Unknown();
        }

        const bool os = !ana_setup.apply_os_cut || electron.GetCharge() * tau.GetCharge() == -1;
        region.SetCharge(os);

        for(auto wp = working_points.rbegin(); wp != working_points.rend(); ++wp) {
            if(tau->Passed(TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017, *wp)) {
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
