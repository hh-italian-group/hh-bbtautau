/*! Final analysis step for the eTau channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseEventAnalyzer.h"

namespace analysis {

class bbetauAnalyzer : public BaseEventAnalyzer {
public:
    bbetauAnalyzer(const AnalyzerArguments& _args) : BaseEventAnalyzer(_args, Channel::ETau) {}

protected:
    virtual EventRegion DetermineEventRegion(EventInfo& eventInfoBase, EventCategory /*eventCategory*/) override
    {
        const LepCandidate& electron = eventInfoBase.GetFirstLeg();
        const LepCandidate& tau = eventInfoBase.GetSecondLeg();

        EventRegion region;

        if(ana_setup.mode == SignalMode::HTT) {
            const double mt = analysis::Calculate_MT(electron.GetMomentum(), eventInfoBase.GetMET().GetMomentum());
            if(mt >= cuts::H_tautau_Run2::ETau::mt) return EventRegion::Unknown();
        }

        const bool os = electron.GetCharge() * tau.GetCharge() == -1;
        region.SetCharge(os);

        if(ana_setup.qcd_method == QCDmethod::invert_muon) {
            const auto& [tau_discr, tau_wp] = signalObjectSelector.GetTauVSjetDiscriminator();
            if(!tau->Passed(tau_discr, tau_wp) || !SetRegionIsoRange(electron, region))
                return EventRegion::Unknown();
        } else {
            if(!electron->passEleIsoId(DiscriminatorWP::Tight) || !SetRegionIsoRange(tau, region))
                return EventRegion::Unknown();
        }

        return region;
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbetauAnalyzer, analysis::AnalyzerArguments)
