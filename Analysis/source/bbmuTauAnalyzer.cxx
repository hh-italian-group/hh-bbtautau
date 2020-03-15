/*! Final analysis step for the muTau channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseEventAnalyzer.h"

namespace analysis {

class bbmutauAnalyzer : public BaseEventAnalyzer {
public:
    bbmutauAnalyzer(const AnalyzerArguments& _args) : BaseEventAnalyzer(_args, Channel::MuTau) {}

protected:
    virtual EventRegion DetermineEventRegion(EventInfo& eventInfoBase, EventCategory /*eventCategory*/) override
    {
        const LepCandidate& muon = eventInfoBase.GetFirstLeg();
        const LepCandidate& tau = eventInfoBase.GetSecondLeg();

        EventRegion region;

      	if(ana_setup.mode == SignalMode::HTT) {
            const double mt = analysis::Calculate_MT(muon.GetMomentum(),eventInfoBase.GetMET().GetMomentum());
            if(mt >= cuts::H_tautau_Run2::MuTau::mt) return EventRegion::Unknown();
        }

        const bool os = muon.GetCharge() * tau.GetCharge() == -1;
        region.SetCharge(os);

        if(ana_setup.qcd_method == QCDmethod::invert_muon) {
            const auto& [tau_discr, tau_wp] = signalObjectSelector.GetTauVSjetDiscriminator();
            if(!tau->Passed(tau_discr, tau_wp) || !SetRegionIsoRange(muon, region))
                return EventRegion::Unknown();
        } else {
            if(muon.GetIsolation() >= cuts::hh_bbtautau_Run2::MuTau::muonID::pfRelIso04
                    || !SetRegionIsoRange(tau, region))
                return EventRegion::Unknown();
        }

        return region;
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbmutauAnalyzer, analysis::AnalyzerArguments)
