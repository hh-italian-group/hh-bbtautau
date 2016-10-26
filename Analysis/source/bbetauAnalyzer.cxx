/*! Analyze flat-tree for etau channel for Htautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/SemileptonicEventAnalyzer.h"

namespace analysis {

class bbetauAnalyzer : public SemileptonicFlatTreeAnalyzer<ElectronCandidate> {
public:
    using SemileptonicFlatTreeAnalyzer<ElectronCandidate>::SemileptonicFlatTreeAnalyzer;

protected:
    virtual std::string TreeName() const override { return "eTau"; }

    virtual EventRegion DetermineEventRegion(EventInfo& event, EventCategory /*eventCategory*/) override
    {
        using namespace cuts::Htautau_2015::ETau;

        const ElectronCandidate& electron = event.GetFirstLeg();
        const TauCandidate& tau = event.GetSecondLeg();

        if(tau->againstElectronMVA6(DiscriminatorWP::Tight) < 0.5
                || tau->againstMuon3(DiscriminatorWP::Loose) < 0.5
                || electron->iso() >= 0.15
                || event->extraelec_veto || event->extramuon_veto)
            return EventRegion::Unknown;

        const bool os = electron.GetCharge() * tau.GetCharge() == -1;
        const bool iso = tau->iso() > 0.2;
        const bool low_mt = true;

        if(iso && os) return low_mt ? EventRegion::OS_Isolated : EventRegion::OS_Iso_HighMt;
        if(iso && !os) return low_mt ? EventRegion::SS_Isolated : EventRegion::SS_Iso_HighMt;
        if(os) return low_mt ? EventRegion::OS_AntiIsolated : EventRegion::OS_AntiIso_HighMt;
        return low_mt ? EventRegion::SS_AntiIsolated : EventRegion::SS_AntiIso_HighMt;
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbetauAnalyzer, analysis::AnalyzerArguments)
