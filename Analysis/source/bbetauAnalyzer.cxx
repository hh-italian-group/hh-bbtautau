/*! Analyze flat-tree for etau channel for Htautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/SemileptonicEventAnalyzer.h"
#include "Analysis/include/BDT_reader_etau.h"

namespace analysis {

class bbetauAnalyzer : public SemileptonicFlatTreeAnalyzer<ElectronCandidate> {
public:
    //using SemileptonicFlatTreeAnalyzer<ElectronCandidate>::SemileptonicFlatTreeAnalyzer;
	bbetauAnalyzer(const AnalyzerArguments& _args) : SemileptonicFlatTreeAnalyzer(_args), MVA_reader(_args.weight_file()) {}

protected:
    virtual std::string TreeName() const override { return "eTau"; }

    virtual EventRegion DetermineEventRegion(EventInfo& event, EventCategory /*eventCategory*/) override
    {
        using namespace cuts::Htautau_2015::ETau;

        const ElectronCandidate& electron = event.GetFirstLeg();
        const TauCandidate& tau = event.GetSecondLeg();
		double BDT_wp =args.BDT_cut(); /*Francesco*/
		double BDT_output = MVA_reader.BDT_score(event);

        if(tau->againstElectronMVA6(DiscriminatorWP::Tight) < 0.5
                || tau->againstMuon3(DiscriminatorWP::Loose) < 0.5
                || electron->iso() >= 0.15
                || event->extraelec_veto || event->extramuon_veto
				|| BDT_output < BDT_wp /*Francesco*/)
            return EventRegion::Unknown;

        const bool os = electron.GetCharge() * tau.GetCharge() == -1;
        const bool iso = tau->iso() > 0.2;
        const bool low_mt = true;

        if(iso && os) return low_mt ? EventRegion::OS_Isolated : EventRegion::OS_Iso_HighMt;
        if(iso && !os) return low_mt ? EventRegion::SS_Isolated : EventRegion::SS_Iso_HighMt;
        if(os) return low_mt ? EventRegion::OS_AntiIsolated : EventRegion::OS_AntiIso_HighMt;
        return low_mt ? EventRegion::SS_AntiIsolated : EventRegion::SS_AntiIso_HighMt;
    }
	
private:
	BDT_reader MVA_reader;
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbetauAnalyzer, analysis::AnalyzerArguments)
