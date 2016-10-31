/*! Analyze flat-tree for mu-tau channel for HHbbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/SemileptonicEventAnalyzer.h"
#include "Analysis/include/BDT_reader_mutau.h" /*Francesco*/

namespace analysis {

class bbmutauAnalyzer : public SemileptonicFlatTreeAnalyzer<MuonCandidate> {
public:
    //using SemileptonicFlatTreeAnalyzer<MuonCandidate>::SemileptonicFlatTreeAnalyzer;
	bbmutauAnalyzer(const AnalyzerArguments& _args) : SemileptonicFlatTreeAnalyzer(_args), MVA_reader(_args.weight_file()) {}

protected:
    virtual std::string TreeName() const override { return "muTau"; }

    virtual EventRegion DetermineEventRegion(EventInfo& event, EventCategory /*eventCategory*/) override
    {
        using namespace cuts::Htautau_2015::MuTau;

        const MuonCandidate& muon = event.GetFirstLeg();
        const TauCandidate& tau = event.GetSecondLeg();
		//double BDT_wp = -0.065; /*Francesco*/
		double BDT_wp =args.BDT_cut(); /*Francesco*/
		double BDT_output = MVA_reader.BDT_score(event);
		//std::cout<<"  BDT score: " << BDT_output << " passed: " << (BDT_output > BDT_wp) <<std::endl;

        if(tau->againstMuon3(DiscriminatorWP::Tight) < tauID::againstMuonTight3
                || tau->againstElectronMVA6(DiscriminatorWP::VLoose) < tauID::againstElectronVLooseMVA6
                || muon->iso() >= muonID::pFRelIso
                || event->dilepton_veto
				|| BDT_output < BDT_wp /*Francesco*/
				//|| event->pfmt_1>BDT_wp
                /*|| (event.extraelec_veto || event.extramuon_veto) */)
            return EventRegion::Unknown;

        const bool os = muon.GetCharge() * tau.GetCharge() == -1;
		//const bool iso = event.byTightIsolationMVArun2v1DBoldDMwLT_2 > 0.5;
        const bool iso = tau->iso() > 0.2;
        //const bool low_mt = event.pfmt_1 < muonID::mt;
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

PROGRAM_MAIN(analysis::bbmutauAnalyzer, analysis::AnalyzerArguments)
