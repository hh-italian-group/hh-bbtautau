/*! Final analysis step for the muTau channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseEventAnalyzer.h"

namespace analysis {

class bbmutauAnalyzer : public BaseEventAnalyzer {
public:
    bbmutauAnalyzer(const AnalyzerArguments& _args) : BaseEventAnalyzer(_args, Channel::MuTau) {}

protected:
    virtual EventRegion DetermineEventRegion(EventInfoBase& eventInfoBase, EventCategory /*eventCategory*/) override
    {
        //static const std::vector<DiscriminatorWP> working_points = {
        //    DiscriminatorWP::VLoose, DiscriminatorWP::Loose, DiscriminatorWP::Medium
        //};

	static const std::map<DiscriminatorWP,double> working_points = {
            {DiscriminatorWP::VVLoose,2.0}, {DiscriminatorWP::Medium,0.15}
        };

        const LepCandidate& muon = eventInfoBase.GetFirstLeg();
        const LepCandidate& tau = eventInfoBase.GetSecondLeg();
        
        EventRegion region;

        if(ana_setup.mode == SignalMode::HTT || ana_setup.mode == SignalMode::HTT_sync){
            double mt = analysis::Calculate_MT(muon.GetMomentum(),eventInfoBase.GetMET().GetMomentum());
            if(mt >= cuts::H_tautau_2016::mt) return EventRegion::Unknown();
	    double m_tt_vis = (muon.GetMomentum() + tau.GetMomentum()).M();
	    if(m_tt_vis <= 50) return EventRegion::Unknown();
	    //cut for tau ID SF compatibility
	    //double pzeta = analysis::Calculate_Pzeta(muon.GetMomentum(),tau.GetMomentum(),eventInfoBase.GetMET().GetMomentum());
	    //if(pzeta <= -25) return EventRegion::Unknown();
        }

        const bool os = !ana_setup.apply_os_cut || muon.GetCharge() * tau.GetCharge() == -1;
        region.SetCharge(os);

	if(!tau->Passed(ana_setup.tauIdDiscriminator, ana_setup.workingPoint)) return EventRegion::Unknown();

        //for(auto wp = working_points.rbegin(); wp != working_points.rend(); ++wp) {
          //  if(tau->Passed(TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017, *wp)) {
            //    region.SetLowerIso(*wp);
            //    if(wp != working_points.rbegin())
            //        region.SetUpperIso(*(--wp));
            //    break;
            //}
        //}

        //return region;

	for(auto wp = working_points.rbegin(); wp != working_points.rend(); ++wp) {
            if(muon.GetIsolation() < wp->second) {
                region.SetLowerIso(wp->first);
                if(wp != working_points.rbegin())
                    region.SetUpperIso((--wp)->first);
                break;
            }
        }

        if(muon.GetIsolation() > 0.15 && muon.GetIsolation() < 0.3) return EventRegion::Unknown();
        return region;
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbmutauAnalyzer, analysis::AnalyzerArguments)
