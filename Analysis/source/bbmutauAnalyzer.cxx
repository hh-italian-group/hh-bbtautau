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
          static const std::vector<DiscriminatorWP> working_points_tau = {
             DiscriminatorWP::VLoose, DiscriminatorWP::Loose, DiscriminatorWP::Medium
          };

	       static const std::map<DiscriminatorWP,double> working_points_mu = {
              {DiscriminatorWP::VVLoose,2.0}, {DiscriminatorWP::Medium,0.15}
            };

        const LepCandidate& muon = eventInfoBase.GetFirstLeg();
        const LepCandidate& tau = eventInfoBase.GetSecondLeg();

        EventRegion region;

      	if(ana_setup.mode != SignalMode::HH && ana_setup.mode != SignalMode::HTT_sync){
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

      	TauIdDiscriminator tauIdDiscriminator = TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017;
      	if(ana_setup.mode == SignalMode::TauPOG_deepTauVsJet || ana_setup.mode == SignalMode::TauPOG_deepTauVsJet_full)
      		tauIdDiscriminator = TauIdDiscriminator::byDeepTau2017v2p1VSjet;

        if(ana_setup.qcd_method == QCDmethod::invert_muon){
          if(!tau->Passed(tauIdDiscriminator, ana_setup.tauID_wp)) return EventRegion::Unknown();

          for(auto wp = working_points_mu.rbegin(); wp != working_points_mu.rend(); ++wp) {
             if(muon.GetIsolation() < wp->second) {
                 region.SetLowerIso(wp->first);
                 if(wp != working_points_mu.rbegin())
                     region.SetUpperIso((--wp)->first);
                 break;
             }
           }

           if(muon.GetIsolation() > 0.15 && muon.GetIsolation() < 0.3) return EventRegion::Unknown();
        }
        else {
          if(muon.GetIsolation() >= cuts::hh_bbtautau_2017::MuTau::muonID::pfRelIso04) return EventRegion::Unknown();

          for(auto wp = working_points_tau.rbegin(); wp != working_points_tau.rend(); ++wp) {
             if(tau->Passed(tauIdDiscriminator, *wp)) {
                 region.SetLowerIso(*wp);
                 if(wp != working_points_tau.rbegin())
                     region.SetUpperIso(*(--wp));
                 break;
              }
          }
        }

        return region;
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbmutauAnalyzer, analysis::AnalyzerArguments)
