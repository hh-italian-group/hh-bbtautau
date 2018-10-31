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
        static const std::vector<DiscriminatorWP> working_points = {
            DiscriminatorWP::VLoose, DiscriminatorWP::Loose, DiscriminatorWP::Medium
        };


        if (event->run == 297050 && event->lumi == 678 && event->evt == 799196447)
            std::cout<<"    ----Inizio Determine Event Region----     "<<std::endl;

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


        if (std::abs(tau_1.GetMomentum().eta())>2.1 || std::abs(tau_2.GetMomentum().eta())>2.1) return EventRegion::Unknown();
        if (event->run == 297050 && event->lumi == 678 && event->evt == 799196447)
            std::cout<<"    ----Taglio Eta passato----     "<<std::endl;

        EventRegion region_tau1, region_tau2;

        const bool os = !ana_setup.apply_os_cut || tau_1.GetCharge() * tau_2.GetCharge() == -1;
        region_tau1.SetCharge(os);
        region_tau2.SetCharge(os);

        for(auto wp_1 = working_points.rbegin(); wp_1 != working_points.rend(); ++wp_1) {
            if(tau_1->tauID(TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017, *wp_1)) {
                region_tau1.SetLowerIso(*wp_1);
                if(wp_1 != working_points.rbegin())
                    region_tau1.SetUpperIso(*(--wp_1));
                break;
            }
        }

        for(auto wp_2 = working_points.rbegin(); wp_2 != working_points.rend(); ++wp_2) {
            if(tau_2->tauID(TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017, *wp_2)) {
                region_tau2.SetLowerIso(*wp_2);
                if(wp_2 != working_points.rbegin())
                    region_tau2.SetUpperIso(*(--wp_2));
                break;
            }
        }

        if (event->run == 297050 && event->lumi == 678 && event->evt == 799196447)
            std::cout<<"    ----Prima di egion_tau1.HasLowerIso()----     "<<std::endl;
        if(!region_tau1.HasLowerIso() || !region_tau2.HasLowerIso()) return EventRegion::Unknown();


        if (event->run == 297050 && event->lumi == 678 && event->evt == 799196447)
            std::cout<<"    ----Prima di region tau2----     "<<std::endl;
        if(region_tau1.GetLowerIso() >= DiscriminatorWP::Medium) return region_tau2;

        if (event->run == 297050 && event->lumi == 678 && event->evt == 799196447)
            std::cout<<"    ----Prima di region tau1----     "<<std::endl;

        if(region_tau2.GetLowerIso() >= DiscriminatorWP::Medium) return region_tau1;

        if (event->run == 297050 && event->lumi == 678 && event->evt == 799196447)
            std::cout<<"    ----Prima di Region Unkown----     "<<std::endl;

        return EventRegion::Unknown();
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbtautauAnalyzer, analysis::AnalyzerArguments)
