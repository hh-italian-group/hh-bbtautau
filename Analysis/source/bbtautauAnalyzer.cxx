/*! Final analysis step for the tauTau channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseEventAnalyzer.h"
#include "h-tautau/McCorrections/include/TauIdWeight.h"

namespace analysis {

class bbtautauAnalyzer : public BaseEventAnalyzer<TauCandidate, TauCandidate> {
public:
    using Base = BaseEventAnalyzer<TauCandidate, TauCandidate>;
    using Base::BaseEventAnalyzer;

    bbtautauAnalyzer(const AnalyzerArguments& _args) :
        Base(_args), gen(174296), dm_prob_distr(0, 1),
        id_weight_calc("h-tautau/McCorrections/data/fitresults_tt_moriond2017.json", DiscriminatorWP::Medium)
    {
    }

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

        EventRegion region_tau1, region_tau2;

        if(!event.GetTriggerResults().AnyAcceptAndMatch(trigger_patterns)) return EventRegion::Unknown();        

        CorrectIdWeight(event);
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

private:
    void CorrectIdWeight(EventInfo& eventInfo)
    {
        auto& event = *const_cast<ntuple::Event*>(&(*eventInfo));
        event.decayMode_1 = GenDecayMode();
        event.decayMode_2 = GenDecayMode();
        const double id_weight = id_weight_calc.Get(event);
        event.weight_total *= id_weight;
    }

    int GenDecayMode()
    {
        static const std::vector<double> dm_prob = { 0.1151, 0.37, 0.1521 };
        static const double total_prob = std::accumulate(dm_prob.begin(), dm_prob.end(), 0.);
        static const std::vector<int> decay_modes = { 0, 1, 10 };
        double prob = dm_prob_distr(gen) * total_prob;
        for(size_t n = 0; n < dm_prob.size(); ++n) {
            if(prob <= dm_prob.at(n))
                return decay_modes.at(n);
            prob -= dm_prob.at(n);
        }
        throw exception("Unable to generate DM.");
    }

private:
    std::mt19937_64 gen;
    std::uniform_real_distribution<double> dm_prob_distr;
    mc_corrections::TauIdWeight id_weight_calc;
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbtautauAnalyzer, analysis::AnalyzerArguments)
