/*! Final analysis step to estimate correction for DY normalization using the muMu channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseEventAnalyzer.h"

namespace analysis {

class bbmumuAnalyzer : public BaseEventAnalyzer<MuonCandidate, MuonCandidate> {
public:

    using Base = BaseEventAnalyzer<MuonCandidate, MuonCandidate>;
    using Base::BaseEventAnalyzer;
    using HiggsBBCandidate = EventInfoBase::HiggsBBCandidate;
protected:
    virtual EventRegion DetermineEventRegion(EventInfo& event, EventCategory /*eventCategory*/) override
    {
        static const std::map<DiscriminatorWP,double> working_points = {
            {DiscriminatorWP::VVLoose,2.0}, {DiscriminatorWP::Medium,0.15}
        };

        const MuonCandidate& muon1 = event.GetFirstLeg();
        const MuonCandidate& muon2 = event.GetSecondLeg();
//        const HiggsBBCandidate& jets  = event.GetHiggsBB();

        EventRegion region_muon1, region_muon2;
        const bool os = !ana_setup.apply_os_cut || muon1.GetCharge() * muon2.GetCharge() == -1;
        region_muon1.SetCharge(os);
        region_muon2.SetCharge(os);

        for(auto wp = working_points.rbegin(); wp != working_points.rend(); ++wp) {
            if(muon1.GetIsolation() < wp->second) {
                region_muon1.SetLowerIso(wp->first);
                if(wp != working_points.rbegin())
                    region_muon1.SetUpperIso((--wp)->first);
                break;
            }
        }

        for(auto wp = working_points.rbegin(); wp != working_points.rend(); ++wp) {
            if(muon2.GetIsolation() < wp->second) {
                region_muon2.SetLowerIso(wp->first);
                if(wp != working_points.rbegin())
                    region_muon2.SetUpperIso((--wp)->first);
                break;
            }
        }

        //region.SetLowerIso(DiscriminatorWP::Medium);
        bool muon1_ex=false;
        bool muon2_ex=false;
        if(muon1.GetIsolation() > 0.15 && muon1.GetIsolation() < 0.5) muon1_ex=true;
        if(muon2.GetIsolation() > 0.15 && muon2.GetIsolation() < 0.5) muon2_ex=true;
        if(muon1_ex || muon2_ex) return EventRegion::Unknown();

        if(!region_muon1.HasLowerIso() || !region_muon2.HasLowerIso()) return EventRegion::Unknown();
        if(region_muon1.GetLowerIso() >= DiscriminatorWP::Medium) return region_muon2;
        if(region_muon2.GetLowerIso() >= DiscriminatorWP::Medium) return region_muon1;
        return EventRegion::Unknown();

        //return region;
    }

    virtual EventSubCategory DetermineEventSubCategory(EventInfo& event, const EventCategory& /*category*/,
                                                       std::map<SelectionCut, double>& /*mva_scores*/) override
    {
        double mass_muMu = event.GetHiggsTTMomentum(false).M();
        //const bool jetMass = mass_jj > 80 && mass_jj < 160;
        //const bool muonMass= mass_muMu > 60;
        /*double pt_jets = event.GetHiggsBB().GetFirstDaughter().GetMomentum().Pt()
                        + event.GetHiggsBB().GetSecondDaughter().GetMomentum().Pt()
                        + event->ht_other_jets;*/



        EventSubCategory sub_category;
        if(event.HasBjetPair()){
            double mass_jj = event.GetHiggsBB().GetMomentum().M();
            sub_category.SetCutResult(SelectionCut::mh, ana_setup.massWindowParams.at(SelectionCut::mh)
                              .IsInside(mass_muMu,mass_jj));
        }
        else sub_category.SetCutResult(SelectionCut::mh, false);
        sub_category.SetCutResult(SelectionCut::lowMET, event.GetMET().GetMomentum().Pt() < 45);
        sub_category.SetCutResult(SelectionCut::lowHT, event->ht_other_jets <= 20);
        sub_category.SetCutResult(SelectionCut::medHT,  event->ht_other_jets > 20 && event->ht_other_jets <= 250);
        sub_category.SetCutResult(SelectionCut::highHT, event->ht_other_jets > 250);
        sub_category.SetCutResult(SelectionCut::vlowPt, event.GetHiggsTTMomentum(false).Pt() <= 20);
        sub_category.SetCutResult(SelectionCut::lowPt,  event.GetHiggsTTMomentum(false).Pt() > 20 &&
                                                        event.GetHiggsTTMomentum(false).Pt() <= 40);
        sub_category.SetCutResult(SelectionCut::medPt, event.GetHiggsTTMomentum(false).Pt() > 40 &&
                                                        event.GetHiggsTTMomentum(false).Pt() <= 100);
        sub_category.SetCutResult(SelectionCut::highPt, event.GetHiggsTTMomentum(false).Pt() > 100);
        sub_category.SetCutResult(SelectionCut::mtt, mass_muMu > 60);

        return sub_category;
    }

};

} // namespace analysis

PROGRAM_MAIN(analysis::bbmumuAnalyzer, analysis::AnalyzerArguments)
