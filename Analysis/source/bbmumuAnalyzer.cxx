/*! Final analysis step to estimate correction for DY normalization using the muMu channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseEventAnalyzer.h"

namespace analysis {

class bbmumuAnalyzer : public BaseEventAnalyzer {
public:
    using EventInfo = ::analysis::EventInfo<MuonCandidate, MuonCandidate>;
    using HiggsBBCandidate = EventInfoBase::HiggsBBCandidate;

    bbmumuAnalyzer(const AnalyzerArguments& _args) : BaseEventAnalyzer(_args, Channel::MuMu) {}

protected:
    virtual EventRegion DetermineEventRegion(EventInfoBase& eventInfoBase, EventCategory /*eventCategory*/) override
    {
        static const std::map<DiscriminatorWP,double> working_points = {
            {DiscriminatorWP::Loose,2.0}, {DiscriminatorWP::Medium,0.15}
        };

        EventInfo& event = *dynamic_cast<EventInfo*>(&eventInfoBase);

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

        if(!region_muon1.HasLowerIso() || !region_muon2.HasLowerIso()) return EventRegion::Unknown();
        if(region_muon1.GetLowerIso() >= DiscriminatorWP::Medium) return region_muon2;
        if(region_muon2.GetLowerIso() >= DiscriminatorWP::Medium) return region_muon1;
        return EventRegion::Unknown();

        //return region;
    }

    virtual EventSubCategory DetermineEventSubCategory(EventInfoBase& event, const EventCategory& /*category*/,
                                                       std::map<SelectionCut, double>& /*mva_scores*/) override
    {
        double mass_muMu = event.GetHiggsTTMomentum(false).M();
        double mass_jj = event.GetHiggsBB().GetMomentum().M();
        //const bool jetMass = mass_jj > 80 && mass_jj < 160;
        //const bool muonMass= mass_muMu > 60;
        /*double pt_jets = event.GetHiggsBB().GetFirstDaughter().GetMomentum().Pt()
                        + event.GetHiggsBB().GetSecondDaughter().GetMomentum().Pt()
                        + event->ht_other_jets;*/



        EventSubCategory sub_category;
        sub_category.SetCutResult(SelectionCut::mh, ana_setup.massWindowParams.at(SelectionCut::mh)
                              .IsInside(mass_muMu,mass_jj));
        sub_category.SetCutResult(SelectionCut::lowMET,event.GetMET().GetMomentum().Pt() < 45);
        sub_category.SetCutResult(SelectionCut::lowHT, event->ht_other_jets <= 20);
        sub_category.SetCutResult(SelectionCut::medHT,event->ht_other_jets > 20 && event->ht_other_jets <= 250);
        sub_category.SetCutResult(SelectionCut::highHT,event->ht_other_jets > 250);

        return sub_category;
    }

};

} // namespace analysis

PROGRAM_MAIN(analysis::bbmumuAnalyzer, analysis::AnalyzerArguments)
