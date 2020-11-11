/*! Final analysis step to estimate correction for DY normalization using the muMu channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseEventAnalyzer.h"

namespace analysis {

class bbmumuAnalyzer : public BaseEventAnalyzer {
public:
    using HiggsBBCandidate = EventInfo::HiggsBBCandidate;

    bbmumuAnalyzer(const AnalyzerArguments& _args) : BaseEventAnalyzer(_args, Channel::MuMu) {}

protected:
    virtual EventRegion DetermineEventRegion(EventInfo& eventInfoBase, EventCategory /*eventCategory*/) override
    {

        const std::array<const LepCandidate*, 2> muons = { &eventInfoBase.GetLeg(1), &eventInfoBase.GetLeg(2) };
        std::array<EventRegion, 2> regions;
        const bool os = muons.at(0)->GetCharge() * muons.at(1)->GetCharge() == -1;

        for(size_t n = 0; n < muons.size(); ++n) {
            if(!SetRegionIsoRange(*muons.at(n), regions.at(n)))
                return EventRegion::Unknown();
            regions.at(n).SetCharge(os);
        }

        if(regions.at(0).GetLowerIso() >= DiscriminatorWP::Medium) return regions.at(1);
        if(regions.at(1).GetLowerIso() >= DiscriminatorWP::Medium) return regions.at(0);
        return EventRegion::Unknown();
    }

    virtual EventSubCategory DetermineEventSubCategory(EventInfo& event, const EventCategory& /*category*/,
                                                       std::map<SelectionCut, double>& /*mva_scores*/) /*override*/
    {
        const double mass_muMu = event.GetHiggsTTMomentum(false)->M();
        const double ht_other_jets = event.GetHT(false, true);
        const double pt_mumu =  event.GetHiggsTTMomentum(false)->pt();
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
        sub_category.SetCutResult(SelectionCut::lowHT, ht_other_jets <= 20);
        sub_category.SetCutResult(SelectionCut::medHT, ht_other_jets > 20 && ht_other_jets <= 250);
        sub_category.SetCutResult(SelectionCut::highHT, ht_other_jets > 250);
        sub_category.SetCutResult(SelectionCut::vlowPtNLO, pt_mumu <= 20);
        sub_category.SetCutResult(SelectionCut::lowPtNLO, pt_mumu > 20 && pt_mumu <= 40);
        sub_category.SetCutResult(SelectionCut::medPtNLO, pt_mumu > 40 && pt_mumu <= 100);
        sub_category.SetCutResult(SelectionCut::highPtNLO, pt_mumu > 100);
        sub_category.SetCutResult(SelectionCut::vlowPtLO, ana_setup.pt_sel_bins.size() > 0
                                                          && pt_mumu <= ana_setup.pt_sel_bins.at(0));
        sub_category.SetCutResult(SelectionCut::lowPtLO, ana_setup.pt_sel_bins.size() > 1
                                                         && pt_mumu > ana_setup.pt_sel_bins.at(0)
                                                         && pt_mumu <= ana_setup.pt_sel_bins.at(1));
        sub_category.SetCutResult(SelectionCut::medPt1LO, ana_setup.pt_sel_bins.size() > 2
                                                          && pt_mumu > ana_setup.pt_sel_bins.at(1)
                                                          && pt_mumu <= ana_setup.pt_sel_bins.at(2));
        sub_category.SetCutResult(SelectionCut::medPt2LO, ana_setup.pt_sel_bins.size() > 3
                                                          && pt_mumu > ana_setup.pt_sel_bins.at(2)
                                                          && pt_mumu <= ana_setup.pt_sel_bins.at(3));
        sub_category.SetCutResult(SelectionCut::highPtLO, ana_setup.pt_sel_bins.size() > 4
                                                          && pt_mumu > ana_setup.pt_sel_bins.at(3)
                                                          && pt_mumu <= ana_setup.pt_sel_bins.at(4));
        sub_category.SetCutResult(SelectionCut::vhighPtLO, ana_setup.pt_sel_bins.size() > 4
                                                           && pt_mumu > ana_setup.pt_sel_bins.at(4));
        sub_category.SetCutResult(SelectionCut::mtt, mass_muMu > 50);

        return sub_category;
    }

};

} // namespace analysis

PROGRAM_MAIN(analysis::bbmumuAnalyzer, analysis::AnalyzerArguments)
