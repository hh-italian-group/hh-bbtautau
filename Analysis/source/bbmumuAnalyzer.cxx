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
        static const std::vector<std::string> trigger_patterns = {
            "HLT_IsoMu22_v"
        };

        const MuonCandidate& muon1 = event.GetFirstLeg();
        const MuonCandidate& muon2 = event.GetSecondLeg();
        const HiggsBBCandidate& jets  = event.GetHiggsBB();

        if(!event.GetTriggerResults().AnyAcceptAndMatch(trigger_patterns)) return EventRegion::Unknown();

        EventRegion region;
        const bool os = !ana_setup.apply_os_cut || muon1.GetCharge() * muon2.GetCharge() == -1;
        region.SetCharge(os);
        region.SetLowerIso(DiscriminatorWP::Medium);


        double mass_muMu = (muon1.GetMomentum()+muon2.GetMomentum()).M();
        double mass_jj = jets.GetMomentum().M();
        const bool jetMass = (mass_jj > 80 && mass_jj < 160);
        const bool muonMass= (mass_muMu > 60);
        if(!jetMass || !muonMass) return EventRegion::Unknown();
        if(event.GetMET().GetMomentum().Pt() > 45 ) return EventRegion::Unknown();

        return region;
    }

    virtual const EventCategorySet& EventCategoriesToProcess() const
        {
            static const EventCategorySet categories = {
                EventCategory::TwoJets_ZeroBtag(),
                EventCategory::TwoJets_OneBtag(), /*EventCategory::TwoJets_OneLooseBtag(),*/
                EventCategory::TwoJets_TwoBtag() /*EventCategory::TwoJets_TwoLooseBtag()*/
            };
            return categories;
        }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbmumuAnalyzer, analysis::AnalyzerArguments)
