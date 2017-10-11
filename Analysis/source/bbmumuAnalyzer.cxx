/*! Final analysis step to estimate correction for DY normalization using the muMu channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/BaseEventAnalyzer.h"

namespace analysis {

class bbmumuAnalyzer : public BaseEventAnalyzer<MuonCandidate, MuonCandidate> {
public:
    using Base = BaseEventAnalyzer<MuonCandidate, MuonCandidate>;
    using Base::BaseEventAnalyzer;

protected:
    virtual EventRegion DetermineEventRegion(EventInfo& event, EventCategory /*eventCategory*/) override
    {
        static const std::vector<std::string> trigger_patterns = {
            "HLT_IsoMu22_v"
        };

        const MuonCandidate& muon1 = event.GetFirstLeg();
        const MuonCandidate& muon2 = event.GetSecondLeg();
        const JetCandidate& jet1  = event.GetJets().at(0);
        const JetCandidate& jet2 = event.GetJets().at(1);

        if(!event.GetTriggerResults().AnyAcceptAndMatch(trigger_patterns)) return EventRegion::Unknown();

        const bool os = !ana_setup.apply_os_cut || muon1.GetCharge() * muon2.GetCharge() == -1;
        const bool iso = true;
        double mass_muMu = (muon1.GetMomentum()+muon2.GetMomentum()).M();
        double mass_jj = (jet1.GetMomentum()+jet2.GetMomentum()).M();
        const bool jetMass = (mass_jj > 80 && mass_jj < 160);
        const bool muonMass= (mass_muMu > 60);
        return EventRegion(os, iso);
    }

    virtual const EventRegionSet& EventRegionsToProcess() const override
    {
        static const EventRegionSet regions = {
            EventRegion::OS_Isolated()
        };
        return regions;
    }

    virtual const EventSubCategorySet& EventSubCategoriesToProcess() const override
    {
        static const EventSubCategorySet sub_categories = {
            EventSubCategory().SetCutResult(SelectionCut::InsideMassWindow, true)
        };
        return sub_categories;
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbmumuAnalyzer, analysis::AnalyzerArguments)
