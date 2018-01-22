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
//        const HiggsBBCandidate& jets  = event.GetHiggsBB();

        if(!event.GetTriggerResults().AnyAcceptAndMatch(trigger_patterns)) return EventRegion::Unknown();

        EventRegion region;
        const bool os = !ana_setup.apply_os_cut || muon1.GetCharge() * muon2.GetCharge() == -1;
        region.SetCharge(os);
        region.SetLowerIso(DiscriminatorWP::Medium);

        return region;
    }

    virtual EventSubCategory DetermineEventSubCategory(EventInfo& event, const EventCategory& /*category*/,
                                                       std::map<SelectionCut, double>& /*mva_scores*/) override
    {
        double mass_muMu = event.GetHiggsTTMomentum(false).M();
        double mass_jj = event.GetHiggsBB().GetMomentum().M();
        const bool jetMass = mass_jj > 80 && mass_jj < 160;
        const bool muonMass= mass_muMu > 60;

        EventSubCategory sub_category;
        sub_category.SetCutResult(SelectionCut::mh, jetMass && muonMass);
        sub_category.SetCutResult(SelectionCut::lowMET,event.GetMET().GetMomentum().Pt() < 45);

        return sub_category;
    }

};

} // namespace analysis

PROGRAM_MAIN(analysis::bbmumuAnalyzer, analysis::AnalyzerArguments)
