/*! Analyze flat-tree for mu-tau channel for HHbbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/SemileptonicFlatTreeAnalyzer.h"

namespace analysis {

class bbmutauAnalyzer : public SemileptonicFlatTreeAnalyzer {
public:
    bbmutauAnalyzer(const AnalyzerArguments& _args)
         : SemileptonicFlatTreeAnalyzer(_args, ChannelId()) {}

protected:
    virtual Channel ChannelId() const override { return Channel::MuTau; }

    virtual EventRegion DetermineEventRegion(const ntuple::Sync& event, EventCategory /*eventCategory*/) override
    {
        using namespace cuts::Htautau_2015::MuTau;

        if(event.againstMuonTight3_2 < tauID::againstMuonTight3
                || event.againstElectronVLooseMVA6_2 < tauID::againstElectronVLooseMVA6
                || event.iso_1 >= muonID::pFRelIso
                || event.dilepton_veto)
            return EventRegion::Unknown;

        const bool os = event.q_1 * event.q_2 == -1;
        const bool iso = event.iso_2 > 0.2;
        const bool low_mt = event.pfmt_1 < muonID::mt;

        if(iso && os) return low_mt ? EventRegion::OS_Isolated : EventRegion::OS_Iso_HighMt;
        if(iso && !os) return low_mt ? EventRegion::SS_Isolated : EventRegion::SS_Iso_HighMt;
        if(os) return low_mt ? EventRegion::OS_AntiIsolated : EventRegion::OS_AntiIso_HighMt;
        return low_mt ? EventRegion::SS_AntiIsolated : EventRegion::SS_AntiIso_HighMt;
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbmutauAnalyzer, analysis::AnalyzerArguments)
