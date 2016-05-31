/*! Analyze flat-tree for etau channel for Htautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/SemileptonicFlatTreeAnalyzer.h"

namespace analysis {

class bbetauAnalyzer : public SemileptonicFlatTreeAnalyzer {
public:
    bbetauAnalyzer(const AnalyzerArguments& _args)
         : SemileptonicFlatTreeAnalyzer(_args, ChannelId()) {}

protected:
    virtual Channel ChannelId() const override { return analysis::Channel::ETau; }
    virtual std::string TreeName() const override { return "etau"; }

    virtual EventRegion DetermineEventRegion(const ntuple::Sync& event,
                                                       analysis::EventCategory /*eventCategory*/) override
    {
        using namespace cuts::Htautau_2015::ETau;

        if(event.againstMuonTight3_2 < tauID::againstMuonTight3
                || event.againstElectronVLooseMVA6_2 < tauID::againstElectronVLooseMVA6
                || event.iso_1 >= electronID::pFRelIso
                || event.dilepton_veto
                /*|| (event.extraelec_veto || event.extramuon_veto) */)
            return EventRegion::Unknown;

        const bool os = event.q_1 * event.q_2 == -1;
        const bool iso = event.iso_1 < electronID::pFRelIso;
//        const bool low_mt = event.mt_1 < electronID::mt;
        const bool low_mt = true;


        if(iso && os) return low_mt ? EventRegion::OS_Isolated : EventRegion::OS_Iso_HighMt;
        if(iso && !os) return low_mt ? EventRegion::SS_Isolated : EventRegion::SS_Iso_HighMt;
        if(os) return low_mt ? EventRegion::OS_AntiIsolated : EventRegion::OS_AntiIso_HighMt;
        return low_mt ? EventRegion::SS_AntiIsolated : EventRegion::SS_AntiIso_HighMt;
    }
};

} // namespace analysis

PROGRAM_MAIN(analysis::bbetauAnalyzer, analysis::AnalyzerArguments)
