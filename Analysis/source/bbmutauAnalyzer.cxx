/*! Analyze flat-tree for mu-tau channel for HHbbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/SemileptonicFlatTreeAnalyzer.h"

class bbmutauAnalyzer : public analysis::SemileptonicFlatTreeAnalyzer {
public:
    bbmutauAnalyzer(const analysis::AnalyzerArguments& _args)
         : SemileptonicFlatTreeAnalyzer(_args, ChannelId()) {}

private:
  static bool IsAntiIsolatedRegion(const ntuple::Sync& event)
    {
        return event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 > 1.6 &&
                event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < 10;
    }

protected:
    virtual analysis::Channel ChannelId() const override { return analysis::Channel::MuTau; }

    virtual analysis::EventRegionSet DetermineEventRegion(const ntuple::Sync& event,
                                                       analysis::EventCategory eventCategory) override
    {
        using analysis::EventRegion;
        using namespace cuts::Htautau_2015::MuTau;

         analysis::EventRegionSet tmpSet;

        const bool againstLeptons = event.againstMuonTight3_2 >= tauID::againstMuonTight3 &&
                                    event.againstElectronVLooseMVA6_2 >= tauID::againstElectronVLooseMVA6;

        if( !againstLeptons
                 ||  event.iso_1 >= muonID::pFRelIso|| event.dilepton_veto
              /*|| (event.mt_1 >= muonID::mt && event.mt_1<=70 && event.mt_1>140) || event.pt_2 <= 30*/ ){
            tmpSet.insert(EventRegion::Unknown);
            return tmpSet;
          }

        const bool os = event.q_1 * event.q_2 == -1;
        const bool iso = event.iso_2 > 0.2 ;
        const bool low_mt = event.pfmt_1 < muonID::mt;
        //const bool low_mt = true;

        if(iso && os) {
            tmpSet.insert(EventRegion::OS_Isolated);
            if(!low_mt && event.pfmt_1 > 70) tmpSet.insert(EventRegion::OS_Iso_HighMt);
            if(low_mt) tmpSet.insert(EventRegion::OS_Iso_LowMt);
        }
        if(iso && !os) {
            tmpSet.insert(EventRegion::SS_Isolated);
            if(!low_mt && event.pfmt_1 > 70) tmpSet.insert(EventRegion::SS_Iso_HighMt);
            if(low_mt) tmpSet.insert(EventRegion::SS_Iso_LowMt);
        }
        if(!iso && os) {
            tmpSet.insert(EventRegion::OS_AntiIsolated);
            if(!low_mt && event.pfmt_1 > 70) tmpSet.insert(EventRegion::OS_AntiIso_HighMt);
            if(low_mt) tmpSet.insert(EventRegion::OS_AntiIso_LowMt);
        }
        if(!iso && !os) {
            tmpSet.insert(EventRegion::SS_AntiIsolated);
            if(!low_mt && event.pfmt_1 > 70) tmpSet.insert(EventRegion::SS_AntiIso_HighMt);
            if(low_mt) tmpSet.insert(EventRegion::SS_AntiIso_LowMt);
        }

        return tmpSet;
    }
};

PROGRAM_MAIN(bbmutauAnalyzer, analysis::AnalyzerArguments)
