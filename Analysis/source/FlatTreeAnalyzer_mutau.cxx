/*! Analyze flat-tree for mu-tau channel for HHbbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/SemileptonicFlatTreeAnalyzer.h"

class FlatTreeAnalyzer_mutau : public analysis::SemileptonicFlatTreeAnalyzer {
public:
    FlatTreeAnalyzer_mutau(const std::string& source_cfg, const std::string& _inputPath,
                           const std::string& outputFileName, const std::string& signal_list,
                           bool applyPostFitCorrections = false, bool saveFullOutput = false)
         : SemileptonicFlatTreeAnalyzer(analysis::DataCategoryCollection(source_cfg, signal_list, ChannelId()),
                                        _inputPath, outputFileName, applyPostFitCorrections, saveFullOutput)
    {
    }

protected:
    virtual analysis::Channel ChannelId() const override { return analysis::Channel::MuTau; }

    virtual analysis::EventRegion DetermineEventRegion(const ntuple::Flat& event,
                                                       analysis::EventCategory eventCategory) override
    {
        using analysis::EventRegion;
        using namespace cuts::Htautau_Summer13::MuTau;

        if(!event.againstMuonTight_2
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                || (event.pfRelIso_1 >= muonID::pFRelIso && !IsAntiIsolatedRegion(event))
                || (event.mt_1 >= muonID::mt && !IsHighMtRegion(event,eventCategory)) /*|| event.pt_2 <= 30*/ )
            return EventRegion::Unknown;

        const bool os = event.q_1 * event.q_2 == -1;
        const bool iso = event.pfRelIso_1 < muonID::pFRelIso;
        const bool low_mt = event.mt_1 < muonID::mt;


        if(iso && os) return low_mt ? EventRegion::OS_Isolated : EventRegion::OS_Iso_HighMt;
        if(iso && !os) return low_mt ? EventRegion::SS_Isolated : EventRegion::SS_Iso_HighMt;
        if(os) return low_mt ? EventRegion::OS_AntiIsolated : EventRegion::OS_AntiIso_HighMt;
        return low_mt ? EventRegion::SS_AntiIsolated : EventRegion::SS_AntiIso_HighMt;
    }
};
