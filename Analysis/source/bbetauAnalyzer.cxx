/*! Analyze flat-tree for etau channel for Htautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/SemileptonicFlatTreeAnalyzer.h"

class FlatTreeAnalyzer_etau : public analysis::SemileptonicFlatTreeAnalyzer {
public:
    FlatTreeAnalyzer_etau(const std::string& source_cfg, const std::string& _inputPath,
                          const std::string& outputFileName, const std::string& signal_list,
                          bool applyPostFitCorrections = false, bool saveFullOutput = false)
        : SemileptonicFlatTreeAnalyzer(analysis::DataCategoryCollection(source_cfg, signal_list, ChannelId()),
                                       _inputPath, outputFileName, applyPostFitCorrections, saveFullOutput)
    {
    }

protected:
    virtual analysis::Channel ChannelId() const override { return analysis::Channel::ETau; }

    virtual analysis::EventRegion DetermineEventRegion(const ntuple::Flat& event,
                                                       analysis::EventCategory eventCategory) override
    {
        using analysis::EventRegion;
        using namespace cuts::Htautau_Summer13::ETau;

        if(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                || (event.pfRelIso_1 >= electronID::pFRelIso && !IsAntiIsolatedRegion(event))
                || (event.mt_1 >= electronID::mt && !IsHighMtRegion(event,eventCategory)) /*|| event.pt_2 <= 30*/ )
            return EventRegion::Unknown;

        const bool os = event.q_1 * event.q_2 == -1;
        const bool iso = event.pfRelIso_1 < electronID::pFRelIso;
        const bool low_mt = event.mt_1 < electronID::mt;


        if(iso && os) return low_mt ? EventRegion::OS_Isolated : EventRegion::OS_Iso_HighMt;
        if(iso && !os) return low_mt ? EventRegion::SS_Isolated : EventRegion::SS_Iso_HighMt;
        if(os) return low_mt ? EventRegion::OS_AntiIsolated : EventRegion::OS_AntiIso_HighMt;
        return low_mt ? EventRegion::SS_AntiIsolated : EventRegion::SS_AntiIso_HighMt;
    }
};
