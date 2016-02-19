/*! Study of possible issues in event selection.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/LightBaseFlatTreeAnalyzer.h"

class EventSelectionStudy : public analysis::LightBaseFlatTreeAnalyzer {
public:
    EventSelectionStudy(const std::string& _inputFileName, const std::string& _outputFileName)
         : LightBaseFlatTreeAnalyzer(_inputFileName,_outputFileName)
    {
        recalc_kinfit = false;
    }

protected:

    virtual void AnalyzeEvent(const analysis::FlatEventInfo& eventInfo, analysis::EventCategory category) override
    {
        using analysis::EventCategory;
        if(DetermineEventRegion(eventInfo,category) != analysis::EventRegion::OS_Isolated) return;
        //if (!analysis::TwoJetsEventCategories_MediumBjets.count(category)) return;

        if(eventInfo.eventEnergyScale == analysis::EventEnergyScale::Central &&
                category == EventCategory::TwoJets_OneBtag &&
                eventInfo.eventType == ntuple::EventType::ZJ )
            std::cout << /*eventInfo.eventType << " " <<*/ eventInfo.event->evt << std::endl;
    }

};
