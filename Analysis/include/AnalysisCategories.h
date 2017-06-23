/*! Definition of data and event categories used in HH->bbTauTau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <iostream>
#include <list>
#include <map>
#include <set>
#include <cmath>

#include <Rtypes.h>

#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/Tools.h"
//#include "AnalysisTools/Print/include/PlotPrimitives.h"
//#include "AnalysisTools/Print/include/RootPrintTools.h"
#include "h-tautau/Analysis/include/EventInfo.h"

namespace analysis {

enum class DataCategoryType { Signal, Signal_SM, Background, DataDrivenBkg, Data};

ENUM_NAMES(DataCategoryType) = {
    { DataCategoryType::Signal, "Signal" }, { DataCategoryType::Signal_SM, "Signal_SM" },
    { DataCategoryType::Background, "Background" }, { DataCategoryType::DataDrivenBkg, "DataDrivenBkg" },
    { DataCategoryType::Data, "Data" },
};


enum class EventRegion { Unknown = 0, OS_Isolated = 1, OS_AntiIsolated = 2, SS_Isolated = 3, SS_AntiIsolated = 4,
                         OS_Iso_HighMt = 5, SS_Iso_HighMt = 6, OS_AntiIso_HighMt = 7, SS_AntiIso_HighMt = 8 };

enum class EventCategory { Inclusive = 0, TwoJets_Inclusive = 1, TwoJets_ZeroBtag = 2, TwoJets_OneBtag = 3,
                           TwoJets_TwoBtag = 4, TwoJets_ZeroLooseBtag = 5, TwoJets_OneLooseBtag = 6,
                           TwoJets_TwoLooseBtag = 7, TwoJets_AtLeastOneBtag = 8, TwoJets_AtLeastOneLooseBtag = 9 };

enum class EventSubCategory { NoCuts = 0, KinematicFitConverged = 1, MassWindow = 2,
                              KinematicFitConvergedWithMassWindow = 3, OutsideMassWindow = 4,
                              KinematicFitConvergedOutsideMassWindow = 5
                            };

ENUM_NAMES(EventRegion) = {
    { EventRegion::Unknown, "Unknown"}, { EventRegion::OS_Isolated, "OS_Isolated"},
    { EventRegion::OS_AntiIsolated, "OS_AntiIsolated"}, { EventRegion::SS_Isolated, "SS_Isolated"},
    { EventRegion::SS_AntiIsolated, "SS_AntiIsolated"},  { EventRegion::OS_Iso_HighMt, "OS_Iso_HighMt"},
    { EventRegion::SS_Iso_HighMt, "SS_Iso_HighMt"} , { EventRegion::OS_AntiIso_HighMt, "OS_AntiIso_HighMt"},
    { EventRegion::SS_AntiIso_HighMt, "SS_AntiIso_HighMt"}
};

ENUM_NAMES(EventCategory) = {
    { EventCategory::Inclusive, "Inclusive" }, { EventCategory::TwoJets_Inclusive, "2jets" },
    { EventCategory::TwoJets_ZeroBtag, "2jets0btag" }, { EventCategory::TwoJets_OneBtag, "2jets1btag"},
    { EventCategory::TwoJets_TwoBtag, "2jets2btag" },
    { EventCategory::TwoJets_ZeroLooseBtag, "2jets0Loosebtag" },
    { EventCategory::TwoJets_OneLooseBtag, "2jets1Loosebtag" },
    { EventCategory::TwoJets_TwoLooseBtag, "2jets2Loosebtag" },
    { EventCategory::TwoJets_AtLeastOneBtag, "2jets_at_least_1btag" },
    { EventCategory::TwoJets_AtLeastOneLooseBtag, "2jets_at_least_1Loosebtag" }
};

ENUM_NAMES(EventSubCategory) = {
    { EventSubCategory::NoCuts, "NoCuts" }, { EventSubCategory::KinematicFitConverged, "KinFitConverged" },
    { EventSubCategory::MassWindow, "MassWindow" },
    { EventSubCategory::KinematicFitConvergedWithMassWindow, "KinFitConvergedWithMassWindow" },
    { EventSubCategory::KinematicFitConvergedOutsideMassWindow, "KinematicFitConvergedOutsideMassWindow" },
    { EventSubCategory::OutsideMassWindow, "OutsideMassWindow" }
};

using EventCategoryVector = std::vector<EventCategory>;
using EventCategorySet = EnumNameMap<EventCategory>::EnumEntrySet;
using EventCategoryMap = std::map<EventCategory, EventCategory>;
using EventSubCategorySet = EnumNameMap<EventSubCategory>::EnumEntrySet;

static const EventCategorySet AllEventCategories = __EventCategory_names<>::names.GetEnumEntries();
static const EventSubCategorySet AllEventSubCategories = __EventSubCategory_names<>::names.GetEnumEntries();

static const EventCategoryMap MediumToLoose_EventCategoryMap = {
    { EventCategory::TwoJets_ZeroBtag, EventCategory::TwoJets_ZeroLooseBtag },
    { EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_OneLooseBtag },
    { EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_TwoLooseBtag },
    { EventCategory::TwoJets_AtLeastOneBtag, EventCategory::TwoJets_AtLeastOneLooseBtag }
};

static const EventCategorySet TwoJetsEventCategories_MediumBjets =
    tools::collect_map_keys<decltype(MediumToLoose_EventCategoryMap), EventCategorySet>(MediumToLoose_EventCategoryMap);

static const EventCategorySet TwoJetsEventCategories_LooseBjets =
    tools::collect_map_values<decltype(MediumToLoose_EventCategoryMap), EventCategorySet>(MediumToLoose_EventCategoryMap);

using EventRegionSet = EnumNameMap<EventRegion>::EnumEntrySet;
using EventRegionMap = std::map<EventRegion, EventRegion>;

static const EventRegionMap HighMt_LowMt_RegionMap =
        { { EventRegion::OS_Iso_HighMt, EventRegion::OS_Isolated },
          { EventRegion::SS_Iso_HighMt, EventRegion::SS_Isolated },
          { EventRegion::OS_AntiIso_HighMt, EventRegion::OS_AntiIsolated },
          { EventRegion::SS_AntiIso_HighMt, EventRegion::SS_AntiIsolated } };

static const EventRegionSet HighMtRegions = { EventRegion::OS_Iso_HighMt, EventRegion::SS_Iso_HighMt,
                                              EventRegion::OS_AntiIso_HighMt, EventRegion::SS_AntiIso_HighMt };

static const EventRegionSet QcdRegions = { EventRegion::OS_Isolated, EventRegion::SS_Isolated,
                                           EventRegion::OS_AntiIsolated, EventRegion::SS_AntiIsolated};

static const EventRegionSet AllEventRegions = __EventRegion_names<>::names.GetEnumEntries();

enum class HTbinning { lt0 = -1, lt100 = 0, f100to200 = 1, f200to400 = 2, f400to600 = 3, gt600 = 4 };
inline HTbinning GetHTbin(double HT)
{
    if(HT < 0) return HTbinning::lt0;
    if(HT < 100) return HTbinning::lt100;
    if(HT < 200) return HTbinning::f100to200;
    if(HT < 400) return HTbinning::f200to400;
    if(HT < 600) return HTbinning::f400to600;
    return HTbinning::gt600;
}

//EventCategoryVector DetermineEventCategories(const std::vector<float>& csv_Bjets,
//                                             const EventInfoBase::BjetPair& selected_bjets, Int_t nBjets_retagged,
//                                             double CSVL, double CSVM, bool doRetag = false)
//{
//    EventCategoryVector categories;
//    categories.push_back(EventCategory::Inclusive);

//    static const std::map< size_t, EventCategory> mediumCategories_map {
//        {  0 , EventCategory::TwoJets_ZeroBtag }, { 1 , EventCategory::TwoJets_OneBtag },
//        {  2 , EventCategory::TwoJets_TwoBtag }
//    };

//    static const std::map< size_t, EventCategory> looseCategories_map {
//        { 0 , EventCategory::TwoJets_ZeroLooseBtag }, { 1, EventCategory::TwoJets_OneLooseBtag },
//        { 2 , EventCategory::TwoJets_TwoLooseBtag }
//    };

//    if (selected_bjets.first < csv_Bjets.size() && selected_bjets.second < csv_Bjets.size()){
//        categories.push_back(EventCategory::TwoJets_Inclusive);

//        size_t n_mediumBtag = 0;
//        if(doRetag) {
//            n_mediumBtag = std::min<size_t>(nBjets_retagged, 2);
//        } else {
//            if(csv_Bjets.at(selected_bjets.first) > CSVM) ++n_mediumBtag;
//            if(csv_Bjets.at(selected_bjets.second) > CSVM) ++n_mediumBtag;
//        }

//        if(mediumCategories_map.count(n_mediumBtag))
//            categories.push_back(mediumCategories_map.at(n_mediumBtag));
//        if(n_mediumBtag > 0)
//            categories.push_back(EventCategory::TwoJets_AtLeastOneBtag);

//        size_t n_looseBtag = 0;
//        if(csv_Bjets.at(selected_bjets.first) > CSVL) ++n_looseBtag;
//        if(csv_Bjets.at(selected_bjets.second) > CSVL) ++n_looseBtag;

//        if(looseCategories_map.count(n_looseBtag))
//            categories.push_back(looseCategories_map.at(n_looseBtag));
//        if(n_looseBtag > 0)
//            categories.push_back(EventCategory::TwoJets_AtLeastOneLooseBtag);
//    }

//    return categories;
//}

} // namespace analysis
