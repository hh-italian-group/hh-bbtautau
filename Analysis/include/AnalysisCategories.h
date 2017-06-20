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

//struct DataCategory {
//    using SFMap = std::map<std::string, double>;

//    std::string name;
//    std::string title;
//    std::string datacard;
//    EColor color;
//    double limits_sf;
//    bool draw;
//    unsigned draw_sf;
//    bool isCategoryToSubtract;

//    std::set<DataCategoryType> types;
//    std::set<Channel> channels;
//    std::set<std::string> sub_categories;
//    SFMap sources_sf;
//    std::map<unsigned, double> exclusive_sf;
//    std::set<std::string> uncertainties;

//    DataCategory()
//        : color(kBlack), limits_sf(1.0), draw(false), draw_sf(1), isCategoryToSubtract(true) {}

//    bool IsSignal() const { return types.count(DataCategoryType::Signal); }
//    bool IsBackground() const { return types.count(DataCategoryType::Background); }
//    bool IsData() const { return types.count(DataCategoryType::Data); }
//    bool IsComposit() const { return types.count(DataCategoryType::Composit); }
//};

//using DataCategoryMap = std::map<std::string, DataCategory>;
//using DataCategoryPtrSet = std::set<const DataCategory*>;
//using DataCategoryPtrMap = std::map<std::string, const DataCategory*>;
//using DataCategoryTypeSet = std::set<DataCategoryType>;
//using DataCategoryPtrVector = std::vector<const DataCategory*>;
//using DataCategoryTypeMap = std::map<DataCategoryType, DataCategoryPtrSet>;

//static const DataCategoryTypeSet dataCategoryTypeForQCD  = { DataCategoryType::QCD,
//                                                             DataCategoryType::QCD_alternative };

//class DataCategoryCollection {
//public:
//    DataCategoryCollection(const std::string& sources_cfg_name, const std::string& signal_list, Channel channel_id)
//    {
//        DataCategory category;
//        std::ifstream cfg(sources_cfg_name);
//        size_t line_number = 0;
//        while(ReadNextCategory(cfg, line_number, category)) {
//            if(category.channels.size() && !category.channels.count(channel_id)) continue;
//            CheckCategoryValidity(category);
//            categories[category.name] = category;
//            all_categories.push_back(&categories[category.name]);
//            for(DataCategoryType type : category.types)
//                categories_by_type[type].insert(&categories[category.name]);
//            for(const auto& source_entry : category.sources_sf)
//                all_sources.insert(source_entry.first);
//            if(category.datacard.size())
//                categories_by_datacard[category.datacard] = &categories[category.name];
//        }
//        const auto& signal_names = ParseSignalList(signal_list);
//        for(const auto& signal_name : signal_names) {
//            if(!signal_name.size()) continue;
//            if(!categories.count(signal_name))
//                throw exception("Undefined signal '%1%'.") % signal_name;
//            categories[signal_name].draw = true;
//        }
//    }

//    DataCategoryCollection(const DataCategoryCollection& other)
//        : all_sources(other.all_sources), categories(other.categories)
//    {
//        for(const auto& category_entry : categories) {
//            for(DataCategoryType type : category_entry.second.types)
//                categories_by_type[type].insert(&category_entry.second);
//        }

//        for(const DataCategory* other_category : other.all_categories)
//            all_categories.push_back(&categories.at(other_category->name));
//    }

//    const DataCategoryPtrVector& GetAllCategories() const { return all_categories; }
//    const DataCategoryPtrSet& GetCategories(DataCategoryType dataCategoryType) const
//    {
//        if(!categories_by_type.count(dataCategoryType))
//            return empty_category_set;
//        return categories_by_type.at(dataCategoryType);
//    }

//    const DataCategory& GetUniqueCategory(DataCategoryType dataCategoryType) const
//    {
//        if(!categories_by_type.at(dataCategoryType).size())
//            throw exception("Unique category for data category type '%1%' not found.") % dataCategoryType;
//        if(categories_by_type.at(dataCategoryType).size() != 1)
//            throw exception("More than one category for data category type '%1%'.") % dataCategoryType;
//        return *(*categories_by_type.at(dataCategoryType).begin());
//    }

//    const DataCategory& FindCategory(const std::string& name) const
//    {
//        if(!categories.count(name))
//            throw exception("Data category '%1%' not found.") % name;
//        return categories.at(name);
//    }

//    const DataCategory& FindCategoryForDatacard(const std::string& datacard) const
//    {
//        if(!categories_by_datacard.count(datacard))
//            throw exception("Data category for datacard '%1%' not found.") % datacard;
//        return *categories_by_datacard.at(datacard);
//    }

//private:
//    void CheckCategoryValidity(const DataCategory& category) const
//    {
//        if(categories.count(category.name))
//            throw exception("Category with name '%1%' is already defined.") % category.name;
//        if(category.sub_categories.size() && !category.types.count(DataCategoryType::Composit))
//            throw exception("Not composit category '%1%' may not contain sub-categories.") % category.name;
//        if(category.types.count(DataCategoryType::Composit) && category.sources_sf.size())
//            throw exception("Composit category '%1%' may not contain direct file definitions.") % category.name;
//        for(const auto& sub_category : category.sub_categories) {
//            if(!categories.count(sub_category))
//                throw exception("Sub-category '%1%' for category '%2%' is not defined.") % sub_category % category.name;
//            if(categories.at(sub_category).types.count(DataCategoryType::Composit))
//                throw exception("Invalid sub-category '%1%' for category '%2%'. Composit category hierarchy is not"
//                                " supported.") % sub_category % category.name;
//        }
////        for(const auto& source_entry : category.sources_sf) {
////            if(all_sources.count(source_entry.first))
////                throw exception("Source '") << source_entry.first << "' is already part of the other data category.";
////        }
//        if(category.datacard.size() && categories_by_datacard.count(category.datacard))
//            throw exception("Category for datacard '%1%' is already defined.") % category.datacard;
//    }

//    static bool ReadNextCategory(std::istream& cfg, size_t& line_number, DataCategory& category)
//    {
//        category = DataCategory();
//        bool category_started = false;
//        while (cfg.good()) {
//            std::string cfgLine;
//            std::getline(cfg,cfgLine);
//            ++line_number;
//            if ((cfgLine.size() && cfgLine.at(0) == '#') || (!cfgLine.size() && !category_started)) continue;
//            if(!cfgLine.size())
//                return true;
//            if(!category_started && cfgLine.at(0) == '[') {
//                const size_t pos = cfgLine.find(']');
//                if(pos == std::string::npos)
//                    throw exception("bad source config syntax in line %1%.") % line_number;
//                category.name = cfgLine.substr(1, pos - 1);
//                category_started = true;
//            } else if(category_started) {
//                ReadParameterLine(cfgLine, line_number, category);
//            } else
//                throw exception("bad source config syntax in line %1%.") % line_number;
//        }
//        return category_started;
//    }

//    static void ReadParameterLine(const std::string& cfgLine, size_t line_number, DataCategory& category)
//    {
//        static const char separator = ':';

//        const size_t pos = cfgLine.find(separator);
//        if(pos == std::string::npos)
//            throw exception("bad source config syntax for a parameter in line %1%.") % line_number;
//        const std::string param_name = cfgLine.substr(0, pos);
//        if(pos + 2 >= cfgLine.size())
//            throw exception("empty parameter value in source config in line %1%.") % line_number;
//        const std::string param_value = cfgLine.substr(pos + 2);
//        std::istringstream ss(param_value);
//        ss >> std::boolalpha;
//        if(param_name == "type") {
//            DataCategoryType type;
//            ss >> type;
//            category.types.insert(type);
//        } else if(param_name == "title") {
//            category.title = param_value;
//        } else if(param_name == "color") {
//            ss >> category.color;
//        } else if(param_name == "file") {
//            std::string file_name;
//            double scale_factor;
//            ss >> file_name;
//            ss >> scale_factor;
//            category.sources_sf[file_name] = scale_factor;
//        } else if(param_name == "limits_sf") {
//            ss >> category.limits_sf;
//        } else if(param_name == "draw_sf") {
//            ss >> category.draw_sf;
//        } else if(param_name == "draw") {
//            ss >> category.draw;
//        } else if(param_name == "isCategoryToSubtract") {
//            ss >> category.isCategoryToSubtract;
//        } else if(param_name == "channel") {
//            Channel channel_id;
//            ss >> channel_id;
//            category.channels.insert(channel_id);
//        } else if(param_name == "datacard") {
//            ss >> category.datacard;
//        } else if(param_name == "subcategory") {
//            category.sub_categories.insert(param_value);
//        } else if(param_name == "uncertainty") {
//            category.uncertainties.insert(param_value);
//        } else if(param_name == "exclusive_sf") {
//            unsigned n_jets;
//            double scale_factor;
//            ss >> n_jets;
//            ss >> scale_factor;
//            category.exclusive_sf[n_jets] = scale_factor;
//        } else
//            throw exception("Unsupported parameter '%1%' in configuration line %2%.") % param_name % line_number;
//    }

//    static std::set<std::string> ParseSignalList(const std::string& signal_list)
//    {
//        static const char separator = ',';

//        std::set<std::string> result;
//        size_t prev_pos = 0;
//        for(bool next = true; next;) {
//            const size_t pos = signal_list.find(separator, prev_pos);
//            next = pos != std::string::npos;
//            const size_t last_pos = next ? pos - 1 : std::string::npos;
//            const std::string signal_name = signal_list.substr(prev_pos, last_pos);
//            result.insert(signal_name);
//            prev_pos = pos + 1;
//        }
//        return result;
//    }

//private:
//    std::set<std::string> all_sources;
//    DataCategoryMap categories;
//    DataCategoryPtrVector all_categories;
//    DataCategoryTypeMap categories_by_type;
//    DataCategoryPtrSet empty_category_set;
//    DataCategoryPtrMap categories_by_datacard;
//};

//std::ostream& operator<<(std::ostream& s, const DataCategory& category){
//    s << "Name: " << category.name << ", Title: '" << category.title << "', Color: " << category.color << std::endl;
//    for(const auto& source : category.sources_sf)
//        s << "File: " << source.first << ", SF: " << source.second << "\n";
//    return s;
//}

//std::ostream& operator<<(std::ostream& s, const DataCategoryPtrSet& dataCategories)
//{
//    bool first_category = true;
//    for(const auto& dataCategory : dataCategories) {
//        if(!first_category)
//            s << ", ";
//        s << dataCategory->name;
//        first_category = false;
//    }
//    return s;
//}

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
