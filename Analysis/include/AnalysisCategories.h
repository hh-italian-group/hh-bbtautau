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
#include "AnalysisTools/Print/include/RootPrintTools.h"
#include "h-tautau/Analysis/include/FlatEventInfo.h"

namespace analysis {

typedef std::map<std::string, double> DataSourceScaleFactorMap;

enum class DataCategoryType { Signal, Background, Data, DYJets, DYJets_incl, DYJets_excl, ZL, ZJ, ZL_MC, ZJ_MC, ZTT,
                              ZTT_MC, ZTT_L, Embedded, TT_Embedded, Limits, Composit, QCD, QCD_alternative, WJets,
                              WJets_MC, DiBoson_MC, DiBoson};
static const std::map<DataCategoryType, std::string> dataCategoryTypeNameMap = {
    { DataCategoryType::Signal, "SIGNAL" }, { DataCategoryType::Background, "BACKGROUND" },
    { DataCategoryType::Data, "DATA" }, { DataCategoryType::DYJets, "DY_JETS" },
    { DataCategoryType::DYJets_incl, "DY_JETS_incl" },{ DataCategoryType::DYJets_excl, "DY_JETS_excl" },
    { DataCategoryType::ZL, "ZL" }, { DataCategoryType::ZJ, "ZJ" },
    { DataCategoryType::ZL_MC, "ZL_MC" }, { DataCategoryType::ZJ_MC, "ZJ_MC" }, { DataCategoryType::ZTT, "ZTT" },
    { DataCategoryType::ZTT_MC, "ZTT_MC" }, { DataCategoryType::ZTT_L, "ZTT_L" }, { DataCategoryType::Embedded, "EMBEDDED" },
    { DataCategoryType::TT_Embedded, "TT_EMBEDDED" },
    { DataCategoryType::Limits, "LIMITS" }, { DataCategoryType::Composit, "COMPOSIT" },
    { DataCategoryType::QCD, "QCD" }, { DataCategoryType::QCD_alternative, "QCD_alternative" },
    { DataCategoryType::WJets, "W_JETS" }, { DataCategoryType::WJets_MC, "W_JETS_MC" },
    { DataCategoryType::DiBoson, "DiBoson" }, { DataCategoryType::DiBoson_MC, "DiBoson_MC" }
};

std::ostream& operator<< (std::ostream& s, const DataCategoryType& dataCategoryType) {
    s << dataCategoryTypeNameMap.at(dataCategoryType);
    return s;
}
std::istream& operator>> (std::istream& s, DataCategoryType& dataCategoryType) {
    std::string name;
    s >> name;
    for(const auto& map_entry : dataCategoryTypeNameMap) {
        if(map_entry.second == name) {
            dataCategoryType = map_entry.first;
            return s;
        }
    }
    throw exception("Unknown data category type '") << name << "'.";
}

struct DataCategory {
    std::string name;
    std::string title;
    std::string datacard;
    EColor color;
    double limits_sf;
    bool draw;
    unsigned draw_sf;
    bool isCategoryToSubtract;

    std::set<DataCategoryType> types;
    std::set<Channel> channels;
    std::set<std::string> sub_categories;
    DataSourceScaleFactorMap sources_sf;
    std::map<unsigned, double> exclusive_sf;
    std::set<std::string> uncertainties;

    DataCategory()
        : color(kBlack), limits_sf(1.0), draw(false), draw_sf(1), isCategoryToSubtract(true) {}

    bool IsSignal() const { return types.count(DataCategoryType::Signal); }
    bool IsBackground() const { return types.count(DataCategoryType::Background); }
    bool IsData() const { return types.count(DataCategoryType::Data); }
    bool IsComposit() const { return types.count(DataCategoryType::Composit); }
};

typedef std::map<std::string, DataCategory> DataCategoryMap;
typedef std::set<const DataCategory*> DataCategoryPtrSet;
typedef std::map<std::string, const DataCategory*> DataCategoryPtrMap;
typedef std::set<DataCategoryType> DataCategoryTypeSet;
typedef std::vector<const DataCategory*> DataCategoryPtrVector;
typedef std::map<DataCategoryType, DataCategoryPtrSet> DataCategoryTypeMap;

static const DataCategoryTypeSet dataCategoryTypeForQCD  = { DataCategoryType::QCD,
                                                             DataCategoryType::QCD_alternative };

class DataCategoryCollection {
public:
    DataCategoryCollection(const std::string& sources_cfg_name, const std::string& signal_list, Channel channel_id)
    {
        DataCategory category;
        std::ifstream cfg(sources_cfg_name);
        size_t line_number = 0;
        while(ReadNextCategory(cfg, line_number, category)) {
            if(category.channels.size() && !category.channels.count(channel_id)) continue;
            CheckCategoryValidity(category);
            categories[category.name] = category;
            all_categories.push_back(&categories[category.name]);
            for(DataCategoryType type : category.types)
                categories_by_type[type].insert(&categories[category.name]);
            for(const auto& source_entry : category.sources_sf)
                all_sources.insert(source_entry.first);
            if(category.datacard.size())
                categories_by_datacard[category.datacard] = &categories[category.name];
        }
        const auto& signal_names = ParseSignalList(signal_list);
        for(const auto& signal_name : signal_names) {
            if(!signal_name.size()) continue;
            if(!categories.count(signal_name))
                throw exception("Undefined signal '") << signal_name << "'.";
            categories[signal_name].draw = true;
        }
    }

    DataCategoryCollection(const DataCategoryCollection& other)
        : all_sources(other.all_sources), categories(other.categories)
    {
        for(const auto& category_entry : categories) {
            for(DataCategoryType type : category_entry.second.types)
                categories_by_type[type].insert(&category_entry.second);
        }

        for(const DataCategory* other_category : other.all_categories)
            all_categories.push_back(&categories.at(other_category->name));
    }

    const DataCategoryPtrVector& GetAllCategories() const { return all_categories; }
    const DataCategoryPtrSet& GetCategories(DataCategoryType dataCategoryType) const
    {
        if(!categories_by_type.count(dataCategoryType))
            return empty_category_set;
        return categories_by_type.at(dataCategoryType);
    }

    const DataCategory& GetUniqueCategory(DataCategoryType dataCategoryType) const
    {
        if(!categories_by_type.at(dataCategoryType).size())
            throw exception("Unique category for data category type '") << dataCategoryType << "' not found.";
        if(categories_by_type.at(dataCategoryType).size() != 1)
            throw exception("More than one category for data category type '") << dataCategoryType << "'.";
        return *(*categories_by_type.at(dataCategoryType).begin());
    }

    const DataCategory& FindCategory(const std::string& name) const
    {
        if(!categories.count(name))
            throw exception("Data category '") << name << "' not found.";
        return categories.at(name);
    }

    const DataCategory& FindCategoryForDatacard(const std::string& datacard) const
    {
        if(!categories_by_datacard.count(datacard))
            throw exception("Data category for datacard '") << datacard << "' not found.";
        return *categories_by_datacard.at(datacard);
    }

private:
    void CheckCategoryValidity(const DataCategory& category) const
    {
        if(categories.count(category.name))
            throw exception("Category with name '") << category.name << "' is already defined.";
        if(category.sub_categories.size() && !category.types.count(DataCategoryType::Composit))
            throw exception("Not composit category '") << category.name << "' may not contain sub-categories.";
        if(category.types.count(DataCategoryType::Composit) && category.sources_sf.size())
            throw exception("Composit category '") << category.name << "' may not contain direct file definitions.";
        for(const auto& sub_category : category.sub_categories) {
            if(!categories.count(sub_category))
                throw exception("Sub-category '") << sub_category << "' for category '"
                                                  << category.name << "' is not defined.";
            if(categories.at(sub_category).types.count(DataCategoryType::Composit))
                throw exception("Invalid sub-category '") << sub_category << "' for category '" << category.name
                                                      << "'. Composit category hierarchy is not supported.";
        }
//        for(const auto& source_entry : category.sources_sf) {
//            if(all_sources.count(source_entry.first))
//                throw exception("Source '") << source_entry.first << "' is already part of the other data category.";
//        }
        if(category.datacard.size() && categories_by_datacard.count(category.datacard))
            throw exception("Category for datacard '") << category.datacard << "' is already defined.";
    }

    static bool ReadNextCategory(std::istream& cfg, size_t& line_number, DataCategory& category)
    {
        category = DataCategory();
        bool category_started = false;
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            ++line_number;
            if ((cfgLine.size() && cfgLine.at(0) == '#') || (!cfgLine.size() && !category_started)) continue;
            if(!cfgLine.size())
                return true;
            if(!category_started && cfgLine.at(0) == '[') {
                const size_t pos = cfgLine.find(']');
                if(pos == std::string::npos)
                    throw exception("bad source config syntax in line ") << line_number;
                category.name = cfgLine.substr(1, pos - 1);
                category_started = true;
            } else if(category_started) {
                ReadParameterLine(cfgLine, line_number, category);
            } else
                throw exception("bad source config syntax in line ") << line_number;
        }
        return category_started;
    }

    static void ReadParameterLine(const std::string& cfgLine, size_t line_number, DataCategory& category)
    {
        static const char separator = ':';

        const size_t pos = cfgLine.find(separator);
        if(pos == std::string::npos)
            throw exception("bad source config syntax for a parameter in line ") << line_number;
        const std::string param_name = cfgLine.substr(0, pos);
        if(pos + 2 >= cfgLine.size())
            throw exception("empty parameter value in source config in line ") << line_number;
        const std::string param_value = cfgLine.substr(pos + 2);
        std::istringstream ss(param_value);
        ss >> std::boolalpha;
        if(param_name == "type") {
            DataCategoryType type;
            ss >> type;
            category.types.insert(type);
        } else if(param_name == "title") {
            category.title = param_value;
        } else if(param_name == "color") {
            ss >> category.color;
        } else if(param_name == "file") {
            std::string file_name;
            double scale_factor;
            ss >> file_name;
            ss >> scale_factor;
            category.sources_sf[file_name] = scale_factor;
        } else if(param_name == "limits_sf") {
            ss >> category.limits_sf;
        } else if(param_name == "draw_sf") {
            ss >> category.draw_sf;
        } else if(param_name == "draw") {
            ss >> category.draw;
        } else if(param_name == "isCategoryToSubtract") {
            ss >> category.isCategoryToSubtract;
        } else if(param_name == "channel") {
            Channel channel_id;
            ss >> channel_id;
            category.channels.insert(channel_id);
        } else if(param_name == "datacard") {
            ss >> category.datacard;
        } else if(param_name == "subcategory") {
            category.sub_categories.insert(param_value);
        } else if(param_name == "uncertainty") {
            category.uncertainties.insert(param_value);
        } else if(param_name == "exclusive_sf") {
            unsigned n_jets;
            double scale_factor;
            ss >> n_jets;
            ss >> scale_factor;
            category.exclusive_sf[n_jets] = scale_factor;
        } else
            throw exception("Unsupported parameter '") << param_name << "' in configuration line " << line_number;
    }

    static std::set<std::string> ParseSignalList(const std::string& signal_list)
    {
        static const char separator = ',';

        std::set<std::string> result;
        size_t prev_pos = 0;
        for(bool next = true; next;) {
            const size_t pos = signal_list.find(separator, prev_pos);
            next = pos != std::string::npos;
            const size_t last_pos = next ? pos - 1 : std::string::npos;
            const std::string signal_name = signal_list.substr(prev_pos, last_pos);
            result.insert(signal_name);
            prev_pos = pos + 1;
        }
        return result;
    }

private:
    std::set<std::string> all_sources;
    DataCategoryMap categories;
    DataCategoryPtrVector all_categories;
    DataCategoryTypeMap categories_by_type;
    DataCategoryPtrSet empty_category_set;
    DataCategoryPtrMap categories_by_datacard;
};

std::ostream& operator<<(std::ostream& s, const DataCategory& category){
    s << "Name: " << category.name << ", Title: '" << category.title << "', Color: " << category.color << std::endl;
    for(const auto& source : category.sources_sf)
        s << "File: " << source.first << ", SF: " << source.second << "\n";
    return s;
}

std::ostream& operator<<(std::ostream& s, const DataCategoryPtrSet& dataCategories)
{
    bool first_category = true;
    for(const auto& dataCategory : dataCategories) {
        if(!first_category)
            s << ", ";
        s << dataCategory->name;
        first_category = false;
    }
    return s;
}

enum class EventRegion { Unknown = 0, OS_Isolated = 1, OS_AntiIsolated = 2, SS_Isolated = 3, SS_AntiIsolated = 4,
                         OS_Iso_HighMt = 5, SS_Iso_HighMt = 6, OS_AntiIso_HighMt = 7, SS_AntiIso_HighMt = 8 };

enum class EventCategory { Inclusive = 0, TwoJets_Inclusive = 1, TwoJets_ZeroBtag = 2, TwoJets_OneBtag = 3,
                           TwoJets_TwoBtag = 4, TwoJets_ZeroLooseBtag = 5, TwoJets_OneLooseBtag = 6,
                           TwoJets_TwoLooseBtag = 7, TwoJets_AtLeastOneBtag = 8, TwoJets_AtLeastOneLooseBtag = 9 };

enum class EventSubCategory { NoCuts = 0, KinematicFitConverged = 1, MassWindow = 2,
                              KinematicFitConvergedWithMassWindow = 3, OutsideMassWindow = 4,
                              KinematicFitConvergedOutsideMassWindow = 5
                            };

namespace detail {
static const std::map<EventCategory, std::string> eventCategoryNamesMap =
          { { EventCategory::Inclusive, "Inclusive" }, { EventCategory::TwoJets_Inclusive, "2jets" },
            { EventCategory::TwoJets_ZeroBtag, "2jets0btag" }, { EventCategory::TwoJets_OneBtag, "2jets1btag"},
            { EventCategory::TwoJets_TwoBtag, "2jets2btag" },
            { EventCategory::TwoJets_ZeroLooseBtag, "2jets0Loosebtag" },
            { EventCategory::TwoJets_OneLooseBtag, "2jets1Loosebtag" },
            { EventCategory::TwoJets_TwoLooseBtag, "2jets2Loosebtag" },
            { EventCategory::TwoJets_AtLeastOneBtag, "2jets_at_least_1btag" },
          { EventCategory::TwoJets_AtLeastOneLooseBtag, "2jets_at_least_1Loosebtag" }};

static const std::map<EventRegion, std::string> eventRegionNamesMap =
          { { EventRegion::Unknown, "Unknown"}, { EventRegion::OS_Isolated, "OS_Isolated"},
            { EventRegion::OS_AntiIsolated, "OS_AntiIsolated"}, { EventRegion::SS_Isolated, "SS_Isolated"},
            { EventRegion::SS_AntiIsolated, "SS_AntiIsolated"}, { EventRegion::OS_Iso_HighMt, "OS_Iso_HighMt"},
            { EventRegion::SS_Iso_HighMt, "SS_Iso_HighMt"} , { EventRegion::OS_AntiIso_HighMt, "OS_AntiIso_HighMt"},
            { EventRegion::SS_AntiIso_HighMt, "SS_AntiIso_HighMt"} };

static const std::map<EventSubCategory, std::string> eventSubCategoryNamesMap =
          { { EventSubCategory::NoCuts, "NoCuts" }, { EventSubCategory::KinematicFitConverged, "KinFitConverged" },
            { EventSubCategory::MassWindow, "MassWindow" },
            { EventSubCategory::KinematicFitConvergedWithMassWindow, "KinFitConvergedWithMassWindow" },
            { EventSubCategory::KinematicFitConvergedOutsideMassWindow, "KinematicFitConvergedOutsideMassWindow" },
            { EventSubCategory::OutsideMassWindow, "OutsideMassWindow" } };
} // namespace detail

typedef std::vector<EventCategory> EventCategoryVector;
typedef std::set<EventCategory> EventCategorySet;
typedef std::map<EventCategory, EventCategory> EventCategoryMap;
typedef std::set<EventSubCategory> EventSubCategorySet;

static const EventCategorySet AllEventCategories = tools::collect_map_keys(detail::eventCategoryNamesMap);
static const EventSubCategorySet AllEventSubCategories = tools::collect_map_keys(detail::eventSubCategoryNamesMap);

static const EventCategoryMap MediumToLoose_EventCategoryMap =
        { { EventCategory::TwoJets_ZeroBtag, EventCategory::TwoJets_ZeroLooseBtag },
          { EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_OneLooseBtag },
          { EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_TwoLooseBtag },
          { EventCategory::TwoJets_AtLeastOneBtag, EventCategory::TwoJets_AtLeastOneLooseBtag }};

static const EventCategorySet TwoJetsEventCategories_MediumBjets =
                                                                tools::collect_map_keys(MediumToLoose_EventCategoryMap);

static const EventCategorySet TwoJetsEventCategories_LooseBjets =
                                                              tools::collect_map_values(MediumToLoose_EventCategoryMap);

typedef std::set<EventRegion> EventRegionSet;
typedef std::map<EventRegion, EventRegion> EventRegionMap;

static const EventRegionMap HighMt_LowMt_RegionMap =
        { { EventRegion::OS_Iso_HighMt, EventRegion::OS_Isolated },
          { EventRegion::SS_Iso_HighMt, EventRegion::SS_Isolated },
          { EventRegion::OS_AntiIso_HighMt, EventRegion::OS_AntiIsolated },
          { EventRegion::SS_AntiIso_HighMt, EventRegion::SS_AntiIsolated } };

static const EventRegionSet HighMtRegions = { EventRegion::OS_Iso_HighMt, EventRegion::SS_Iso_HighMt,
                                              EventRegion::OS_AntiIso_HighMt, EventRegion::SS_AntiIso_HighMt };

static const EventRegionSet QcdRegions = { EventRegion::OS_Isolated, EventRegion::SS_Isolated,
                                           EventRegion::OS_AntiIsolated, EventRegion::SS_AntiIsolated};

static const EventRegionSet AllEventRegions = tools::collect_map_keys(detail::eventRegionNamesMap);

std::ostream& operator<<(std::ostream& s, const EventCategory& eventCategory) {
    s << detail::eventCategoryNamesMap.at(eventCategory);
    return s;
}

std::wostream& operator<<(std::wostream& s, const EventCategory& eventCategory) {
    const std::string str = detail::eventCategoryNamesMap.at(eventCategory);
    s << std::wstring(str.begin(), str.end());
    return s;
}

std::istream& operator>> (std::istream& s, EventCategory& eventCategory)
{
    std::string name;
    s >> name;
    for(const auto& map_entry : detail::eventCategoryNamesMap) {
        if(map_entry.second == name) {
            eventCategory = map_entry.first;
            return s;
        }
    }
    throw exception("Unknown event category '") << name << "'.";
}

std::ostream& operator<<(std::ostream& s, const EventRegion& eventRegion) {
    s << detail::eventRegionNamesMap.at(eventRegion);
    return s;
}

std::istream& operator>> (std::istream& s, EventRegion& eventRegion)
{
    std::string name;
    s >> name;
    for(const auto& map_entry : detail::eventRegionNamesMap) {
        if(map_entry.second == name) {
            eventRegion = map_entry.first;
            return s;
        }
    }
    throw exception("Unknown event region '") << name << "'.";
}

std::ostream& operator<<(std::ostream& s, const EventSubCategory& eventSubCategory) {
    s << detail::eventSubCategoryNamesMap.at(eventSubCategory);
    return s;
}

std::istream& operator>> (std::istream& s, EventSubCategory& eventSubCategory)
{
    std::string name;
    s >> name;
    for(const auto& map_entry : detail::eventSubCategoryNamesMap) {
        if(map_entry.second == name) {
            eventSubCategory = map_entry.first;
            return s;
        }
    }
    throw exception("Unknown event sub-category '") << name << "'.";
}

EventCategoryVector DetermineEventCategories(const std::vector<float>& csv_Bjets,
                                             const FlatEventInfo::BjetPair& selected_bjets, Int_t nBjets_retagged,
                                             double CSVL, double CSVM, bool doRetag = false)
{
    EventCategoryVector categories;
    categories.push_back(EventCategory::Inclusive);

    static const std::map< size_t, EventCategory> mediumCategories_map {
        {  0 , EventCategory::TwoJets_ZeroBtag }, { 1 , EventCategory::TwoJets_OneBtag },
        {  2 , EventCategory::TwoJets_TwoBtag }
    };

    static const std::map< size_t, EventCategory> looseCategories_map {
        { 0 , EventCategory::TwoJets_ZeroLooseBtag }, { 1, EventCategory::TwoJets_OneLooseBtag },
        { 2 , EventCategory::TwoJets_TwoLooseBtag }
    };

    if (selected_bjets.first < csv_Bjets.size() && selected_bjets.second < csv_Bjets.size()){
        categories.push_back(EventCategory::TwoJets_Inclusive);

        size_t n_mediumBtag = 0;
        if(doRetag) {
            n_mediumBtag = std::min<size_t>(nBjets_retagged, 2);
        } else {
            if(csv_Bjets.at(selected_bjets.first) > CSVM) ++n_mediumBtag;
            if(csv_Bjets.at(selected_bjets.second) > CSVM) ++n_mediumBtag;
        }

        if(mediumCategories_map.count(n_mediumBtag))
            categories.push_back(mediumCategories_map.at(n_mediumBtag));
        if(n_mediumBtag > 0)
            categories.push_back(EventCategory::TwoJets_AtLeastOneBtag);

        size_t n_looseBtag = 0;
        if(csv_Bjets.at(selected_bjets.first) > CSVL) ++n_looseBtag;
        if(csv_Bjets.at(selected_bjets.second) > CSVL) ++n_looseBtag;

        if(looseCategories_map.count(n_looseBtag))
            categories.push_back(looseCategories_map.at(n_looseBtag));
        if(n_looseBtag > 0)
            categories.push_back(EventCategory::TwoJets_AtLeastOneLooseBtag);
    }

    return categories;
}

} // namespace analysis
