/*! Collection of histogram containers for flat tree analyzers.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "FlatAnalyzerData.h"
#include "h-tautau/Analysis/include/Htautau_2015.h"
#include "custom_cuts.h"

namespace analysis {

struct FlatAnalyzerDataId;
typedef std::set<FlatAnalyzerDataId> FlatAnalyzerDataIdSet;

struct FlatAnalyzerDataId {
    EventCategory eventCategory;
    EventSubCategory eventSubCategory;
    EventRegion eventRegion;
    EventEnergyScale eventEnergyScale;
    std::string dataCategoryName;

    FlatAnalyzerDataId()
        : eventCategory(EventCategory::Inclusive), eventSubCategory(EventSubCategory::NoCuts),
          eventRegion(EventRegion::OS_Isolated), eventEnergyScale(EventEnergyScale::Central)
    {}

    FlatAnalyzerDataId(EventCategory _eventCategory, EventSubCategory _eventSubCategory, EventRegion _eventRegion,
                       EventEnergyScale _eventEnergyScale, const std::string& _dataCategoryName)
        : eventCategory(_eventCategory), eventSubCategory(_eventSubCategory), eventRegion(_eventRegion),
          eventEnergyScale(_eventEnergyScale), dataCategoryName(_dataCategoryName)
    {}

    bool operator< (const FlatAnalyzerDataId& other) const
    {
        if(eventCategory < other.eventCategory) return true;
        if(eventCategory > other.eventCategory) return false;
        if(eventSubCategory < other.eventSubCategory) return true;
        if(eventSubCategory > other.eventSubCategory) return false;
        if(eventRegion < other.eventRegion) return true;
        if(eventRegion > other.eventRegion) return false;
        if(eventEnergyScale < other.eventEnergyScale) return true;
        if(eventEnergyScale > other.eventEnergyScale) return false;
        return dataCategoryName < other.dataCategoryName;
    }

    std::string GetName() const
    {
        static const std::string separator = "/";
        std::ostringstream ss;
        ss << eventCategory << separator << eventSubCategory << separator << eventRegion << separator
           << eventEnergyScale << separator << dataCategoryName;
        return ss.str();
    }
};

std::ostream& operator<< (std::ostream& s, const FlatAnalyzerDataId& id)
{
    s << id.GetName();
    return s;
}

template<typename ...Types>
struct FlatAnalyzerDataMetaId;

template<>
struct FlatAnalyzerDataMetaId<EventCategory, EventRegion, EventEnergyScale, std::string> {
    EventCategory eventCategory;
    EventRegion eventRegion;
    EventEnergyScale eventEnergyScale;
    std::string dataCategoryName;

    FlatAnalyzerDataMetaId()
        : eventCategory(EventCategory::Inclusive), eventRegion(EventRegion::OS_Isolated),
          eventEnergyScale(EventEnergyScale::Central) {}

    FlatAnalyzerDataMetaId(EventCategory _eventCategory, EventRegion _eventRegion, EventEnergyScale _eventEnergyScale,
                           const std::string& _dataCategoryName)
        : eventCategory(_eventCategory), eventRegion(_eventRegion),
          eventEnergyScale(_eventEnergyScale), dataCategoryName(_dataCategoryName)  {}

    FlatAnalyzerDataId MakeId(EventSubCategory eventSubCategory) const
    {
        return FlatAnalyzerDataId(eventCategory, eventSubCategory, eventRegion, eventEnergyScale, dataCategoryName);
    }

    bool operator<(const FlatAnalyzerDataMetaId<EventCategory, EventRegion, EventEnergyScale, std::string>& other) const
    {
        if(eventCategory < other.eventCategory) return true;
        if(eventCategory > other.eventCategory) return false;
        if(eventRegion < other.eventRegion) return true;
        if(eventRegion > other.eventRegion) return false;
        if(eventEnergyScale < other.eventEnergyScale) return true;
        if(eventEnergyScale > other.eventEnergyScale) return false;
        return dataCategoryName < other.dataCategoryName;
    }

    std::string GetName() const
    {
        static const std::string separator = "/";
        static const std::string wildcard = "*";
        std::ostringstream ss;
        ss << eventCategory << separator << wildcard << separator << eventRegion << separator
           << eventEnergyScale << separator << dataCategoryName;
        return ss.str();
    }
};

typedef FlatAnalyzerDataMetaId<EventCategory, EventRegion, EventEnergyScale, std::string> FlatAnalyzerDataMetaId_noSub;

template<>
struct FlatAnalyzerDataMetaId<EventCategory, EventRegion, std::string> {
    EventCategory eventCategory;
    EventRegion eventRegion;
    std::string dataCategoryName;

    FlatAnalyzerDataMetaId()
        : eventCategory(EventCategory::Inclusive), eventRegion(EventRegion::OS_Isolated) {}

    FlatAnalyzerDataMetaId(EventCategory _eventCategory, EventRegion _eventRegion, const std::string& _dataCategoryName)
        : eventCategory(_eventCategory), eventRegion(_eventRegion), dataCategoryName(_dataCategoryName)  {}

    FlatAnalyzerDataId MakeId(EventSubCategory eventSubCategory, EventEnergyScale eventEnergyScale) const
    {
        return FlatAnalyzerDataId(eventCategory, eventSubCategory, eventRegion, eventEnergyScale, dataCategoryName);
    }

    FlatAnalyzerDataMetaId_noSub MakeMetaId(EventEnergyScale eventEnergyScale) const
    {
        return FlatAnalyzerDataMetaId_noSub(eventCategory, eventRegion, eventEnergyScale, dataCategoryName);
    }

    bool operator< (const FlatAnalyzerDataMetaId<EventCategory, EventRegion, std::string>& other) const
    {
        if(eventCategory < other.eventCategory) return true;
        if(eventCategory > other.eventCategory) return false;
        if(eventRegion < other.eventRegion) return true;
        if(eventRegion > other.eventRegion) return false;
        return dataCategoryName < other.dataCategoryName;
    }

    std::string GetName() const
    {
        static const std::string separator = "/";
        static const std::string wildcard = "*";
        std::ostringstream ss;
        ss << eventCategory << separator << wildcard << separator << eventRegion << separator
           << wildcard << separator << dataCategoryName;
        return ss.str();
    }
};

typedef FlatAnalyzerDataMetaId<EventCategory, EventRegion, std::string> FlatAnalyzerDataMetaId_noSub_noES;

template<>
struct FlatAnalyzerDataMetaId<EventCategory, EventSubCategory, EventRegion, EventEnergyScale> {
    EventCategory eventCategory;
    EventSubCategory eventSubCategory;
    EventRegion eventRegion;
    EventEnergyScale eventEnergyScale;

    FlatAnalyzerDataMetaId()
        : eventCategory(EventCategory::Inclusive), eventSubCategory(EventSubCategory::NoCuts),
          eventRegion(EventRegion::OS_Isolated), eventEnergyScale(EventEnergyScale::Central) {}

    FlatAnalyzerDataMetaId(EventCategory _eventCategory, EventSubCategory _eventSubCategory, EventRegion _eventRegion,
                           EventEnergyScale _eventEnergyScale)
        : eventCategory(_eventCategory), eventSubCategory(_eventSubCategory), eventRegion(_eventRegion),
          eventEnergyScale(_eventEnergyScale) {}

    FlatAnalyzerDataId MakeId(const std::string& dataCategoryName) const
    {
        return FlatAnalyzerDataId(eventCategory, eventSubCategory, eventRegion, eventEnergyScale, dataCategoryName);
    }

    bool operator< (const FlatAnalyzerDataMetaId<EventCategory, EventSubCategory, EventRegion,
                                                 EventEnergyScale>& other) const
    {
        if(eventCategory < other.eventCategory) return true;
        if(eventCategory > other.eventCategory) return false;
        if(eventSubCategory < other.eventSubCategory) return true;
        if(eventSubCategory > other.eventSubCategory) return false;
        if(eventRegion < other.eventRegion) return true;
        if(eventRegion > other.eventRegion) return false;
        return eventEnergyScale < other.eventEnergyScale;
    }

    std::string GetName() const
    {
        static const std::string separator = "/";
        std::ostringstream ss;
        ss << eventCategory << separator << eventSubCategory << separator << eventRegion << separator
           << eventEnergyScale;
        return ss.str();
    }
};

typedef FlatAnalyzerDataMetaId<EventCategory, EventSubCategory, EventRegion, EventEnergyScale>
    FlatAnalyzerDataMetaId_noName;

template<>
struct FlatAnalyzerDataMetaId<EventCategory, EventSubCategory, EventEnergyScale> {
    EventCategory eventCategory;
    EventSubCategory eventSubCategory;
    EventEnergyScale eventEnergyScale;

    FlatAnalyzerDataMetaId()
        : eventCategory(EventCategory::Inclusive), eventSubCategory(EventSubCategory::NoCuts),
          eventEnergyScale(EventEnergyScale::Central) {}

    FlatAnalyzerDataMetaId(EventCategory _eventCategory, EventSubCategory _eventSubCategory,
                           EventEnergyScale _eventEnergyScale)
        : eventCategory(_eventCategory), eventSubCategory(_eventSubCategory), eventEnergyScale(_eventEnergyScale) {}

    FlatAnalyzerDataId MakeId(EventRegion eventRegion, const std::string& dataCategoryName) const
    {
        return FlatAnalyzerDataId(eventCategory, eventSubCategory, eventRegion, eventEnergyScale, dataCategoryName);
    }

    FlatAnalyzerDataMetaId_noName MakeMetaId(EventRegion eventRegion) const
    {
        return FlatAnalyzerDataMetaId_noName(eventCategory, eventSubCategory, eventRegion, eventEnergyScale);
    }

    bool operator< (const FlatAnalyzerDataMetaId<EventCategory, EventSubCategory, EventEnergyScale>& other) const
    {
        if(eventCategory < other.eventCategory) return true;
        if(eventCategory > other.eventCategory) return false;
        if(eventSubCategory < other.eventSubCategory) return true;
        if(eventSubCategory > other.eventSubCategory) return false;
        return eventEnergyScale < other.eventEnergyScale;
    }

    std::string GetName() const
    {
        static const std::string separator = "/";
        static const std::string wildcard = "*";
        std::ostringstream ss;
        ss << eventCategory << separator << eventSubCategory << separator << wildcard << separator
           << eventEnergyScale;
        return ss.str();
    }
};

typedef FlatAnalyzerDataMetaId<EventCategory, EventSubCategory, EventEnergyScale>
    FlatAnalyzerDataMetaId_noRegion_noName;

template<typename ...Types>
std::ostream& operator<< (std::ostream& s, const FlatAnalyzerDataMetaId<Types...>& meta_id)
{
    s << meta_id.GetName();
    return s;
}

typedef std::shared_ptr<FlatAnalyzerData> FlatAnalyzerDataPtr;
typedef std::map<FlatAnalyzerDataId, FlatAnalyzerDataPtr> FlatAnalyzerDataMap;

class FlatAnalyzerDataCollection {
public:
    FlatAnalyzerDataCollection(const std::string& outputFileName, bool store)
    {
        if(store)
            outputFile = root_ext::CreateRootFile(outputFileName);

        TH1::SetDefaultSumw2();
        TH1::AddDirectory(kFALSE);
        TH2::AddDirectory(kFALSE);
    }

    FlatAnalyzerData& Get(const FlatAnalyzerDataId& id, Channel channel)
    {
        auto& anaData = anaDataMap[id];
        if(!anaData)
            anaData = MakeAnaData(id, channel);
        return *anaData;
    }

    template<typename GenericEventInfo>
    void Fill(const FlatAnalyzerDataId& id, Channel channel, const GenericEventInfo& eventInfo, double weight)
    {
        auto& anaData = Get(id, channel);
        anaData.Fill(eventInfo, weight);
//        if (id.eventRegion == EventRegion::OS_Isolated && id.eventCategory == EventCategory::TwoJets_TwoBtag &&
//                id.eventSubCategory == EventSubCategory::KinematicFitConvergedWithMassWindow &&
//                id.eventEnergyScale == EventEnergyScale::Central)
//            std::cout << "SIGNAL final eventID: " << eventInfo.event->evt << ", run: " <<
//                         eventInfo.event->run << ", lumi: " << eventInfo.event->lumi << std::endl;
    }

    template<typename GenericEventInfo>
    void FillSubCategories(const FlatAnalyzerDataMetaId_noSub& meta_id, Channel channel,
                           const GenericEventInfo& eventInfo, double weight)
    {
        using namespace cuts::massWindow;

        const double mass_tautau = eventInfo.event->m_sv_MC;
        Fill<GenericEventInfo>(meta_id.MakeId(EventSubCategory::NoCuts), channel, eventInfo, weight);
        const bool inside_mass_window = mass_tautau > m_tautau_low && mass_tautau < m_tautau_high
                && eventInfo.Hbb.M() > m_bb_low && eventInfo.Hbb.M() < m_bb_high;
        if(inside_mass_window)
            Fill<GenericEventInfo>(meta_id.MakeId(EventSubCategory::MassWindow), channel, eventInfo, weight);
        else
            Fill<GenericEventInfo>(meta_id.MakeId(EventSubCategory::OutsideMassWindow), channel, eventInfo, weight);
        if(eventInfo.fitResults.has_valid_mass)
            Fill<GenericEventInfo>(meta_id.MakeId(EventSubCategory::KinematicFitConverged), channel, eventInfo, weight);
        if(eventInfo.fitResults.has_valid_mass && inside_mass_window)
            Fill<GenericEventInfo>(meta_id.MakeId(EventSubCategory::KinematicFitConvergedWithMassWindow), channel, eventInfo, weight);
        if(eventInfo.fitResults.has_valid_mass && !inside_mass_window)
            Fill<GenericEventInfo>(meta_id.MakeId(EventSubCategory::KinematicFitConvergedOutsideMassWindow), channel, eventInfo, weight);
    }

    template<typename GenericEventInfo>
    void FillEnergyScales(const FlatAnalyzerDataMetaId_noSub_noES& meta_id, Channel channel,
                          const GenericEventInfo& eventInfo, double weight, const std::set<EventEnergyScale>& energyScales)
    {
        for (EventEnergyScale energyScale : energyScales)
            FillSubCategories<GenericEventInfo>(meta_id.MakeMetaId(energyScale), channel, eventInfo, weight);
    }

    template<typename GenericEventInfo>
    void FillCentralAndJetRelatedScales(const FlatAnalyzerDataMetaId_noSub_noES& meta_id, Channel channel,
                                        const GenericEventInfo& eventInfo, double weight)
    {
        static const std::set<EventEnergyScale> energyScales =
            { EventEnergyScale::Central, EventEnergyScale::JetUp, EventEnergyScale::JetDown,
              EventEnergyScale::BtagEfficiencyUp, EventEnergyScale::BtagEfficiencyDown,
              EventEnergyScale::BtagFakeUp, EventEnergyScale::BtagFakeDown };
        FillEnergyScales<GenericEventInfo>(meta_id, channel, eventInfo, weight, energyScales);
    }

    template<typename GenericEventInfo>
    void FillAllEnergyScales(const FlatAnalyzerDataMetaId_noSub_noES& meta_id, Channel channel,
                             const GenericEventInfo& eventInfo, double weight)
    {
        FillEnergyScales<GenericEventInfo>(meta_id, channel, eventInfo, weight, AllEventEnergyScales);
    }

private:
    FlatAnalyzerDataPtr MakeAnaData(const FlatAnalyzerDataId& id, Channel channel) const
    {
        if(channel == Channel::ETau || channel == Channel::MuTau) {
            if(id.eventCategory == EventCategory::TwoJets_TwoBtag
                    || id.eventCategory == EventCategory::TwoJets_TwoLooseBtag)
                return Make<FlatAnalyzerData_semileptonic_2tag>(id);
            return Make<FlatAnalyzerData_semileptonic>(id);
        }

        if(channel == Channel::TauTau) {
            if(id.eventCategory == analysis::EventCategory::TwoJets_TwoBtag
                    || id.eventCategory == analysis::EventCategory::TwoJets_TwoLooseBtag)
                return Make<FlatAnalyzerData_tautau_2tag>(id);
            return Make<FlatAnalyzerData_tautau_other_tag>(id);
        }

        throw exception("Can't make analyzer data for %1% channel.") % channel;
    }

    template<typename AnaData>
    FlatAnalyzerDataPtr Make(const FlatAnalyzerDataId& id) const
    {
        const bool fill_all = id.eventEnergyScale == EventEnergyScale::Central;
        if(outputFile) {
            const std::string dir_name = id.GetName();
            return FlatAnalyzerDataPtr(new AnaData(outputFile, dir_name, fill_all));
        }
        return FlatAnalyzerDataPtr(new AnaData(fill_all));
    }

private:
    std::shared_ptr<TFile> outputFile;
    FlatAnalyzerDataMap anaDataMap;
};

class FlatAnalyzerDataCollectionReader {
public:

    typedef std::map<std::string, const root_ext::AbstractHistogram*> HistogramMap;

    FlatAnalyzerDataCollectionReader(const std::string& file_name)
        : file(root_ext::OpenRootFile(file_name)), anaDataCollection("", false) {}

    template<typename Histogram>
    const root_ext::SmartHistogram<Histogram>* GetHistogram(const FlatAnalyzerDataId& id, Channel channel,
                                                      const std::string& name)
    {
        const std::string full_name = id.GetName() + "/" + name;
        if(!histograms.count(full_name)) {
            Histogram* original_histogram = root_ext::TryReadObject<Histogram>(*file, full_name);
            if(!original_histogram)
            histograms[full_name] = nullptr;
            FlatAnalyzerData& anaData = anaDataCollection.Get(id, channel);
            anaData.CreateAll();
            root_ext::SmartHistogram<Histogram>* smart_hist = anaData.GetPtr<Histogram>(name);
            if(!smart_hist)
                throw exception("Histogram '%1%' not found.") % name;
            smart_hist->CopyContent(*original_histogram);
            histograms[full_name] = smart_hist;
        }
        return dynamic_cast< const root_ext::SmartHistogram<Histogram>* >(histograms.at(full_name));
    }

private:
    std::shared_ptr<TFile> file;
    FlatAnalyzerDataCollection anaDataCollection;
    HistogramMap histograms;
};

} // namespace analysis
