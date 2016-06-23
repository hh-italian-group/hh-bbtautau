/*! Collection of histogram containers for event analyzers.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "EventAnalyzerData.h"
#include "h-tautau/Analysis/include/Htautau_2015.h"
#include "custom_cuts.h"

namespace analysis {

struct EventAnalyzerDataId;
using EventAnalyzerDataIdSet = std::set<EventAnalyzerDataId>;

struct EventAnalyzerDataId {
    EventCategory eventCategory;
    EventSubCategory eventSubCategory;
    EventRegion eventRegion;
    EventEnergyScale eventEnergyScale;
    std::string dataCategoryName;

    EventAnalyzerDataId()
        : eventCategory(EventCategory::Inclusive), eventSubCategory(EventSubCategory::NoCuts),
          eventRegion(EventRegion::OS_Isolated), eventEnergyScale(EventEnergyScale::Central)
    {}

    EventAnalyzerDataId(EventCategory _eventCategory, EventSubCategory _eventSubCategory, EventRegion _eventRegion,
                       EventEnergyScale _eventEnergyScale, const std::string& _dataCategoryName)
        : eventCategory(_eventCategory), eventSubCategory(_eventSubCategory), eventRegion(_eventRegion),
          eventEnergyScale(_eventEnergyScale), dataCategoryName(_dataCategoryName)
    {}

    bool operator< (const EventAnalyzerDataId& other) const
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

std::ostream& operator<< (std::ostream& s, const EventAnalyzerDataId& id)
{
    s << id.GetName();
    return s;
}

template<typename ...Types>
struct EventAnalyzerDataMetaId;

template<>
struct EventAnalyzerDataMetaId<EventCategory, EventRegion, EventEnergyScale, std::string> {
    EventCategory eventCategory;
    EventRegion eventRegion;
    EventEnergyScale eventEnergyScale;
    std::string dataCategoryName;

    EventAnalyzerDataMetaId()
        : eventCategory(EventCategory::Inclusive), eventRegion(EventRegion::OS_Isolated),
          eventEnergyScale(EventEnergyScale::Central) {}

    EventAnalyzerDataMetaId(EventCategory _eventCategory, EventRegion _eventRegion, EventEnergyScale _eventEnergyScale,
                           const std::string& _dataCategoryName)
        : eventCategory(_eventCategory), eventRegion(_eventRegion),
          eventEnergyScale(_eventEnergyScale), dataCategoryName(_dataCategoryName)  {}

    EventAnalyzerDataId MakeId(EventSubCategory eventSubCategory) const
    {
        return EventAnalyzerDataId(eventCategory, eventSubCategory, eventRegion, eventEnergyScale, dataCategoryName);
    }

    bool operator<(const EventAnalyzerDataMetaId<EventCategory, EventRegion, EventEnergyScale, std::string>& other) const
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

using EventAnalyzerDataMetaId_noSub = EventAnalyzerDataMetaId<EventCategory, EventRegion, EventEnergyScale, std::string>;

template<>
struct EventAnalyzerDataMetaId<EventCategory, EventRegion, std::string> {
    EventCategory eventCategory;
    EventRegion eventRegion;
    std::string dataCategoryName;

    EventAnalyzerDataMetaId()
        : eventCategory(EventCategory::Inclusive), eventRegion(EventRegion::OS_Isolated) {}

    EventAnalyzerDataMetaId(EventCategory _eventCategory, EventRegion _eventRegion, const std::string& _dataCategoryName)
        : eventCategory(_eventCategory), eventRegion(_eventRegion), dataCategoryName(_dataCategoryName)  {}

    EventAnalyzerDataId MakeId(EventSubCategory eventSubCategory, EventEnergyScale eventEnergyScale) const
    {
        return EventAnalyzerDataId(eventCategory, eventSubCategory, eventRegion, eventEnergyScale, dataCategoryName);
    }

    EventAnalyzerDataMetaId_noSub MakeMetaId(EventEnergyScale eventEnergyScale) const
    {
        return EventAnalyzerDataMetaId_noSub(eventCategory, eventRegion, eventEnergyScale, dataCategoryName);
    }

    bool operator< (const EventAnalyzerDataMetaId<EventCategory, EventRegion, std::string>& other) const
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

using EventAnalyzerDataMetaId_noSub_noES = EventAnalyzerDataMetaId<EventCategory, EventRegion, std::string>;

template<>
struct EventAnalyzerDataMetaId<EventCategory, EventSubCategory, EventRegion, EventEnergyScale> {
    EventCategory eventCategory;
    EventSubCategory eventSubCategory;
    EventRegion eventRegion;
    EventEnergyScale eventEnergyScale;

    EventAnalyzerDataMetaId()
        : eventCategory(EventCategory::Inclusive), eventSubCategory(EventSubCategory::NoCuts),
          eventRegion(EventRegion::OS_Isolated), eventEnergyScale(EventEnergyScale::Central) {}

    EventAnalyzerDataMetaId(EventCategory _eventCategory, EventSubCategory _eventSubCategory, EventRegion _eventRegion,
                           EventEnergyScale _eventEnergyScale)
        : eventCategory(_eventCategory), eventSubCategory(_eventSubCategory), eventRegion(_eventRegion),
          eventEnergyScale(_eventEnergyScale) {}

    EventAnalyzerDataId MakeId(const std::string& dataCategoryName) const
    {
        return EventAnalyzerDataId(eventCategory, eventSubCategory, eventRegion, eventEnergyScale, dataCategoryName);
    }

    bool operator< (const EventAnalyzerDataMetaId<EventCategory, EventSubCategory, EventRegion,
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

using EventAnalyzerDataMetaId_noName = EventAnalyzerDataMetaId<EventCategory, EventSubCategory, EventRegion, EventEnergyScale>;

template<>
struct EventAnalyzerDataMetaId<EventCategory, EventSubCategory, EventEnergyScale> {
    EventCategory eventCategory;
    EventSubCategory eventSubCategory;
    EventEnergyScale eventEnergyScale;

    EventAnalyzerDataMetaId()
        : eventCategory(EventCategory::Inclusive), eventSubCategory(EventSubCategory::NoCuts),
          eventEnergyScale(EventEnergyScale::Central) {}

    EventAnalyzerDataMetaId(EventCategory _eventCategory, EventSubCategory _eventSubCategory,
                           EventEnergyScale _eventEnergyScale)
        : eventCategory(_eventCategory), eventSubCategory(_eventSubCategory), eventEnergyScale(_eventEnergyScale) {}

    EventAnalyzerDataId MakeId(EventRegion eventRegion, const std::string& dataCategoryName) const
    {
        return EventAnalyzerDataId(eventCategory, eventSubCategory, eventRegion, eventEnergyScale, dataCategoryName);
    }

    EventAnalyzerDataMetaId_noName MakeMetaId(EventRegion eventRegion) const
    {
        return EventAnalyzerDataMetaId_noName(eventCategory, eventSubCategory, eventRegion, eventEnergyScale);
    }

    bool operator< (const EventAnalyzerDataMetaId<EventCategory, EventSubCategory, EventEnergyScale>& other) const
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

using EventAnalyzerDataMetaId_noRegion_noName = EventAnalyzerDataMetaId<EventCategory, EventSubCategory, EventEnergyScale>;

template<typename ...Types>
std::ostream& operator<< (std::ostream& s, const EventAnalyzerDataMetaId<Types...>& meta_id)
{
    s << meta_id.GetName();
    return s;
}

using EventAnalyzerDataPtr = std::shared_ptr<BaseEventAnalyzerData>;
using EventAnalyzerDataMap = std::map<EventAnalyzerDataId, EventAnalyzerDataPtr>;

class EventAnalyzerDataCollection {
public:
    EventAnalyzerDataCollection(const std::string& outputFileName, bool store)
    {
        if(store)
            outputFile = root_ext::CreateRootFile(outputFileName);

        TH1::SetDefaultSumw2();
        TH1::AddDirectory(kFALSE);
        TH2::AddDirectory(kFALSE);
    }

    template<typename FirstLeg>
    EventAnalyzerData<FirstLeg>& Get(const EventAnalyzerDataId& id)
    {
        auto& anaData = anaDataMap[id];
        if(!anaData)
            anaData = MakeAnaData<FirstLeg>(id);
        return *dynamic_cast<EventAnalyzerData<FirstLeg>*>(anaData.get());
    }

    template<typename EventInfo>
    void Fill(const EventAnalyzerDataId& id, EventInfo& event, double weight)
    {
        auto& anaData = Get<typename EventInfo::FirstLeg>(id);
        anaData.Fill(event, weight);
    }

private:
    template<typename FirstLeg>
    EventAnalyzerDataPtr MakeAnaData(const EventAnalyzerDataId& id) const
    {
        const Channel channel = ChannelInfo::IdentifyChannel<FirstLeg>();
        if(channel == Channel::ETau || channel == Channel::MuTau) {
            if(id.eventCategory == EventCategory::TwoJets_TwoBtag
                    || id.eventCategory == EventCategory::TwoJets_TwoLooseBtag)
                return Make<EventAnalyzerData_semileptonic_2tag<FirstLeg>>(id);
            return Make<EventAnalyzerData<FirstLeg>>(id);
        }

        if(channel == Channel::TauTau) {
            if(id.eventCategory == analysis::EventCategory::TwoJets_TwoBtag
                    || id.eventCategory == analysis::EventCategory::TwoJets_TwoLooseBtag)
                return Make<EventAnalyzerData_tautau_2tag>(id);
            return Make<EventAnalyzerData_tautau_other_tag>(id);
        }

        throw exception("Can't make analyzer data for %1% channel.") % channel;
    }

    template<typename AnaData>
    EventAnalyzerDataPtr Make(const EventAnalyzerDataId& id) const
    {
        const bool fill_all = id.eventEnergyScale == EventEnergyScale::Central;
        if(outputFile) {
            const std::string dir_name = id.GetName();
            return EventAnalyzerDataPtr(new AnaData(outputFile, dir_name, fill_all));
        }
        return EventAnalyzerDataPtr(new AnaData(fill_all));
    }

private:
    std::shared_ptr<TFile> outputFile;
    EventAnalyzerDataMap anaDataMap;
};

class EventAnalyzerDataCollectionReader {
public:
    using HistogramMap = std::map<std::string, const root_ext::AbstractHistogram*>;

    EventAnalyzerDataCollectionReader(const std::string& file_name)
        : file(root_ext::OpenRootFile(file_name)), anaDataCollection("", false) {}

    template<typename FirstLeg, typename Histogram>
    const root_ext::SmartHistogram<Histogram>* GetHistogram(const EventAnalyzerDataId& id, const std::string& name)
    {
        const std::string full_name = id.GetName() + "/" + name;
        if(!histograms.count(full_name)) {
            Histogram* original_histogram = root_ext::TryReadObject<Histogram>(*file, full_name);
            if(!original_histogram)
            histograms[full_name] = nullptr;
            auto& anaData = anaDataCollection.Get<FirstLeg>(id);
            anaData.CreateAll();
            root_ext::SmartHistogram<Histogram>* smart_hist = anaData. template GetPtr<Histogram>(name);
            if(!smart_hist)
                throw exception("Histogram '%1%' not found.") % name;
            smart_hist->CopyContent(*original_histogram);
            histograms[full_name] = smart_hist;
        }
        return dynamic_cast< const root_ext::SmartHistogram<Histogram>* >(histograms.at(full_name));
    }

private:
    std::shared_ptr<TFile> file;
    EventAnalyzerDataCollection anaDataCollection;
    HistogramMap histograms;
};

} // namespace analysis
