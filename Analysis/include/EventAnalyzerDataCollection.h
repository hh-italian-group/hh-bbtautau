/*! Collection of histogram containers for event analyzers.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <tuple>
#include "EventAnalyzerData.h"

namespace analysis {

//namespace detail {
//template<size_t n, typename = typename std::enable_if< n < TupleSize >::type>
//std::string GetEventAnalyzerDataIdSubName() const
//{
//    static const std::string separator = "/", wildcard = "*";
//    const auto& value = std::get<n>(id_tuple);
//    std::string name = value.is_initialized() ? ToString(*value) : wildcard;
//    const std::string others_name = GetSubName<n+1>();
//    return others_name.size() ? name + separator + others_name : name;
//}

//template<size_t n, typename = typename std::enable_if< n == TupleSize >::type>
//std::string GetSubName() const { return ""; }

//} // namespace detail

struct EventAnalyzerDataId {
public:
    template<typename T> using Optional = boost::optional<T>;
    using Tuple = std::tuple<Optional<EventCategory>, Optional<EventSubCategory>, Optional<EventRegion>,
                             Optional<EventEnergyScale>, Optional<Dataset>>;
    static constexpr size_t TupleSize = std::tuple_size<Tuple>::value;

    EventAnalyzerDataId() {}
    template<typename ...Args>
    EventAnalyzerDataId(Args&&... args) { Initialize(std::forward<Args>(args)...); }

    template<typename T>
    bool Has() const { return std::get<Optional<T>>(id_tuple).is_initialized(); }
    template<typename T>
    void Set(const T& value) { std::get<Optional<T>>(id_tuple) = value; }
    template<typename T>
    const T& Get() const
    {
        if(!std::get<Optional<T>>(id_tuple).is_initialized())
            throw exception("%1% is not specified in EventAnalyzerDataId = '%2%'.") % typeid(T).name() % *this;
        return *std::get<Optional<T>>(id_tuple);
    }

    bool operator< (const EventAnalyzerDataId& other) const { return id_tuple < other.id_tuple; }
    std::string GetName() const { return GetSubName(std::make_index_sequence<TupleSize>{}); }
    bool IsComplete() const { return ItemsAreInitialized(std::make_index_sequence<TupleSize>{}); }

private:
    void Initialize() {}

    template<typename T, typename ...Args>
    void Initialize(T&& value, Args&&... other_values)
    {
        using ValueType = typename std::remove_const<typename std::remove_reference<T>::type>::type;
        std::get<Optional<ValueType>>(id_tuple) = value;
        Initialize(std::forward<Args>(other_values)...);
    }

    template<size_t n>
    std::string ItemToString() const
    {
        static const std::string wildcard = "*";
        const auto& value = std::get<n>(id_tuple);
        return value.is_initialized() ? ToString(*value) : wildcard;
    }

    template<size_t ...n>
    std::string GetSubName(std::index_sequence<n...>) const
    {
        static const std::string separator = "/";
        const std::vector<std::string> item_names = { ItemToString<n>()... };
        std::string sub_name;
        for(const auto& name : item_names)
            sub_name += name + separator;
        if(sub_name.size())
            sub_name.erase(sub_name.size() - 1);
        return sub_name;
    }

    template<size_t ...n>
    bool ItemsAreInitialized(std::index_sequence<n...>) const
    {
        const std::vector<bool> is_initialized = { std::get<n>(id_tuple).is_initialized()... };
        for(bool init : is_initialized)
            if(!init) return false;
        return true;
    }

private:
    Tuple id_tuple;
};

std::ostream& operator<< (std::ostream& s, const EventAnalyzerDataId& id)
{
    s << id.GetName();
    return s;
}

class BaseEventAnalyzerDataCollection {
public:
    virtual ~BaseEventAnalyzerDataCollection() {}
    virtual BaseEventAnalyzerData& GetBase(const EventAnalyzerDataId& id) = 0;

};

template<typename AnaData>
class EventAnalyzerDataCollection : public BaseEventAnalyzerDataCollection {
public:
    using Data = AnaData;
    using DataPtr = std::shared_ptr<Data>;
    using DataId = EventAnalyzerDataId;
    using DataMap = std::map<DataId, DataPtr>;

    explicit EventAnalyzerDataCollection(const std::string& inputFileName) :
        file(root_ext::OpenRootFile(inputFileName)), readMode(true)
    {}

    EventAnalyzerDataCollection(const std::string& outputFileName, bool store) :
        readMode(false)
    {
        if(store)
            file = root_ext::CreateRootFile(outputFileName);
    }

    Data& Get(const DataId& id)
    {
        auto& anaData = anaDataMap[id];
        if(!anaData)
            anaData = Make(id);
        return *anaData;
    }

    virtual BaseEventAnalyzerData& GetBase(const EventAnalyzerDataId& id) override { return Get(id); }

    template<typename EventInfo>
    void Fill(const EventAnalyzerDataId& id, EventInfo& event, double weight)
    {
        Get(id).Fill(event, weight);
    }

private:
    DataPtr Make(const EventAnalyzerDataId& id) const
    {
        const bool fill_all = id.Get<EventEnergyScale>() == EventEnergyScale::Central;
        if(file) {
            if(!id.IsComplete())
                throw exception("EventAnalyzerDataId is not complete.");
            const std::string dir_name = id.GetName();
            return std::make_shared<Data>(file, dir_name, fill_all, readMode);
        }
        return std::make_shared<Data>(fill_all);
    }

private:
    std::shared_ptr<TFile> file;
    bool readMode;
    DataMap anaDataMap;
};

} // namespace analysis
