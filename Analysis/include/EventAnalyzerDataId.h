/*! Definition of identifier of event analyzer data container.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisCategories.h"

namespace analysis {

namespace detail {

template<typename CollectionItem>
struct EventAnalyzerDataId_IdElementAccessor {
    using ElementType = CollectionItem;
    static ElementType Get(const CollectionItem& item) { return item; }
};

template<typename Key, typename Value>
struct EventAnalyzerDataId_IdElementAccessor<std::pair<const Key, Value>> {
    using ElementType = Key;
    static ElementType Get(const std::pair<Key, Value>& item) { return item.first; }
};

} // namesapce detail

struct EventAnalyzerDataId {
public:
    template<typename T> using Optional = boost::optional<T>;
    using Tuple = std::tuple<Optional<EventCategory>, Optional<EventSubCategory>, Optional<EventRegion>,
                             Optional<UncertaintySource>, Optional<UncertaintyScale>, Optional<Dataset>>;
    static constexpr size_t TupleSize = std::tuple_size<Tuple>::value;

    EventAnalyzerDataId() {}
    EventAnalyzerDataId(const EventAnalyzerDataId& other) = default;
    EventAnalyzerDataId(EventAnalyzerDataId&& other) = default;
    EventAnalyzerDataId& operator=(const EventAnalyzerDataId&) = default;
    EventAnalyzerDataId(EventAnalyzerDataId& other) { *this = other; }

    template<typename ...Args>
    EventAnalyzerDataId(Args&&... args) { Initialize(std::forward<Args>(args)...); }

    template<typename T>
    bool Has() const { return std::get<Optional<T>>(id_tuple).is_initialized(); }
    template<typename T>
    EventAnalyzerDataId Set(const T& value) const
    {
        EventAnalyzerDataId new_id(*this);
        std::get<Optional<T>>(new_id.id_tuple) = value;
        return new_id;
    }
    template<typename T>
    const T& Get() const
    {
        if(!std::get<Optional<T>>(id_tuple).is_initialized())
            throw exception("%1% is not specified in EventAnalyzerDataId = '%2%'.") % typeid(T).name() % *this;
        return *std::get<Optional<T>>(id_tuple);
    }

    bool operator< (const EventAnalyzerDataId& other) const;
    std::string GetName(const std::string& separator = "/") const;
    bool IsComplete() const;

    template<typename ...Collections>
    static std::vector<EventAnalyzerDataId> MetaLoop(Collections&&... cols)
    {
        std::vector<EventAnalyzerDataId> result;
        AddMetaLoopLevels(result, std::forward<Collections>(cols)...);
        return result;
    }

    static EventAnalyzerDataId Parse(const std::string& str);

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
    std::string GetSubName(std::index_sequence<n...>, const std::string& separator) const
    {
        const std::vector<std::string> item_names = { ItemToString<n>()... };
        std::string sub_name;
        for(const auto& name : item_names)
            sub_name += name + separator;
        if(sub_name.size())
            sub_name.erase(sub_name.size() - 1);
        return sub_name;
    }

    template<size_t n>
    bool ParseItem(const std::vector<std::string>& item_strings)
    {
        using TupleElement = typename std::tuple_element<n, Tuple>::type;
        using ValueType = typename TupleElement::value_type;
        static const std::string wildcard = "*";
        const std::string& str = item_strings.at(n);
        if(str != wildcard)
            std::get<n>(id_tuple) = ::analysis::Parse<ValueType>(str);
        return true;
    }

    template<size_t ...n>
    void ParseAllItems(const std::vector<std::string>& item_strings, std::index_sequence<n...>)
    {
        const std::vector<bool> res = { ParseItem<n>(item_strings)... };
    }

    template<size_t ...n>
    bool ItemsAreInitialized(std::index_sequence<n...>) const
    {
        const std::vector<bool> is_initialized = { std::get<n>(id_tuple).is_initialized()... };
        for(bool init : is_initialized)
            if(!init) return false;
        return true;
    }

    static void AddMetaLoopLevels(std::vector<EventAnalyzerDataId>&) {}

    template<typename Collection, typename ...OtherCollections>
    static void AddMetaLoopLevels(std::vector<EventAnalyzerDataId>& result, Collection&& collection,
                         OtherCollections&&... other_collections)
    {
        using Item = typename std::remove_reference<Collection>::type::value_type;
        using ElementAccessor = typename detail::EventAnalyzerDataId_IdElementAccessor<Item>;
        using ElementType = typename ElementAccessor::ElementType;

        const size_t N = std::max<size_t>(result.size(), 1) * collection.size();
        std::vector<EventAnalyzerDataId> new_result;
        new_result.reserve(N);
        if(!result.size()) {
            for(const auto& item : collection) {
                const auto& id_element = ElementAccessor::Get(item);
                new_result.emplace_back(id_element);
            }
        } else {
            for(const EventAnalyzerDataId& ref_id : result) {
                if(ref_id.Has<ElementType>())
                    throw exception("Duplicated element type '%1%' in meta loop.") % typeid(ElementType).name();
                for(const auto& item : collection) {
                    const auto& id_element = ElementAccessor::Get(item);
                    new_result.push_back(ref_id.Set(id_element));
                }
            }
        }

        result = std::move(new_result);
        AddMetaLoopLevels(result, std::forward<OtherCollections>(other_collections)...);
    }

private:
    Tuple id_tuple;
};

std::ostream& operator<<(std::ostream& s, const EventAnalyzerDataId& id);
std::istream& operator>>(std::istream& s, EventAnalyzerDataId& id);

} // namespace analysis
