/*! Definition of identifier of event analyzer data container.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "hh-bbtautau/Analysis/include/EventAnalyzerDataId.h"

namespace analysis {

bool EventAnalyzerDataId::operator< (const EventAnalyzerDataId& other) const
{
    return id_tuple < other.id_tuple;
}

bool EventAnalyzerDataId::operator== (const EventAnalyzerDataId& other) const
{
    return id_tuple == other.id_tuple;
}

bool EventAnalyzerDataId::operator!= (const EventAnalyzerDataId& other) const
{
    return id_tuple != other.id_tuple;
}

std::string EventAnalyzerDataId::GetName(const std::string& separator) const
{
    return GetSubName(std::make_index_sequence<TupleSize>{}, separator);
}

bool EventAnalyzerDataId::IsComplete() const
{
    return ItemsAreInitialized(std::make_index_sequence<TupleSize>{});
}

EventAnalyzerDataId EventAnalyzerDataId::Parse(const std::string& str)
{
    static const std::string separator = "/";
    try {
        std::vector<std::string> item_strings = SplitValueList(str, true, separator, false);
        if(item_strings.size() != std::tuple_size<Tuple>::value)
            throw exception("Number of elements != %1%.") % std::tuple_size<Tuple>::value;
        EventAnalyzerDataId id;
        id.ParseAllItems(item_strings, std::make_index_sequence<TupleSize>{});
        return id;
    } catch(std::exception& e) {
        throw exception("Invalid EventAnalyzerDataId = '%1%'. %2%") % str % e.what();
    }
}

std::ostream& operator<<(std::ostream& s, const EventAnalyzerDataId& id)
{
    s << id.GetName();
    return s;
}

std::istream& operator>>(std::istream& s, EventAnalyzerDataId& id)
{
    std::string str;
    s >> str;
    id = EventAnalyzerDataId::Parse(str);
    return s;
}

} // namespace analysis
