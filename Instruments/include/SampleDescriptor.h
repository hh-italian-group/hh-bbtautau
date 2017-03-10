/*! Definition of the sample descriptor used for sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <list>
#include <istream>
#include <ostream>
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/PhysicalValue.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"

namespace analysis {

namespace sample_merging{
#ifndef __APPLE__
    ENUM_OSTREAM_OPERATOR()
    ENUM_ISTREAM_OPERATORS()
#endif
enum class FileType { inclusive, exclusive };
ENUM_NAMES(FileType) = {
    { FileType::inclusive, "inclusive" },
    { FileType::exclusive, "exclusive" }
};

template<typename BinDescriptor, typename GenMap>
struct SampleDescriptor {
    BinDescriptor bin;
    GenMap gen_counts;

    size_t Integral() const
    {
        size_t total_n_events = 0;
        for (const auto& genEventCount : gen_counts){
            const size_t nevents = genEventCount.second;
            total_n_events += nevents;
        }
        return total_n_events;
    }
    size_t Integral(const BinDescriptor& binDescriptor_range) const
    {
        size_t totalEvents = 0;
        for (const auto& genEventCount : gen_counts){
            const auto& genId = genEventCount.first;
            const size_t nevents = genEventCount.second;
            if(binDescriptor_range.Contains(genId))
                totalEvents += nevents;
        }
        return totalEvents;
    }

    double Integral_double(const BinDescriptor& binDescriptor_range) const
    {
        double totalEvents = 0;
        for (const auto& genEventCount : gen_counts){
            const auto& genId = genEventCount.first;
            const double nevents = genEventCount.second;
            if(binDescriptor_range.Contains(genId))
                totalEvents += nevents;
        }
        return totalEvents;
    }
};

} // namespace sample_merging

} // namespace analysis
