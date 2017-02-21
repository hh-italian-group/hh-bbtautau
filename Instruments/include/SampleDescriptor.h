/*! Definition of the file descriptor.
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
};


} // namespace analysis
