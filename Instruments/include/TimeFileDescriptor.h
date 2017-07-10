/*! Definition of the file descriptor for TTbar sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <list>
#include <istream>
#include <ostream>
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/PhysicalValue.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/ConfigReader.h"

namespace analysis {

struct TimeFileDescriptor {
    std::string name;
    std::vector<std::string> file_paths;


    TimeFileDescriptor(){}


};

using TimeFileDescriptorCollection = std::unordered_map<std::string, TimeFileDescriptor>;

} // namespace analysis
