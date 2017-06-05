/*! Definition of the analyzer descriptor.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <list>
#include <istream>
#include <ostream>
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/PhysicalValue.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisCategories.h"

namespace analysis {

struct AnalyzerDescriptor {
    std::string name;
    std::vector<std::string> file_paths;
    DataCategoryType categoryType;
    double cross_section;
    double int_lumi;
    std::string channel;

    AnalyzerDescriptor()
        : {}

};

using AnalyzerDescriptorCollection = std::unordered_map<std::string, AnalyzerDescriptor>;

} // namespace analysis
