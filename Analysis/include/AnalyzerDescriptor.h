/*! Definition of the analyzer descriptor.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <list>
#include <istream>
#include <ostream>
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Print/include/PlotPrimitives.h"
#include "AnalysisTools/Core/include/PhysicalValue.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisCategories.h"

namespace analysis {

struct AnalyzerDescriptor {
    std::string name;
    std::map<std::string, double> file_cross_section_map;
    DataCategoryType categoryType;
//    double cross_section;
    double int_lumi;
    std::string channel;
    root_ext::Color color;

    AnalyzerDescriptor()
        : color(kBlack) {}

};

using AnalyzerDescriptorCollection = std::unordered_map<std::string, AnalyzerDescriptor>;

} // namespace analysis
