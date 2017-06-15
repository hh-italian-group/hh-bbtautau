
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

class AnalyzerDescriptor {
public:
    AnalyzerDescriptor()
        : apply_mass_cut(false) {}

private:
    std::string name;
    double int_lumi;
    std::vector<std::string> final_variables;
    bool apply_mass_cut;

};

using AnalyzerDescriptorCollection = std::unordered_map<std::string, AnalyzerDescriptor>;

} // namespace analysis
