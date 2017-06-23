
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
#include "h-tautau/Analysis/include/AnalysisTypes.h"

namespace analysis {

class AnalyzerSetup { //modifica
public:
    AnalyzerSetup()
        : apply_mass_cut(false) {}

    std::string name;
    double int_lumi;
    std::vector<std::string> final_variables;
    bool apply_mass_cut;
    std::vector<EventEnergyScale> energy_scales;

};

using AnalyzerSetupCollection = std::unordered_map<std::string, AnalyzerSetup>;

} // namespace analysis
