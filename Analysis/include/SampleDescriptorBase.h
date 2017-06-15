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

class SampleDescriptorBase {
public:
    SampleDescriptorBase()
        : color(kBlack), draw(false), create_hist(false) {}

private:
    std::string name;
    std::string title;
    root_ext::Color color;
    bool draw;
    std::string channel;
    DataCategoryType categoryType;
    bool create_hist;
    std::string datacard_name;

};

using SampleDescriptorBaseCollection = std::unordered_map<std::string, SampleDescriptorBase>;

} // namespace analysis
