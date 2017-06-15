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
#include "SampleDescriptorBase.h"

namespace analysis {

class SampleDescriptor : public analysis::SampleDescriptorBase
{
public:
    SampleDescriptor() {}

private:
    std::string name;
    std::vector<std::string> file_paths;
    double cross_section;
    std::string weight_file;
    std::vector<std::string> listSignalPoints; //mass for resonant, radion or graviton, nodes for non-resonant

};

using SampleDescriptorCollection = std::unordered_map<std::string, SampleDescriptor>;

} // namespace analysis
