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
#include "SampleDescriptor.h"

namespace analysis {

class CombineSampleDescriptor : public analysis::SampleDescriptor
{
public:
    CombineSampleDescriptor() {}

private:
    std::string name;
    std::vector<*SampleDescriptor> sample_descriptors;

};

using CombineSampleDescriptorCollection = std::unordered_map<std::string, CombineSampleDescriptor>;

} // namespace analysis
