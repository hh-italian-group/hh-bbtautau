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
#include "Instruments/include/SampleDescriptor.h"

namespace analysis {

namespace sample_merging{

struct SMBinDescriptor {
    std::string name;
    std::string file_path;
    FileType fileType;

    PhysicalValue nu;
    PhysicalValue weight;
    size_t inclusive_integral;

    SMBinDescriptor()
        : nu(0.0, std::numeric_limits<double>::infinity()),
          weight(std::numeric_limits<double>::quiet_NaN()) {}


};

using SMBinDescriptorCollection = std::unordered_map<std::string, SMBinDescriptor>;

} //namespace sample_merging

} // namespace analysis
