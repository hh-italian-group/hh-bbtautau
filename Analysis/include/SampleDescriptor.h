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
#include <boost/algorithm/string.hpp>

namespace analysis {

class SampleDescriptor : public analysis::SampleDescriptorBase
{
public:
    SampleDescriptor() {}

    std::string GetFileName(size_t signal_point) const
    {
        std::string file_name;
        for (auto signal_point_iter : listSignalPoints){
            std::string point_prefix = signal_point_iter.first;
            std::vector<std::string> signal_points = signal_point_iter.second;
            std::string point_value = signal_points.at(signal_point);
            file_name = boost::algorithm::replace_all_copy(file_path_pattern, point_prefix, point_value);
//            file_name = file_path_pattern;
        }
        return file_name;
    }

    std::string name;
    std::vector<std::string> file_paths;
    std::string file_path_pattern;
    double cross_section;
    std::string weight_file;
    std::map<std::string, std::vector<std::string>> listSignalPoints; //mass for resonant, radion or graviton, nodes for non-resonant
    std::map<std::string, root_ext::Color> draw_ex;
    std::vector<double> norm_sf;
    std::vector<std::string> datacard_name_ex;
};

using SampleDescriptorCollection = std::unordered_map<std::string, SampleDescriptor>;

} // namespace analysis
