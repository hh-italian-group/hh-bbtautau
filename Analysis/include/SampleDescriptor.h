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
#include <boost/algorithm/string.hpp>
#include "AnalysisCategories.h"

namespace analysis {

class SampleDescriptorBase {
public:
    SampleDescriptorBase()
        : color(kBlack), draw(false) {}

    std::string name;
    std::string title;
    root_ext::Color color;
    bool draw;
    std::vector<Channel> channels;
    DataCategoryType categoryType;
    std::string datacard_name;

    bool CreateDatacard()
    {
        if (!datacard_name.size())
            return false;
        return true;
    }

};


class SampleDescriptor : public analysis::SampleDescriptorBase
{
public:
    SampleDescriptor() {}

    std::string GetFileName(size_t signal_point) const
    {
        std::string file_name = file_path_pattern;
        for (auto signal_point_iter : signal_points){
            std::string point_prefix = signal_point_iter.first;
            std::vector<std::string> signal_points_list = signal_point_iter.second;
            std::string point_value = signal_points_list.at(signal_point);
            boost::algorithm::replace_all(file_name, point_prefix, point_value);
        }
        return file_name;
    }

    std::string name;
    std::vector<std::string> file_paths;
    std::string file_path_pattern;
    double cross_section;
    std::map<std::string, std::vector<std::string>> signal_points; //mass for resonant, radion or graviton, nodes for non-resonant
    std::map<std::string, std::string> draw_ex;
    std::vector<double> norm_sf;
    std::vector<std::string> datacard_name_ex;
};

using SampleDescriptorCollection = std::unordered_map<std::string, SampleDescriptor>;

class CombineSampleDescriptor : public analysis::SampleDescriptorBase
{
public:

    CombineSampleDescriptor() {}

    std::string name;
    std::vector<std::string> sample_descriptors;

};

using CombineSampleDescriptorCollection = std::unordered_map<std::string, CombineSampleDescriptor>;

} // namespace analysis

