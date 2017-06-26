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
    std::set<Channel> channels;
    DataCategoryType categoryType;
    std::string datacard_name;

    bool CreateDatacard()
    {
        return datacard_name.size();
    }

};


class SampleDescriptor : public analysis::SampleDescriptorBase
{
public:
    SampleDescriptor() : cross_section(0) {}

    std::map<std::string, std::vector<std::string>> GetMapOfVectorOfString() const
    {
        std::map<std::string, std::vector<std::string>> signal_points_map;
        for (auto signal_point_iter : signal_points){
            std::string point_prefix = signal_point_iter.first;
            std::string signal_points_list_str = signal_point_iter.second;
            std::vector<std::string> signal_points_list = SplitValueList(signal_points_list_str);
            signal_points_map[point_prefix] = signal_points_list;
        }
        return signal_points_map;
    }

    std::string GetFileName(size_t signal_point) const
    {
        std::string file_name = file_path_pattern;
        for (auto signal_point_iter : signal_points){
            std::string point_prefix = signal_point_iter.first;
            std::string signal_points_list_str = signal_point_iter.second;
            std::vector<std::string> signal_points_list = SplitValueList(signal_points_list_str);
            std::string point_value = signal_points_list.at(signal_point);
            boost::algorithm::replace_all(file_name, point_prefix, point_value);
        }
        return file_name;
    }

    std::string name;
    std::vector<std::string> file_paths;
    std::string file_path_pattern;
    double cross_section;
    std::map<std::string, std::string> signal_points; //mass for resonant, radion or graviton, nodes for non-resonant
    std::map<std::string, std::string> draw_ex; //should be std::map<std::string, Color> not working
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

