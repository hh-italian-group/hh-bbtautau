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
#include "h-tautau/Analysis/include/AnalysisTypes.h"

namespace analysis {

class AnalyzerSetup {
public:

    std::string name;
    double int_lumi{0};
    std::vector<std::string> final_variables;
    bool apply_mass_cut{false};
    std::set<EventEnergyScale> energy_scales;

};

using AnalyzerSetupCollection = std::unordered_map<std::string, AnalyzerSetup>;

class SampleDescriptorBase {
public:

    SampleDescriptorBase() {} //constructor
    SampleDescriptorBase(const SampleDescriptorBase& ) = default; //copy constructor
    virtual ~SampleDescriptorBase(){} //destructor

    SampleDescriptorBase& operator= ( const SampleDescriptorBase& ) = default; //assignment

    std::string name;
    std::string title;
    root_ext::Color color{kBlack};
    bool draw{false};
    std::set<Channel> channels;
    DataCategoryType categoryType;
    std::string datacard_name;

    bool CreateDatacard() const
    {
        return datacard_name.size();
    }

};


class SampleDescriptor : public analysis::SampleDescriptorBase
{
public:
    using SampleDescriptorBase::SampleDescriptorBase; //to inherit the constructor, copy constructor and assignment of parent class
                                                      // otherwise you should declare them as well as parent class, except the destructor

    std::string GetFileName(size_t signal_point) const
    {
        if(signal_point >= GetNSignalPoints())
            throw analysis::exception("Signal point chosen is bigger than the size of signal points.");
        std::string file_name = file_path_pattern;
        for (const auto& signal_point_iter : signal_points){
            const auto& point_prefix = signal_point_iter.first;
            const auto& signal_points_list = signal_point_iter.second;
            const std::string& point_value = signal_points_list.at(signal_point);
            boost::algorithm::replace_all(file_name, "{" + point_prefix + "}", point_value);
        }
        return file_name;
    }

    std::string name;
    std::vector<std::string> file_paths;
    std::string file_path_pattern;
    double cross_section{0};
    std::map<std::string, std::string> signal_points_raw; //mass for resonant, radion or graviton, nodes for non-resonant
    std::map<std::string, std::vector<std::string>> signal_points;
    std::map<std::string, root_ext::Color> draw_ex;
    std::vector<double> norm_sf;
    std::vector<std::string> datacard_name_ex;

    void UpdateSignalPoints()
    {
        signal_points.clear();
        size_t n_points = 0;
        for (const auto& signal_point_iter : signal_points_raw){
            const std::string& point_prefix = signal_point_iter.first;
            const std::string& signal_points_list_str = signal_point_iter.second;
            std::vector<std::string> signal_points_list = SplitValueList(signal_points_list_str);
            if(!signal_points_list.size())
                throw analysis::exception("Empty signal points vector.");
            if(n_points && signal_points_list.size() != n_points)
                throw analysis::exception("signal_points_list has different size from n_points.");
            n_points = signal_points_list.size();

            signal_points[point_prefix] = signal_points_list;
        }

    }

    size_t GetNSignalPoints() const { return signal_points.size() ? signal_points.begin()->second.size() : 0; }

};

using SampleDescriptorCollection = std::unordered_map<std::string, SampleDescriptor>;

class CombineSampleDescriptor : public analysis::SampleDescriptorBase
{
public:

    using SampleDescriptorBase::SampleDescriptorBase;   //to inherit the constructor, copy constructor and assignment of parent class
                                                        // otherwise you should declare them as well as parent class, except the destructor

    std::string name;
    std::vector<std::string> sample_descriptors;

};

using CombineSampleDescriptorCollection = std::unordered_map<std::string, CombineSampleDescriptor>;

} // namespace analysis

