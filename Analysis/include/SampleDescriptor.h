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

struct AnalyzerSetup {
    std::string name;
    double int_lumi{0};
    std::vector<std::string> final_variables;
    bool apply_mass_cut{false}, apply_os_cut{true}, apply_iso_cut{true};
    std::set<EventEnergyScale> energy_scales;
    std::vector<std::string> data, signals, backgrounds, cmb_samples;
    std::vector<std::string> draw_sequence;
};

using AnalyzerSetupCollection = std::unordered_map<std::string, AnalyzerSetup>;

struct SampleDescriptorBase {
    struct Point {
        std::string name, full_name, title, file_path, datacard_name;
        double norm_sf{1}, datacard_sf{1}, draw_sf{1};
        bool draw{false};
        root_ext::Color color{kBlack};
        std::vector<double> param_values;
    };
    using PointCollection = std::vector<Point>;

    std::string name;
    std::string title;
    root_ext::Color color{kBlack};
    double datacard_sf{1}, draw_sf{1};
    std::set<Channel> channels;
    SampleType sampleType{SampleType::MC};
    std::string datacard_name;

    PointCollection working_points;

    SampleDescriptorBase() {} //constructor
    SampleDescriptorBase(const SampleDescriptorBase& ) = default; //copy constructor
    virtual ~SampleDescriptorBase(){} //destructor

    SampleDescriptorBase& operator= ( const SampleDescriptorBase& ) = default; //assignment

    bool HasDatacardName() const { return datacard_name.size(); }

    virtual void CreateWorkingPoints()
    {
        working_points.clear();
        Point point;
        point.full_name = name;
        point.title = title.size() ? title : point.full_name;
        ReplacePatternItem(point.title, "factor", draw_sf);
        point.datacard_name = datacard_name;
        point.draw_sf = draw_sf;
        point.datacard_sf = datacard_sf;
        point.draw = true;
        point.color = color;
        working_points.push_back(point);
    }

    virtual std::map<std::string, size_t> GetModelParameterNames() const { return {}; }

protected:
    template<typename T>
    static void ReplacePatternItem(std::string& str, const std::string& pattern_item, const T& value)
    {
        const auto value_str = ToString(value);
        boost::algorithm::replace_all(str, "{" + pattern_item + "}", value_str);
    }
};

struct SampleDescriptor : SampleDescriptorBase
{
public:
    using SampleDescriptorBase::SampleDescriptorBase; //to inherit the constructor, copy constructor and assignment of parent class
                                                      // otherwise you should declare them as well as parent class, except the destructor

    std::string name_suffix;
    std::string file_path;
    double cross_section{1};
    std::map<std::string, std::vector<double>> points; // mass for resonant, radion or graviton, nodes for non-resonant
    std::map<std::string, root_ext::Color> draw_ex;
    std::vector<double> norm_sf;

    virtual void CreateWorkingPoints() override
    {
        if(!points.size()) {
            SampleDescriptorBase::CreateWorkingPoints();
            working_points.at(0).file_path = file_path;
            return;
        }
        working_points.clear();
        const size_t N = GetNWorkingPoints();
        for(size_t n = 0; n < N; ++n) {
            Point point;
            point.name = ResolvePattern(name_suffix, n);
            point.full_name = name + "_" + point.name;
            point.title = title.size() ? ResolvePattern(title, n) : point.full_name;
            ReplacePatternItem(point.title, "factor", draw_sf);
            point.file_path = ResolvePattern(file_path, n);
            point.datacard_name = ResolvePattern(datacard_name, n);
            point.norm_sf = n < norm_sf.size() ? norm_sf.at(n) : 1;
            point.draw_sf = draw_sf;
            point.datacard_sf = datacard_sf;
            if(draw_ex.count(point.name)) {
                point.draw = true;
                point.color = draw_ex.at(point.name);
            }
            for(const auto& param_values : points)
                point.param_values.push_back(param_values.second.at(n));
            working_points.push_back(point);
        }
    }

    virtual std::map<std::string, size_t> GetModelParameterNames() const override
    {
        std::map<std::string, size_t> param_names;
        size_t param_id = 0;
        for (auto& p : points)
            param_names[p.first] = param_id++;
        return param_names;
    }

    size_t GetNWorkingPoints() const { return points.size() ? points.begin()->second.size() : 1; }

private:
    std::string ResolvePattern(const std::string& pattern, size_t point_index) const
    {
        if(point_index >= GetNWorkingPoints())
            throw analysis::exception("Signal point chosen is bigger than the size of signal points.");
        std::string result = pattern;
        for (const auto& point_iter : points){
            const auto& point_prefix = point_iter.first;
            const auto& points_list = point_iter.second;
            ReplacePatternItem(result, point_prefix, points_list.at(point_index));
        }
        return result;
    }
};

using SampleDescriptorCollection = std::unordered_map<std::string, SampleDescriptor>;

struct CombinedSampleDescriptor : public SampleDescriptorBase
{
    using SampleDescriptorBase::SampleDescriptorBase;

    std::vector<std::string> sample_descriptors;
};

using CombinedSampleDescriptorCollection = std::unordered_map<std::string, CombinedSampleDescriptor>;

} // namespace analysis

