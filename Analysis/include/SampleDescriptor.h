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
#include "AnalysisTools/Core/include/PropertyConfigReader.h"
#include <boost/algorithm/string.hpp>
#include "AnalysisCategories.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "EventAnalyzerDataId.h"

namespace analysis {

struct AnalyzerSetup {
    std::string name;
    double int_lumi{0};
    Period period;
    std::vector<std::string> final_variables;
    bool apply_mass_cut{false}, apply_os_cut{true}, apply_iso_cut{true};
    EventEnergyScaleSet energy_scales;
    EventCategorySet categories;
    EventSubCategorySet sub_categories;
    EventRegionSet regions;
    std::vector<std::string> data, signals, backgrounds, cmb_samples;
    std::vector<std::string> draw_sequence;
    std::map<EventCategory, std::string> limit_categories;
    std::string mva_setup, hist_cfg;
    std::vector<EventAnalyzerDataId> syncDataIds;
    std::string plot_cfg, plot_page_opt, unc_cfg;

    std::map<SelectionCut,analysis::EllipseParameters> massWindowParams;

    bool IsSignal(const std::string& sample_name) const
    {
        const auto iter = std::find(signals.begin(), signals.end(), sample_name);
        return iter != signals.end();
    }
};

using AnalyzerSetupCollection = std::unordered_map<std::string, AnalyzerSetup>;

struct MvaReaderSetup {
    using Range = ::analysis::Range<unsigned>;
    using StrSet = std::unordered_set<std::string>;

    struct Params {
        std::string name;
        int spin; double mass;
        double cut;
        boost::optional<Range> training_range;
        StrSet samples;
    };

    std::string name;
    std::map<std::string, std::string> trainings;
    std::map<std::string, StrSet> variables;
    std::map<std::string, std::vector<double>> masses;
    std::map<std::string, std::vector<int>> spins;
    std::map<std::string, std::vector<double>> cuts;
    std::map<std::string, std::string> legacy;
    std::map<std::string, Range> training_ranges;
    std::map<std::string, StrSet> samples;

    std::map<SelectionCut, Params> selections;

    void CreateSelections()
    {
        selections.clear();

        const size_t n_min = static_cast<size_t>(SelectionCut::MVA_first);
        const size_t n_max = static_cast<size_t>(SelectionCut::MVA_last);
        size_t n = n_min;
        for(const auto& training : trainings) {
            const std::string tr_name = training.first;
            Params params;
            params.name = tr_name;
            if(!spins.count(tr_name))
                throw exception("Missing spin working points for MVA training %1%.") % tr_name;
            const auto& tr_spins = spins.at(tr_name);
            if(!masses.count(tr_name))
                throw exception("Missing mass working points for MVA training %1%.") % tr_name;
            const auto& tr_masses = masses.at(tr_name);
            if(tr_spins.size() != tr_masses.size())
                throw exception("Inconsistend number of spin and mass working points for MVA training %1%.") % tr_name;
            if(!cuts.count(tr_name))
                throw exception("Missing cut working points for MVA training %1%.") % tr_name;
            const auto& tr_cuts = cuts.at(tr_name);
            if(training_ranges.count(tr_name)) {
                params.training_range = training_ranges.at(tr_name);
                if(!samples.count(tr_name))
                    throw exception("Used samples not specified for MVA training %1%.") % tr_name;
                params.samples = samples.at(tr_name);
            }

            for(size_t wp_index = 0; wp_index < tr_spins.size(); ++wp_index) {
                params.spin = tr_spins.at(wp_index);
                params.mass = tr_masses.at(wp_index);
                for(double cut_wp : tr_cuts) {
                    params.cut = cut_wp;
                    if(n > n_max)
                        throw exception("Exceded the maximum number (%1%) of the available MVA selection points")
                            % (n_max - n_min + 1);
                    SelectionCut sel = static_cast<SelectionCut>(n);
                    selections[sel] = params;
                    ++n;
                }
            }
        }
    }

    std::string GetMvaWPSuffix(const EventSubCategory& sub_category) const
    {
        std::ostringstream ss;
        for(const auto& sel : selections) {
            if(!sub_category.HasCut(sel.first)) continue;
            const Params& params = sel.second;
            ss << params.name << "_spin" << params.spin << "_M" << params.mass << "_cut" << params.cut << "_";
        }
        std::string suffix = ss.str();
        if(suffix.size())
            suffix.erase(suffix.size() - 1);
        return suffix;
    }

    static MvaReaderSetup Join(const std::vector<MvaReaderSetup>& setups, const std::string& name = "")
    {
        MvaReaderSetup joined;
        joined.name = name;
        for(const auto& setup : setups) {
            if(!setup.trainings.size()) continue;
            const bool need_suffix = setup.trainings.size() > 1;
            InsertMap(joined.trainings, setup.trainings, setup.name, need_suffix);
            InsertMap(joined.variables, setup.variables, setup.name, need_suffix);
            InsertMap(joined.masses, setup.masses, setup.name, need_suffix);
            InsertMap(joined.spins, setup.spins, setup.name, need_suffix);
            InsertMap(joined.cuts, setup.cuts, setup.name, need_suffix);
            InsertMap(joined.legacy, setup.legacy, setup.name, need_suffix);
            InsertMap(joined.training_ranges, setup.training_ranges, setup.name, need_suffix);
            InsertMap(joined.samples, setup.samples, setup.name, need_suffix);
        }
        joined.CreateSelections();
        return joined;
    }

private:
    template<typename Value>
    static void InsertMap(std::map<std::string, Value>& target, const std::map<std::string, Value>& original,
                          const std::string& setup_name, bool need_suffix)
    {
        for(const auto& item : original) {
            std::string item_name = setup_name;
            if(need_suffix)
                item_name += "_" + item.first;
            if(target.count(item_name))
                throw exception("Unable to merge MvaReader setups. Duplicated entries are found.");
            target[item_name] = item.second;
        }
    }
};

using MvaReaderSetupCollection = std::unordered_map<std::string, MvaReaderSetup>;

struct SampleDescriptorBase {
    struct Point {
        std::string name, full_name, title, file_path, datacard_name;
        SampleType sampleType;
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
        point.sampleType = sampleType;
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
            point.sampleType = sampleType;
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
            throw exception("Signal point chosen is bigger than the size of signal points.");
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

struct ModellingUncertainty {
    using ValueType = double;
    struct SampleUnc { ValueType unc{0}, sf{1}; };
    using SampleUncMap = std::map<std::string, SampleUnc>;

    PropertyList uncertainties, scale_factors;
    std::string ref_category;
    SampleUncMap samples;

    void CreateSampleUncMap()
    {
        samples.clear();
        for(const auto& item : uncertainties) {
            const auto& name = item.first;
            samples[name].unc = uncertainties.Get<ValueType>(name);
            scale_factors.Read(name, samples[name].sf);
        }
    }
};

class ModellingUncertaintyCollection {
private:
    struct Key {
        static constexpr char wildcard = '*', separator = '/';
        boost::optional<Channel> channel;
        EventCategory category;

        Key() {}
        Key(EventCategory _category) : category(_category) {}
        Key(Channel _channel, EventCategory _category) : channel(_channel), category(_category) {}

        bool operator<(const Key& other) const
        {
            if(channel != other.channel) return channel < other.channel;
            return category < other.category;
        }

        std::string ToString() const
        {
            std::ostringstream ss;
            if(channel)
                ss << *channel;
            else
                ss << wildcard;
            ss << separator;
            ss << category;
            return ss.str();
        }

        static Key Parse(const std::string& str)
        {
            static const std::string separators(1, separator);
            const auto items = SplitValueList(str, true, separators, false);
            if(items.size() != 2)
                throw exception("Invalid modelling uncertainty key = '%1%'.") % str;
            Key key;
            if(items.at(0).size() != 1 || items.at(0).at(0) != wildcard)
                key.channel = ::analysis::Parse<Channel>(items.at(0));
            key.category = ::analysis::Parse<EventCategory>(items.at(1));
            return key;
        }
    };

    using UncMap = std::map<Key, ModellingUncertainty>;

public:
    void Add(const std::string& key_str, const ModellingUncertainty& modelling_unc)
    {
        const Key key = Key::Parse(key_str);
        if(unc_map.count(key))
            throw exception("Duplicated modelling unc key = '%1%'.") % key_str;
        if(!key.channel && modelling_unc.ref_category.empty())
            throw exception("Reference category is not specified for '%1%'.") % key_str;
        unc_map[key] = modelling_unc;
    }

    const ModellingUncertainty& Get(Channel channel, EventCategory category) const
    {
        const Key key(channel, category);
        if(unc_map.count(key))
            return unc_map.at(key);
        const Key meta_key(category);
        if(unc_map.count(meta_key)) {
            const auto ref_category = Parse<EventCategory>(unc_map.at(meta_key).ref_category);
            const Key ref_key(channel, ref_category);
            if(unc_map.count(ref_key))
                return unc_map.at(ref_key);
        }

        throw exception("Modelling uncertainties not found for '%1%'.") % key.ToString();
    }

private:
    UncMap unc_map;
};

} // namespace analysis

