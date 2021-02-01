/*! Definition of the analyzer descriptor.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/PropertyConfigReader.h"
#include "AnalysisTools/Print/include/PlotPrimitives.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Analysis/include/SignalObjectSelector.h"
#include "h-tautau/JetTools/include/BTagger.h"
#include "AnalysisCategories.h"
#include "EventAnalyzerDataId.h"

namespace analysis {

struct AnalyzerSetup {
    std::string name;
    double int_lumi{0};
    Period period;
    SignalMode mode;
    QCDmethod qcd_method;
    std::string qcd_shape_str;
    std::string xs_cfg;
    EventRegion qcd_shape;
    DiscriminatorWP tauID_wp;
    std::vector<double> pt_sel_bins;
    bool use_kinFit{false}, use_svFit{false}, allow_calc_svFit{false}, use_IterativeFit{false};
    std::map<std::string, std::set<UncertaintySource>> unc_sources;
    std::set<UncertaintySource> norm_unc_sources;
    EventCategorySet categories;
    EventSubCategorySet sub_categories;
    EventCategorySet categories_base;
    EventSubCategorySet sub_categories_base;
    EventRegionSet regions;
    std::vector<std::string> regions_str;
    std::vector<std::string> data, signals, backgrounds, cmb_samples;
    std::vector<std::string> draw_sequence;
    std::map<EventCategory, std::string> limit_categories;
    std::string mva_setup, trigger_path;
    std::vector<std::string> syncDataIds;
    std::vector<std::string> hist_cfg;
    std::string plot_cfg, plot_page_opt, unc_cfg;
    std::map<Channel, std::string> r_factors_file;
    BTaggerKind jet_ordering;
    double qcd_ss_os_sf{0};
    double qcd_ss_os_err{0};
    std::string mdnn_version;
    std::map<Channel, std::vector<std::string>> trigger;
    std::map<Channel, std::vector<std::string>> trigger_vbf;
    std::map<SelectionCut,analysis::EllipseParameters> massWindowParams;
    std::map<std::string, std::vector<std::string>> limit_setup_raw;
    std::map<std::string, std::map<EventCategory, std::string>> limit_setup;

    bool IsSignal(const std::string& sample_name) const;
    void CreateLimitSetups();
    void ConvertToEventRegion();
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

    void CreateSelections();
    std::string GetMvaWPSuffix(const EventSubCategory& sub_category) const;

    static MvaReaderSetup Join(const std::vector<MvaReaderSetup>& setups, const std::string& name = "");

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
        std::string name, full_name, title, file_path, datacard_name, trigger, trigger_vbf;
        SampleType sampleType;
        double norm_sf{1}, datacard_sf{1}, draw_sf{1};
        std::string cross_section;
        bool draw{false};
        root_ext::Color color{kBlack};
        std::vector<std::string> param_values;
    };
    using PointCollection = std::vector<Point>;

    std::string name;
    std::string title;
    root_ext::Color color{kBlack};
    double datacard_sf{1}, draw_sf{1};
    std::set<Channel> channels;
    SampleType sampleType{SampleType::MC};
    std::string datacard_name;
    std::string postfit_name;
    std::string norm_sf_file;
    std::string NLO_weight_file;
    std::string sampleOrder{"LO"};
    DYFitModel fit_method{DYFitModel::None};
    bool apply_top_pt_unc{false}, is_TuneCP5{false};  /* is_TuneCP5 only valid for 2016 */ \

    PointCollection working_points;

    SampleDescriptorBase() {} //constructor
    SampleDescriptorBase(const SampleDescriptorBase& ) = default; //copy constructor
    virtual ~SampleDescriptorBase(){} //destructor

    SampleDescriptorBase& operator= ( const SampleDescriptorBase& ) = default; //assignment

    bool HasDatacardName() const;
    virtual void CreateWorkingPoints();
    virtual std::map<std::string, size_t> GetModelParameterNames() const;

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
    using SampleDescriptorBase::SampleDescriptorBase; //to inherit the constructor, copy constructor and assignment of
                                                      // parent class, otherwise you should declare them as well as
                                                      // parent class, except the destructor

    std::string name_suffix, reference_pu_sample;
    std::string file_path;
    double cross_section{1};
    std::map<std::string, std::vector<std::string>> points; // mass for resonant, radion or graviton, nodes for non-resonant
    std::map<std::string, root_ext::Color> draw_ex;
    std::vector<double> norm_sf;
    std::vector<std::string> point_xs;
    bool create_orthogonal_points{false};

    virtual void CreateWorkingPoints() override;
    virtual std::map<std::string, size_t> GetModelParameterNames() const override;
    size_t GetNWorkingPoints() const;

private:
    std::string ResolvePattern(const std::string& pattern, size_t point_index) const;
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

    void CreateSampleUncMap();
};

class ModellingUncertaintyCollection {
private:
    struct Key {
        static constexpr char wildcard = '*', separator = '/';
        boost::optional<Channel> channel;
        EventCategory category;

        Key() {}
        Key(EventCategory _category);
        Key(Channel _channel, EventCategory _category);
        bool operator<(const Key& other) const;
        std::string ToString() const;

        static Key Parse(const std::string& str);
    };

    using UncMap = std::map<Key, ModellingUncertainty>;

public:
    void Add(const std::string& key_str, const ModellingUncertainty& modelling_unc);
    const ModellingUncertainty& Get(Channel channel, EventCategory category) const;

private:
    UncMap unc_map;
};

} // namespace analysis
