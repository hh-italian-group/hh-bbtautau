/*! Definition of the analyzer descriptor.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/SampleDescriptor.h"

namespace analysis {

bool AnalyzerSetup::IsSignal(const std::string& sample_name) const
{
    const auto iter = std::find(signals.begin(), signals.end(), sample_name);
    return iter != signals.end();
}

void AnalyzerSetup::CreateLimitSetups()
{
    for(const auto& item: limit_setup_raw) {
        const auto& setup_name = item.first;
        for (size_t i = 0; i < item.second.size(); i++){
            auto variable_categories = SplitValueList(item.second.at(i), false, ":");
            if(variable_categories.size() != 2)
                throw exception("The Number of parameters is %1%, only 2 are allowed") % item.second.size() ;

            const auto categories_str = SplitValueList(variable_categories.at(1), false, ",");
            for (size_t n = 0; n < categories_str.size(); n++){
                const EventCategory categories = ::analysis::Parse<EventCategory>(categories_str.at(n));
                limit_setup[setup_name][categories] = variable_categories.at(0);
            }
        }
    }
}

void AnalyzerSetup::ConvertToEventRegion()
{
    qcd_shape = analysis::Parse<analysis::EventRegion>(qcd_shape_str);

    for(const auto& region : regions_str)
        regions.insert(analysis::Parse<analysis::EventRegion>(region));
}

void MvaReaderSetup::CreateSelections()
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

std::string MvaReaderSetup::GetMvaWPSuffix(const EventSubCategory& sub_category) const
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

MvaReaderSetup MvaReaderSetup::Join(const std::vector<MvaReaderSetup>& setups, const std::string& name)
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

bool SampleDescriptorBase::HasDatacardName() const { return datacard_name.size(); }

void SampleDescriptorBase::CreateWorkingPoints()
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

std::map<std::string, size_t> SampleDescriptorBase::GetModelParameterNames() const { return {}; }


void SampleDescriptor::CreateWorkingPoints()
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

std::map<std::string, size_t> SampleDescriptor::GetModelParameterNames() const
{
    std::map<std::string, size_t> param_names;
    size_t param_id = 0;
    for (auto& p : points)
        param_names[p.first] = param_id++;
    return param_names;
}

size_t SampleDescriptor::GetNWorkingPoints() const { return points.size() ? points.begin()->second.size() : 1; }

std::string SampleDescriptor::ResolvePattern(const std::string& pattern, size_t point_index) const
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

void ModellingUncertainty::CreateSampleUncMap()
{
    samples.clear();
    for(const auto& item : uncertainties) {
        const auto& name = item.first;
        samples[name].unc = uncertainties.Get<ValueType>(name);
        scale_factors.Read(name, samples[name].sf);
    }
}

ModellingUncertaintyCollection::Key::Key(EventCategory _category) : category(_category) {}
ModellingUncertaintyCollection::Key::Key(Channel _channel, EventCategory _category) :
    channel(_channel), category(_category)
{
}

bool ModellingUncertaintyCollection::Key::operator<(const Key& other) const
{
    if(channel != other.channel) return channel < other.channel;
    return category < other.category;
}

std::string ModellingUncertaintyCollection::Key::ToString() const
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

ModellingUncertaintyCollection::Key ModellingUncertaintyCollection::Key::Parse(const std::string& str)
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

void ModellingUncertaintyCollection::Add(const std::string& key_str, const ModellingUncertainty& modelling_unc)
{
    const Key key = Key::Parse(key_str);
    if(unc_map.count(key))
        throw exception("Duplicated modelling unc key = '%1%'.") % key_str;
    if(!key.channel && modelling_unc.ref_category.empty())
        throw exception("Reference category is not specified for '%1%'.") % key_str;
    unc_map[key] = modelling_unc;
}

const ModellingUncertainty& ModellingUncertaintyCollection::Get(Channel channel, EventCategory category) const
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

} // namespace analysis
