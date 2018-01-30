/*! Definition of the file configuration entry reader for analyzer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "SampleDescriptor.h"

namespace analysis {

class AnalyzerConfigEntryReader : public ConfigEntryReaderT<AnalyzerSetup> {
public:
    using Condition = ConfigEntryReader::Condition;
    using ConfigEntryReaderT<AnalyzerSetup>::ConfigEntryReaderT;

    virtual void EndEntry() override
    {
        CheckReadParamCounts("int_lumi", 1, Condition::less_equal);
        CheckReadParamCounts("period", 1, Condition::less_equal);
        CheckReadParamCounts("final_variables", 1, Condition::less_equal);
        CheckReadParamCounts("apply_mass_cut", 1, Condition::less_equal);
        CheckReadParamCounts("apply_os_cut", 1, Condition::less_equal);
        CheckReadParamCounts("apply_iso_cut", 1, Condition::less_equal);
        CheckReadParamCounts("energy_scales", 1, Condition::less_equal);
        CheckReadParamCounts("categories", 1, Condition::less_equal);
        CheckReadParamCounts("sub_categories", 1, Condition::less_equal);
        CheckReadParamCounts("regions", 1, Condition::less_equal);
        CheckReadParamCounts("data", 1, Condition::less_equal);
        CheckReadParamCounts("signals", 1, Condition::less_equal);
        CheckReadParamCounts("backgrounds", 1, Condition::less_equal);
        CheckReadParamCounts("cmb_samples", 1, Condition::less_equal);
        CheckReadParamCounts("draw_sequence", 1, Condition::less_equal);
        CheckReadParamCounts("limit_category", 0, Condition::greater_equal);
        CheckReadParamCounts("mva_setup", 1, Condition::less_equal);
        CheckReadParamCounts("hist_cfg", 1, Condition::less_equal);
        CheckReadParamCounts("syncDataIds", 1, Condition::less_equal);
        CheckReadParamCounts("plot_cfg", 1, Condition::less_equal);
        CheckReadParamCounts("plot_page_opt", 1, Condition::less_equal);
        CheckReadParamCounts("massWindowParams", 0, Condition::greater_equal);
        CheckReadParamCounts("unc_cfg", 1, Condition::less_equal);

        ConfigEntryReaderT<AnalyzerSetup>::EndEntry();
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("int_lumi", current.int_lumi);
        ParseEntry("period", current.period);
        ParseEntryList("final_variables", current.final_variables);
        ParseEntry("apply_mass_cut", current.apply_mass_cut);
        ParseEntry("apply_os_cut", current.apply_os_cut);
        ParseEntry("apply_iso_cut", current.apply_iso_cut);
        ParseEnumList("energy_scales", current.energy_scales);
        ParseEnumList("categories", current.categories);
        ParseEnumList("sub_categories", current.sub_categories);
        ParseEnumList("regions", current.regions);
        ParseEntryList("data", current.data);
        ParseEntryList("signals", current.signals);
        ParseEntryList("backgrounds", current.backgrounds);
        ParseEntryList("cmb_samples", current.cmb_samples);
        ParseEntryList("draw_sequence", current.draw_sequence);
        ParseEntry("limit_category", current.limit_categories);
        ParseEntry("mva_setup", current.mva_setup);
        ParseEntry("hist_cfg", current.hist_cfg);
        ParseEntryList("syncDataIds", current.syncDataIds);
        ParseEntry("plot_cfg", current.plot_cfg);
        ParseEntry("plot_page_opt", current.plot_page_opt);
        ParseEntry("massWindowParams", current.massWindowParams);
        ParseEntry("unc_cfg", current.unc_cfg);
    }
};

class MvaReaderSetupEntryReader : public ConfigEntryReaderT<MvaReaderSetup>  {
public:
    using Condition = ConfigEntryReader::Condition;
    using ConfigEntryReaderT<MvaReaderSetup>::ConfigEntryReaderT;

    virtual void EndEntry() override
    {
        CheckReadParamCounts("training", 0, Condition::greater_equal);
        CheckReadParamCounts("variables", 0, Condition::greater_equal);
        CheckReadParamCounts("masses", 0, Condition::greater_equal);
        CheckReadParamCounts("spins", 0, Condition::greater_equal);
        CheckReadParamCounts("cuts", 0, Condition::greater_equal);
        CheckReadParamCounts("legacy", 0, Condition::greater_equal);
        CheckReadParamCounts("training_range", 0, Condition::greater_equal);
        CheckReadParamCounts("samples", 0, Condition::greater_equal);

        current.CreateSelections();
        ConfigEntryReaderT<MvaReaderSetup>::EndEntry();
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("training", current.trainings);
        ParseMappedEntryList("variables", current.variables, false);
        ParseMappedEntryList("masses", current.masses, true);
        ParseMappedEntryList("spins", current.spins, true);
        ParseMappedEntryList("cuts", current.cuts, false);
        ParseEntry("training_range", current.training_ranges);
        ParseMappedEntryList("samples", current.samples, false);
        ParseEntry("legacy", current.legacy);
    }
};

template<typename Descriptor>
class SampleDescriptorBaseConfigEntryReader : public ConfigEntryReaderT<Descriptor>, public virtual ConfigEntryReader {
public:
    using Condition = ConfigEntryReader::Condition;
    using ConfigEntryReaderT<Descriptor>::ConfigEntryReaderT;

    virtual void EndEntry() override
    {
        CheckReadParamCounts("title", 1, Condition::less_equal);
        CheckReadParamCounts("color", 1, Condition::less_equal);
        CheckReadParamCounts("draw_sf", 1, Condition::less_equal);
        CheckReadParamCounts("channels", 1, Condition::less_equal);
        CheckReadParamCounts("sample_type", 1, Condition::less_equal);
        CheckReadParamCounts("datacard_name", 1, Condition::less_equal);

        this->current.CreateWorkingPoints();
        ConfigEntryReaderT<Descriptor>::EndEntry();
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("title", this->current.title);
        ParseEntry("color", this->current.color);
        ParseEntry<double, NumericalExpression>("datacard_sf", this->current.datacard_sf,
                                                [](double sf){ return sf > 0; });
        ParseEntry<double, NumericalExpression>("draw_sf", this->current.draw_sf, [](double sf){ return sf > 0; });
        ParseEntryList("channels", this->current.channels);
        ParseEntry("sample_type", this->current.sampleType);
        ParseEntry("datacard_name", this->current.datacard_name);
    }
};

class SampleDescriptorConfigEntryReader : public SampleDescriptorBaseConfigEntryReader<SampleDescriptor> {
public:
    using Base = SampleDescriptorBaseConfigEntryReader<SampleDescriptor>;
    using Base::SampleDescriptorBaseConfigEntryReader;


    virtual void EndEntry() override
    {
        CheckReadParamCounts("name_suffix", 1, Condition::less_equal);
        CheckReadParamCounts("file_path", 1, Condition::less_equal);
        CheckReadParamCounts("cross_section", 1, Condition::less_equal);
        CheckReadParamCounts("points", 0, Condition::greater_equal);
        CheckReadParamCounts("draw_ex", 0, Condition::greater_equal);
        CheckReadParamCounts("norm_sf", 1, Condition::less_equal);

        Base::EndEntry();
    }

    virtual void ReadParameter(const std::string& param_name, const std::string& param_value,
                               std::istringstream& ss) override
    {
        ParseEntry("name_suffix", current.name_suffix);
        ParseEntry("file_path", current.file_path);
        ParseEntry<double, NumericalExpression>("cross_section", current.cross_section,
                                                [](double xs){ return xs > 0; });
        ParseMappedEntryList("points", current.points, true);
        ParseEntry("draw_ex", current.draw_ex);
        ParseEntryList("norm_sf", current.norm_sf, true);

        Base::ReadParameter(param_name,param_value,ss);
    }

};

class CombinedSampleDescriptorConfigEntryReader :
        public SampleDescriptorBaseConfigEntryReader<CombinedSampleDescriptor> {
public:
    using Base = SampleDescriptorBaseConfigEntryReader<CombinedSampleDescriptor>;
    using Base::SampleDescriptorBaseConfigEntryReader;

    CombinedSampleDescriptorConfigEntryReader(CombinedSampleDescriptorCollection& _descriptors,
                                             const SampleDescriptorCollection& _sampleDescriptors) :
       Base(_descriptors), sampleDescriptorCollection(&_sampleDescriptors) {}

    virtual void EndEntry() override
    {
        CheckReadParamCounts("sample_descriptors", 1, Condition::equal_to);
        Base::EndEntry();
    }

    virtual void ReadParameter(const std::string& param_name, const std::string& param_value,
                               std::istringstream& ss) override
    {
        ParseEntryList("sample_descriptors", current.sample_descriptors, false, " \t",
                   [&](const std::string& name){return sampleDescriptorCollection->count(name);});

        Base::ReadParameter(param_name,param_value,ss);
    }

private:
    const SampleDescriptorCollection* sampleDescriptorCollection;
};

class ModellingUncertaintyEntryReader : public ConfigEntryReader {
public:
    ModellingUncertaintyEntryReader(ModellingUncertaintyCollection& _items) : items(&_items) {}

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override
    {
        ConfigEntryReader::StartEntry(name, reference_name);
        current = ModellingUncertainty();
        current_name = name;
    }

    virtual void EndEntry() override
    {
        CheckReadParamCounts("unc", 1, Condition::less_equal);
        CheckReadParamCounts("sf", 1, Condition::less_equal);
        CheckReadParamCounts("ref_category", 1, Condition::less_equal);

        current.CreateSampleUncMap();
        items->Add(current_name, current);
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("unc", current.uncertainties);
        ParseEntry("sf", current.scale_factors);
        ParseEntry("ref_category", current.ref_category);
    }

private:
    std::string current_name;
    ModellingUncertainty current;
    ModellingUncertaintyCollection* items;
};

} // namespace analysis





