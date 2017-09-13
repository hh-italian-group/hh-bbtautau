/*! Definition of the file configuration entry reader for analyzer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "SampleDescriptor.h"

namespace analysis {

class AnalyzerConfigEntryReader : public ConfigEntryReaderT<AnalyzerSetup>, public virtual ConfigEntryReader  {
public:
    using Condition = ConfigEntryReader::Condition;
    using ConfigEntryReaderT<AnalyzerSetup>::ConfigEntryReaderT;

    virtual void EndEntry() override
    {
        CheckReadParamCounts("int_lumi", 1, Condition::less_equal);
        CheckReadParamCounts("final_variable", 0, Condition::greater_equal);
        CheckReadParamCounts("apply_mass_cut", 1, Condition::less_equal);
        CheckReadParamCounts("apply_os_cut", 1, Condition::less_equal);
        CheckReadParamCounts("apply_iso_cut", 1, Condition::less_equal);
        CheckReadParamCounts("energy_scales", 1, Condition::less_equal);
        CheckReadParamCounts("signals", 1, Condition::less_equal);
        CheckReadParamCounts("backgrounds", 1, Condition::less_equal);
        CheckReadParamCounts("data", 1, Condition::less_equal);
        CheckReadParamCounts("data_driven_bkg", 1, Condition::less_equal);
        CheckReadParamCounts("draw_sequence", 1, Condition::less_equal);

        ConfigEntryReaderT<AnalyzerSetup>::EndEntry();
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("int_lumi", current.int_lumi);
        ParseEntry("final_variable", current.final_variables);
        ParseEntry("apply_mass_cut", current.apply_mass_cut);
        ParseEntry("apply_os_cut", current.apply_os_cut);
        ParseEntry("apply_iso_cut", current.apply_iso_cut);
        ParseEnumList("energy_scales", current.energy_scales);
        ParseEntryList("signals", current.signals);
        ParseEntryList("backgrounds", current.backgrounds);
        ParseEntryList("data", current.data);
        ParseEntryList("data_driven_bkg", current.data_driven_bkg);
        ParseEntryList("draw_sequence", current.draw_sequence);
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
        CheckReadParamCounts("draw", 1, Condition::less_equal);
        CheckReadParamCounts("draw_sf", 1, Condition::less_equal);
        CheckReadParamCounts("channels", 1, Condition::less_equal);
        CheckReadParamCounts("categoryType", 1, Condition::less_equal);
        CheckReadParamCounts("datacard_name", 1, Condition::less_equal);

        ConfigEntryReaderT<Descriptor>::EndEntry();
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("title", this->current.title);
        ParseEntry("color", this->current.color);
        ParseEntry("draw", this->current.draw);
        ParseEntry("draw_sf", this->current.draw_sf);
        ParseEntryList("channels", this->current.channels);
        ParseEntry("categoryType", this->current.categoryType);
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
        CheckReadParamCounts("file_path", 0, Condition::greater_equal);
        CheckReadParamCounts("file_path_pattern", 1, Condition::less_equal);
        CheckReadParamCounts("cross_section", 1, Condition::less_equal);
        CheckReadParamCounts("signal_points", 0, Condition::greater_equal);
        CheckReadParamCounts("draw_ex", 0, Condition::greater_equal);
        CheckReadParamCounts("norm_sf", 1, Condition::less_equal);
        CheckReadParamCounts("datacard_name_ex", 0, Condition::greater_equal);

        current.UpdateSignalPoints();

        Base::EndEntry();
    }

    virtual void ReadParameter(const std::string& param_name, const std::string& param_value,
                               std::istringstream& ss) override
    {
        ParseEntry("name_suffix", current.name_suffix);
        ParseEntry("file_path", current.file_paths);
        ParseEntry("file_path_pattern", current.file_path_pattern);
        ParseEntry<double,NumericalExpression>("cross_section", current.cross_section, [](double xs){return xs > 0;});
        ParseEntry("signal_points", current.signal_points_raw);
        ParseEntry("draw_ex", current.draw_ex);
        ParseEntry("norm_sf", current.norm_sf);
        ParseEntry("datacard_name_ex", current.datacard_name_ex);

        Base::ReadParameter(param_name,param_value,ss);
    }

};

class CombineSampleDescriptorConfigEntryReader : public SampleDescriptorBaseConfigEntryReader<CombineSampleDescriptor> {
public:
    using Base = SampleDescriptorBaseConfigEntryReader<CombineSampleDescriptor>;
    using Base::SampleDescriptorBaseConfigEntryReader;

    CombineSampleDescriptorConfigEntryReader(CombineSampleDescriptorCollection& _descriptors,
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

} // namespace analysis





