/*! Definition of the file configuration entry reader for analyzer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "SampleDescriptor.h"

namespace analysis {

template<typename Descriptor>
class SampleDescriptorBaseConfigEntryReader : public ConfigEntryReaderT<Descriptor>, public virtual ConfigEntryReader  {
public:
    using Condition = ConfigEntryReader::Condition;
    using ConfigEntryReaderT<Descriptor>::ConfigEntryReaderT;

    virtual void EndEntry() override
    {
        CheckReadParamCounts("title", 0, Condition::greater_equal);
        CheckReadParamCounts("color", 0, Condition::greater_equal);
        CheckReadParamCounts("draw", 0, Condition::greater_equal);
        CheckReadParamCounts("channels", 0, Condition::greater_equal);
        CheckReadParamCounts("categoryType", 0, Condition::greater_equal);
        CheckReadParamCounts("datacard_name", 0, Condition::greater_equal);

        ConfigEntryReaderT<Descriptor>::EndEntry();
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("title", this->current.title);
        ParseEntry("color", this->current.color);
        ParseEntry("draw", this->current.draw);
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
        CheckReadParamCounts("file_path", 0, Condition::greater_equal);
        CheckReadParamCounts("file_path_pattern", 1, Condition::less_equal);
        CheckReadParamCounts("cross_section", 0, Condition::greater_equal);
        CheckReadParamCounts("signal_points_raw", 0, Condition::greater_equal);
        CheckReadParamCounts("signal_points", 0, Condition::greater_equal);
        CheckReadParamCounts("draw_ex", 0, Condition::greater_equal);
        CheckReadParamCounts("norm_sf", 0, Condition::greater_equal);
        CheckReadParamCounts("datacard_name_ex", 0, Condition::greater_equal);

        current.UpdateSignalPoints();

        Base::EndEntry();


    }

    virtual void ReadParameter(const std::string& param_name, const std::string& param_value,
                               std::istringstream& ss) override
    {
        ParseEntry("file_path", current.file_paths);
        ParseEntry("file_path_pattern", current.file_path_pattern);
        ParseEntry<double,NumericalExpression>("cross_section", current.cross_section, [](double xs){return xs > 0;});
        ParseEntry("signal_points_raw", current.signal_points_raw);
        ParseEntry("signal_points", current.signal_points);
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
        CheckReadParamCounts("sample_descriptor", 0, Condition::greater_equal);

        Base::EndEntry();

    }

    virtual void ReadParameter(const std::string& param_name, const std::string& param_value,
                               std::istringstream& ss) override
    {
        ParseEntry("sample_descriptor", current.sample_descriptors,
                   [&](const std::string& name){return sampleDescriptorCollection->count(name);});

        Base::ReadParameter(param_name,param_value,ss);
    }

private:
    const SampleDescriptorCollection* sampleDescriptorCollection;
};

} // namespace analysis





