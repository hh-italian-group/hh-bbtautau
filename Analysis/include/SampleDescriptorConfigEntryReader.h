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

    virtual void EndEntry() override;
    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override;
};

class MvaReaderSetupEntryReader : public ConfigEntryReaderT<MvaReaderSetup>  {
public:
    using Condition = ConfigEntryReader::Condition;
    using ConfigEntryReaderT<MvaReaderSetup>::ConfigEntryReaderT;

    virtual void EndEntry() override;
    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override;
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
        CheckReadParamCounts("postfit_name", 1, Condition::less_equal);
        CheckReadParamCounts("norm_sf_file", 1, Condition::less_equal);
        CheckReadParamCounts("NLO_weight_file", 1, Condition::less_equal);
        CheckReadParamCounts("fit_method", 1, Condition::less_equal);
        CheckReadParamCounts("sample_order", 1, Condition::less_equal);

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
        ParseEntry("postfit_name", this->current.postfit_name);
        ParseEntry("norm_sf_file", this->current.norm_sf_file);
        ParseEntry("NLO_weight_file", this->current.NLO_weight_file);
        ParseEntry("fit_method", this->current.fit_method);
        ParseEntry("sample_order", this->current.sampleOrder);
    }
};

class SampleDescriptorConfigEntryReader : public SampleDescriptorBaseConfigEntryReader<SampleDescriptor> {
public:
    using Base = SampleDescriptorBaseConfigEntryReader<SampleDescriptor>;
    using Base::SampleDescriptorBaseConfigEntryReader;

    virtual void EndEntry() override;
    virtual void ReadParameter(const std::string& param_name, const std::string& param_value,
                               std::istringstream& ss) override;
};

class CombinedSampleDescriptorConfigEntryReader :
        public SampleDescriptorBaseConfigEntryReader<CombinedSampleDescriptor> {
public:
    using Base = SampleDescriptorBaseConfigEntryReader<CombinedSampleDescriptor>;
    using Base::SampleDescriptorBaseConfigEntryReader;

    CombinedSampleDescriptorConfigEntryReader(CombinedSampleDescriptorCollection& _descriptors,
                                             const SampleDescriptorCollection& _sampleDescriptors);

    virtual void EndEntry() override;
    virtual void ReadParameter(const std::string& param_name, const std::string& param_value,
                               std::istringstream& ss) override;

private:
    const SampleDescriptorCollection* sampleDescriptorCollection;
};

class ModellingUncertaintyEntryReader : public ConfigEntryReader {
public:
    ModellingUncertaintyEntryReader(ModellingUncertaintyCollection& _items);
    virtual void StartEntry(const std::string& name, const std::string& reference_name) override;
    virtual void EndEntry() override;
    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override;

private:
    std::string current_name;
    ModellingUncertainty current;
    ModellingUncertaintyCollection* items;
};

} // namespace analysis
