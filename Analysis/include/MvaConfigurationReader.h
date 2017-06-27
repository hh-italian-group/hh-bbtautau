/*! Definition of MvaVariablesStudy class, the main class for Mva studies. 
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "hh-bbtautau/Analysis/include/MvaConfiguration.h"

namespace  analysis {
namespace mva_study{
class MvaConfigReader : public analysis::ConfigEntryReader {
public:
    MvaConfigReader(MvaSetupCollection& _descriptors) : descriptors(&_descriptors) {}

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override
    {
        ConfigEntryReader::StartEntry(name, reference_name);
        current = reference_name.size() ? descriptors->at(reference_name) : MvaSetup();
        current.name = name;
    }

    virtual void EndEntry() override
    {
        CheckReadParamCounts("channels", 0, Condition::greater_equal);
        CheckReadParamCounts("mass_range", 0, Condition::greater_equal);
        CheckReadParamCounts("param_list", 0, Condition::greater_equal);
        CheckReadParamCounts("param", 0, Condition::greater_equal);
        CheckReadParamCounts("param_range", 0, Condition::greater_equal);
        CheckReadParamCounts("disabled_params", 0, Condition::greater_equal);
        CheckReadParamCounts("significant_params", 0, Condition::greater_equal);
        CheckReadParamCounts("variables", 0, Condition::greater_equal);
        CheckReadParamCounts("use_mass_var", 0, Condition::greater_equal);

        (*descriptors)[current.name] = current;
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntryList("channels", current.channels, false);
        ParseEntry("mass_range", current.mass_range);
        ParseEntry("param_list", current.param_list);
        ParseEntry("param", current.params);
        ParseEntry("param_range", current.param_range);
        ParseEntry("disabled_params", current.disabled_params);
        ParseEntryList("significant_params", current.significant_params, false);
        ParseEntryList("variables", current.variables, false);
        ParseEntry("use_mass_var", current.use_mass_var);
    }

private:
    MvaSetup current;
    MvaSetupCollection* descriptors;
};

class SampleConfigReader : public analysis::ConfigEntryReader {
public:
    SampleConfigReader(SampleEntryListCollection& _descriptors) : descriptors(&_descriptors) {}

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override
    {
        ConfigEntryReader::StartEntry(name, reference_name);
        current = reference_name.size() ? descriptors->at(reference_name) : SampleEntryList();
        current.name = name;
    }

    virtual void EndEntry() override
    {
        CheckReadParamCounts("file", 1, Condition::greater_equal);

        (*descriptors)[current.name] = current;
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("file", current.files);
    }

private:
    SampleEntryList current;
    SampleEntryListCollection* descriptors;
};

}
}
