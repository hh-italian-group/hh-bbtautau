/*! Definition of MvaVariablesStudy class, the main class for Mva studies. 
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "hh-bbtautau/Analysis/include/MvaConfiguration.h"

namespace  analysis {
namespace mva_study{
class MvaConfigReader : public ConfigEntryReaderT<MvaSetup> {
public:
    using ConfigEntryReaderT<MvaSetup>::ConfigEntryReaderT;

    virtual void EndEntry() override
    {
        CheckReadParamCounts("channels", 1, Condition::less_equal);
        CheckReadParamCounts("mass_range", 1, Condition::less_equal);
        CheckReadParamCounts("param_list", 0, Condition::greater_equal);
        CheckReadParamCounts("param", 0, Condition::greater_equal);
        CheckReadParamCounts("param_range", 0, Condition::greater_equal);
        CheckReadParamCounts("disabled_params", 0, Condition::greater_equal);
        CheckReadParamCounts("significant_params", 0, Condition::greater_equal);
        CheckReadParamCounts("variables", 1, Condition::less_equal);
        CheckReadParamCounts("use_mass_var", 1, Condition::less_equal);

        ConfigEntryReaderT<MvaSetup>::EndEntry();
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
};

class SampleConfigReader : public analysis::ConfigEntryReaderT<SampleEntryList> {
public:
    using ConfigEntryReaderT<SampleEntryList>::ConfigEntryReaderT;

    virtual void EndEntry() override
    {
        CheckReadParamCounts("file", 1, Condition::greater_equal);

        ConfigEntryReaderT<SampleEntryList>::EndEntry();
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("file", current.files);
    }

};

}
}
