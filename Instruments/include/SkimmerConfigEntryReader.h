/*! Definition of classes to read TupleSkimmer configuration.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "SkimmerConfig.h"

namespace analysis {
namespace tuple_skimmer {

class SetupEntryReader : public ConfigEntryReaderT<Setup> {
public:
    using ConfigEntryReaderT<Setup>::ConfigEntryReaderT;

    virtual void EndEntry() override
    {
        CheckReadParamCounts("energy_scales", 1, Condition::equal_to);
        CheckReadParamCounts("channels", 1, Condition::equal_to);
        CheckReadParamCounts("tau_ids", 1, Condition::equal_to);

        ConfigEntryReaderT<Setup>::EndEntry();
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEnumList("energy_scales", current.energy_scales);
        ParseEnumList("channels", current.channels);
        ParseEntryList("tau_ids", current.tau_ids);
    }
};

class SkimJobEntryReader : public ConfigEntryReaderT<SkimJob> {
public:
    using ConfigEntryReaderT<SkimJob>::ConfigEntryReaderT;

    virtual void EndEntry() override
    {
        CheckReadParamCounts("energy_scales", 1, Condition::equal_to);
        CheckReadParamCounts("channels", 1, Condition::equal_to);
        CheckReadParamCounts("tau_ids", 1, Condition::equal_to);

        ConfigEntryReaderT<SkimJob>::EndEntry();
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEnumList("dy_weight", current.dy_weight);
        ParseEnumList("sm_weight", current.sm_weight);
        ParseEntryList("tt_weight", current.tt_weight);
        ParseEntryList("wjets_weight", current.wjets_weight);
        ParseEntryList("merged_output", current.merged_output);
    }

private:
    void ParseFileDescriptor(const std::string& param_name, const std::string& param_value)
    {
        std::vector<std::string> inputs;
        if(param_name == "file") {
            inputs = SplitValueList(param_value, false);
            current.files.emplace_back(inputs);
        } else if(param_name == "file_ex") {
            auto columns = SplitValueList(param_value, false);
            if(columns.size() < 2)
                throw exception("Invalid extended file description.");
            inputs.insert(inputs.end(), columns.begin() + 1, columns.end());
            current.files.emplace_back(inputs, columns.at(0));
        } else if(param_name == "file_xs") {
            auto columns = SplitValueList(param_value, false);
            if(columns.size() < 2)
                throw exception("Invalid description for file with cross-section.");
            inputs.insert(inputs.end(), columns.begin() + 1, columns.end());
            const double xs = Parse<double>(columns.at(0));
            current.files.emplace_back(inputs, xs);
        } else {
            return;
        }
        throw param_parsed_exception();
    }
};

} // namespace tuple_skimmer
} // namespace analysis
