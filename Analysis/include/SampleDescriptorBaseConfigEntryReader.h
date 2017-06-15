/*! Definition of the file configuration entry reader for analyzer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "SampleDescriptorBase.h"

namespace analysis {


class SampleDescriptorBaseConfigEntryReader : public analysis::ConfigEntryReader {
public:
    SampleDescriptorBaseConfigEntryReader(SampleDescriptorBaseCollection& _descriptors) : descriptors(&_descriptors) {}

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override
    {
        ConfigEntryReader::StartEntry(name, reference_name);
        current = reference_name.size() ? descriptors->at(reference_name) : SampleDescriptorBase();
        current.name = name;
    }

    virtual void EndEntry() override
    {
        CheckReadParamCounts("title", 1, Condition::equal_to);
        CheckReadParamCounts("color", 1, Condition::equal_to);
        CheckReadParamCounts("draw", 1, Condition::equal_to);
        CheckReadParamCounts("channel", 1, Condition::equal_to);
        CheckReadParamCounts("categoryType", 1, Condition::equal_to);
        CheckReadParamCounts("create_hist", 1, Condition::equal_to);
        CheckReadParamCounts("datacard_name", 1, Condition::equal_to);

        (*descriptors)[current.name] = current;
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("title", current.title);
        ParseEntry("color", current.color);
        ParseEntry("draw", current.draw);
        ParseEntry("channel", current.channel);
        ParseEntry("categoryType", current.categoryType);
        ParseEntry("create_hist", current.create_hist);
        ParseEntry("datacard_name", current.datacard_name);
    }

private:
    SampleDescriptorBase current;
    SampleDescriptorBaseCollection* descriptors;
};

} // namespace analysis
