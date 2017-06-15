/*! Definition of the file configuration entry reader for analyzer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "CombineSampleDescriptor.h"
#include "SampleDescriptorConfigEntryReader.h"

namespace analysis {


class CombineSampleDescriptorConfigEntryReader : public analysis::SampleDescriptorConfigEntryReader {
public:
    CombineSampleDescriptorConfigEntryReader(CombineSampleDescriptorCollection& _descriptors) : descriptors(&_descriptors) {}

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override
    {
        ConfigEntryReader::StartEntry(name, reference_name);
        current = reference_name.size() ? descriptors->at(reference_name) : CombineSampleDescriptor();
        current.name = name;
    }

    virtual void EndEntry() override
    {
        CheckReadParamCounts("sample_descriptor", 0, Condition::greater_equal);

        (*descriptors)[current.name] = current;
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("sample_descriptors", current.sample_descriptors);
    }

private:
    CombineSampleDescriptor current;
    CombineSampleDescriptorCollection* descriptors;
};

} // namespace analysis
