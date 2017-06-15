/*! Definition of the file configuration entry reader for analyzer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "SampleDescriptor.h"
#include "SampleDescriptorBaseConfigEntryReader.h"

namespace analysis {


class SampleDescriptorConfigEntryReader : public analysis::SampleDescriptorBaseConfigEntryReader {
public:
    SampleDescriptorConfigEntryReader(SampleDescriptorCollection& _descriptors) : descriptors(&_descriptors) {}

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override
    {
        ConfigEntryReader::StartEntry(name, reference_name);
        current = reference_name.size() ? descriptors->at(reference_name) : SampleDescriptor();
        current.name = name;
    }

    virtual void EndEntry() override
    {
        CheckReadParamCounts("file_path", 0, Condition::greater_equal);
        CheckReadParamCounts("cross_section", 1, Condition::equal_to);
        CheckReadParamCounts("weight_file", 1, Condition::equal_to);
        CheckReadParamCounts("listSignalPoint", 0, Condition::greater_equal);

        (*descriptors)[current.name] = current;
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("file_path", current.file_paths);
        ParseEntry("cross_section", current.cross_section);
        ParseEntry("weight_file", current.weight_file);
        ParseEntry("listSignalPoint", current.listSignalPoints);
    }

private:
    SampleDescriptor current;
    SampleDescriptorCollection* descriptors;
};

} // namespace analysis
