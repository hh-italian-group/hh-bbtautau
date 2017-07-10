/*! Definition of the file configuration entry reader for TTbar sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "TimeFileDescriptor.h"

namespace analysis {

class TimeFileConfigEntryReader : public analysis::ConfigEntryReader {
public:
    TimeFileConfigEntryReader(TimeFileDescriptorCollection& _descriptors) : descriptors(&_descriptors) {}

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override
    {
        ConfigEntryReader::StartEntry(name, reference_name);
        current = reference_name.size() ? descriptors->at(reference_name) : TimeFileDescriptor();
        current.name = name;
    }

    virtual void EndEntry() override
    {
        CheckReadParamCounts("file_path", 0, Condition::greater_equal);

        (*descriptors)[current.name] = current;
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("file_path", current.file_paths);
    }

private:
    TimeFileDescriptor current;
    TimeFileDescriptorCollection* descriptors;
};

} // namespace analysis
