/*! Definition of the file configuration entry reader for TTbar sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "TTFileDescriptor.h"

namespace analysis {

namespace sample_merging{

class TTFileConfigEntryReader : public analysis::ConfigEntryReader {
public:
    TTFileConfigEntryReader(TTBinDescriptorCollection& _descriptors) : descriptors(&_descriptors) {}

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override
    {
        ConfigEntryReader::StartEntry(name, reference_name);
        current = reference_name.size() ? descriptors->at(reference_name) : TTBinDescriptor();
        current.name = name;
    }

    virtual void EndEntry() override
    {
        CheckReadParamCounts("file_path", 0, Condition::greater_equal);
        CheckReadParamCounts("genType", 1, Condition::equal_to);
        CheckReadParamCounts("fileType", 1, Condition::equal_to);

        (*descriptors)[current.name] = current;
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("file_path", current.file_paths);
        ParseEntry("genType", current.genType);
        ParseEntry("fileType", current.fileType);
    }

private:
    TTBinDescriptor current;
    TTBinDescriptorCollection* descriptors;
};

} //namespace sample_merging

} // namespace analysis
