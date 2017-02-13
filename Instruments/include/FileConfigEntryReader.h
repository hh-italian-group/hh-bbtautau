/*! Definition of the file configuration entry reader.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "FileDescriptor.h"

namespace analysis {

class FileConfigEntryReader : public analysis::ConfigEntryReader {
public:
    FileConfigEntryReader(DYBinDescriptorCollection& _descriptors) : descriptors(&_descriptors) {}

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override
    {
        ConfigEntryReader::StartEntry(name, reference_name);
        current = reference_name.size() ? descriptors->at(reference_name) : DYBinDescriptor();
        current.name = name;
    }

    virtual void EndEntry() override
    {
        CheckReadParamCounts("file_path", 1, Condition::less_equal);
        CheckReadParamCounts("n_jet", 1, Condition::less_equal);
        CheckReadParamCounts("n_bjet", 1, Condition::less_equal);
        CheckReadParamCounts("n_ht", 1, Condition::less_equal);
        CheckReadParamCounts("fileType", 1, Condition::less_equal);

        (*descriptors)[current.name] = current;
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("file_path", current.file_path);
        ParseEntry("n_jet", current.n_jet);
        ParseEntry("n_bjet", current.n_bjet);
        ParseEntry("n_ht", current.n_ht);
        ParseEntry("fileType", current.fileType);
    }

private:
    DYBinDescriptor current;
    DYBinDescriptorCollection* descriptors;
};

} // namespace analysis
