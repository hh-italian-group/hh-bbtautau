/*! Definition of the file configuration entry reader for analyzer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalyzerDescriptor.h"

namespace analysis {


class AnalyzerConfigEntryReader : public analysis::ConfigEntryReader {
public:
    AnalyzerConfigEntryReader(AnalyzerDescriptorCollection& _descriptors) : descriptors(&_descriptors) {}

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override
    {
        ConfigEntryReader::StartEntry(name, reference_name);
        current = reference_name.size() ? descriptors->at(reference_name) : AnalyzerDescriptor();
        current.name = name;
    }

    virtual void EndEntry() override
    {
        CheckReadParamCounts("file_cross_section_map", 0, Condition::greater_equal);
        CheckReadParamCounts("categoryType", 1, Condition::equal_to);
//        CheckReadParamCounts("cross_section", 1, Condition::equal_to);
        CheckReadParamCounts("int_lumi", 1, Condition::equal_to);
        CheckReadParamCounts("channel", 1, Condition::equal_to);
        CheckReadParamCounts("color", 1, Condition::less_equal);

        (*descriptors)[current.name] = current;
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("file_cross_section_map", current.file_cross_section_map);
        ParseEntry("fileType", current.categoryType);
//        ParseEntry("cross_section", current.cross_section);
        ParseEntry("int_lumi", current.int_lumi);
        ParseEntry("channel", current.channel);
        ParseEntry("color", current.color);
    }

private:
    AnalyzerDescriptor current;
    AnalyzerDescriptorCollection* descriptors;
};

} // namespace analysis
