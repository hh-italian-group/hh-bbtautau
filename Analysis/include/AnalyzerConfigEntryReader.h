/*! Definition of the file configuration entry reader for analyzer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalyzerSetup.h"

namespace analysis {


class AnalyzerConfigEntryReader : public analysis::ConfigEntryReader {
public:
    AnalyzerConfigEntryReader(AnalyzerSetupCollection& _descriptors) : descriptors(&_descriptors) {}

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override
    {
        ConfigEntryReader::StartEntry(name, reference_name);
        current = reference_name.size() ? descriptors->at(reference_name) : AnalyzerSetup();
        current.name = name;
    }

    virtual void EndEntry() override
    {
        CheckReadParamCounts("int_lumi", 1, Condition::equal_to);
        CheckReadParamCounts("final_variable", 0, Condition::greater_equal);
        CheckReadParamCounts("apply_mass_cut", 1, Condition::equal_to);
        CheckReadParamCounts("energy_scales", 0, Condition::greater_equal);

        (*descriptors)[current.name] = current;
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("int_lumi", current.int_lumi);
        ParseEntry("final_variable", current.final_variables);
        ParseEntry("apply_mass_cut", current.apply_mass_cut);
        ParseEnumList("energy_scales", current.energy_scales);
    }

private:
    AnalyzerSetup current;
    AnalyzerSetupCollection* descriptors;
};

} // namespace analysis
