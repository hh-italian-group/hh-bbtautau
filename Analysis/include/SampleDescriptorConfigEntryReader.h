/*! Definition of the file configuration entry reader for analyzer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "SampleDescriptor.h"
#include "SampleDescriptorBase.h"
#include "SampleDescriptorBaseConfigEntryReader.h"

namespace analysis {


class SampleDescriptorConfigEntryReader : public SampleDescriptorBaseConfigEntryReader<SampleDescriptor> {
public:
    using Base = SampleDescriptorBaseConfigEntryReader<SampleDescriptor>;
    using Base::SampleDescriptorBaseConfigEntryReader;


    virtual void EndEntry() override
    {
        CheckReadParamCounts("file_path", 0, Condition::greater_equal);
        CheckReadParamCounts("file_path_pattern", 1, Condition::less_equal);
        CheckReadParamCounts("cross_section", 1, Condition::equal_to);
        CheckReadParamCounts("weight_file", 1, Condition::equal_to);
        CheckReadParamCounts("listSignalPoint", 0, Condition::greater_equal);
        CheckReadParamCounts("draw_ex", 0, Condition::greater_equal);
        CheckReadParamCounts("norm_sf", 0, Condition::greater_equal);
        CheckReadParamCounts("datacard_name_ex", 0, Condition::greater_equal);

        Base::EndEntry();
    }

    virtual void ReadParameter(const std::string& param_name, const std::string& param_value,
                               std::istringstream& ss) override
    {
        ParseEntry("file_path", current.file_paths);
        ParseEntry("file_path_pattern", current.file_path_pattern);
        ParseEntry<double,NumericalExpression>("cross_section", current.cross_section, [](double xs){return xs > 0;});
        ParseEntry("weight_file", current.weight_file);
        ParseEntry("listSignalPoint", current.listSignalPoints);
        ParseEntry("draw_ex", current.draw_ex);
        ParseEntry("norm_sf", current.norm_sf);
        ParseEntry("datacard_name_ex", current.datacard_name_ex);

        Base::ReadParameter(param_name,param_value,ss);
    }

};

} // namespace analysis
