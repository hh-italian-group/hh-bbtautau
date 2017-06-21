/*! Definition of the file configuration entry reader for analyzer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "SampleDescriptor.h"
#include "SampleDescriptorBase.h"

namespace analysis {


class SampleDescriptorConfigEntryReader {
public:
    SampleDescriptorConfigEntryReader(SampleDescriptorCollection& _descriptors,
                                      SampleDescriptorBaseCollection& _baseDescriptors) :
        descriptorsBase(&_baseDescriptors) ,
        descriptors(&_descriptors) {}

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override
    {
        ConfigEntryReader::StartEntry(name, reference_name);
        current = reference_name.size() ? descriptors->at(reference_name) : SampleDescriptor();
        current.name = name;

        currentBase = reference_name.size() ? descriptorsBase->at(reference_name) : SampleDescriptorBase();
        currentBase.name = name;
    }

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

        (*descriptors)[current.name] = current;

        CheckReadParamCounts("title", 1, Condition::equal_to);
        CheckReadParamCounts("color", 1, Condition::equal_to);
        CheckReadParamCounts("draw", 1, Condition::equal_to);
        CheckReadParamCounts("channel", 1, Condition::equal_to);
        CheckReadParamCounts("categoryType", 1, Condition::equal_to);
        CheckReadParamCounts("create_hist", 1, Condition::equal_to);
        CheckReadParamCounts("datacard_name", 1, Condition::equal_to);

        (*descriptorsBase)[currentBase.name] = currentBase;
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("file_path", current.file_paths);
        ParseEntry("file_path_pattern", current.file_path_pattern);
        ParseEntry<double,NumericalExpression>("cross_section", current.cross_section, [](double xs){return xs > 0;});
        ParseEntry("weight_file", current.weight_file);
        ParseEntry("listSignalPoint", current.listSignalPoints);
        ParseEntry("draw_ex", current.draw_ex);
        ParseEntry("norm_sf", current.norm_sf);
        ParseEntry("datacard_name_ex", current.datacard_name_ex);

        ParseEntry("title", currentBase.title);
        ParseEntry("color", currentBase.color);
        ParseEntry("draw", currentBase.draw);
        ParseEntry("channel", currentBase.channel);
        ParseEntry("categoryType", currentBase.categoryType);
        ParseEntry("create_hist", currentBase.create_hist);
        ParseEntry("datacard_name", currentBase.datacard_name);
    }

private:
    SampleDescriptor current;
    SampleDescriptorCollection* descriptors;
    SampleDescriptorBase currentBase;
    SampleDescriptorBaseCollection* descriptorsBase;
};

} // namespace analysis
