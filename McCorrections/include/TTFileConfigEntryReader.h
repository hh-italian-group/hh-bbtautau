/*! Definition of the file configuration entry reader for TTbar sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "hh-bbtautau/McCorrections/include/TTFileDescriptor.h"

namespace analysis {
namespace sample_merging{

class TTFileConfigEntryReader : public analysis::ConfigEntryReader {
public:
    TTFileConfigEntryReader(TTBinDescriptorCollection& _descriptors);

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override;
    virtual void EndEntry() override;
    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override;

private:
    TTBinDescriptor current;
    TTBinDescriptorCollection* descriptors;
};

} //namespace sample_merging
} // namespace analysis
