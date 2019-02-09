/*! Definition of the file configuration entry reader for DY sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "hh-bbtautau/McCorrections/include/NJets_HT_BinFileDescriptor.h"

namespace analysis {

namespace sample_merging{

class NJets_HT_BinFileConfigEntryReader : public analysis::ConfigEntryReader {
public:
    NJets_HT_BinFileConfigEntryReader(NJets_HT_BinDescriptorCollection& _descriptors);

    virtual void StartEntry(const std::string& name, const std::string& reference_name) override;
    virtual void EndEntry() override;
    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override;

private:
    NJets_HT_BinFileDescriptor current;
    NJets_HT_BinDescriptorCollection* descriptors;
};

} //namespace sample_merging

} // namespace analysis
