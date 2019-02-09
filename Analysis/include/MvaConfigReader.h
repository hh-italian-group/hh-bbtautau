/*! Definition of MvaVariablesStudy class, the main class for Mva studies.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "hh-bbtautau/Analysis/include/MvaConfiguration.h"

namespace  analysis {
namespace mva_study {
class MvaConfigReader : public ConfigEntryReaderT<MvaSetup> {
public:
    using ConfigEntryReaderT<MvaSetup>::ConfigEntryReaderT;

    virtual void EndEntry() override;
    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override;
};

class SampleConfigReader : public analysis::ConfigEntryReaderT<SampleEntryList> {
public:
    using ConfigEntryReaderT<SampleEntryList>::ConfigEntryReaderT;

    virtual void EndEntry() override;
    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override;
};

}
}
