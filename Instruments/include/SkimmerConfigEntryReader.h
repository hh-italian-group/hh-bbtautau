/*! Definition of classes to read TupleSkimmer configuration.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "SkimmerConfig.h"

namespace analysis {
namespace tuple_skimmer {

class SetupEntryReader : public ConfigEntryReaderT<Setup> {
public:
    using ConfigEntryReaderT<Setup>::ConfigEntryReaderT;

    virtual void EndEntry() override;
    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override;
};

class SkimJobEntryReader : public ConfigEntryReaderT<SkimJob> {
public:
    using ConfigEntryReaderT<SkimJob>::ConfigEntryReaderT;

    virtual void EndEntry() override;
    virtual void ReadParameter(const std::string& param_name, const std::string& param_value,
                               std::istringstream& /*ss*/) override;

private:
    void ParseFileDescriptor(const std::string& param_name, const std::string& param_value);
};

} // namespace tuple_skimmer
} // namespace analysis
