/*! Definition of the file configuration entry reader for TTbar sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/McCorrections/include/TTFileConfigEntryReader.h"

namespace analysis {
namespace sample_merging{

TTFileConfigEntryReader::TTFileConfigEntryReader(TTBinDescriptorCollection& _descriptors)
    : descriptors(&_descriptors) {}

void TTFileConfigEntryReader::StartEntry(const std::string& name, const std::string& reference_name)
{
    ConfigEntryReader::StartEntry(name, reference_name);
    current = reference_name.size() ? descriptors->at(reference_name) : TTBinDescriptor();
    current.name = name;
}

void TTFileConfigEntryReader::EndEntry()
{
    CheckReadParamCounts("file_path", 0, Condition::greater_equal);
    CheckReadParamCounts("genType", 1, Condition::equal_to);
    CheckReadParamCounts("fileType", 1, Condition::equal_to);

    (*descriptors)[current.name] = current;
}

void TTFileConfigEntryReader::ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                                            std::istringstream& /*ss*/)
{
    ParseEntry("file_path", current.file_paths);
    ParseEntry("genType", current.genType);
    ParseEntry("fileType", current.fileType);
}

} //namespace sample_merging
} // namespace analysis
