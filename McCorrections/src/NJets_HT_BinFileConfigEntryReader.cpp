/*! Definition of the file configuration entry reader for DY sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/McCorrections/include/NJets_HT_BinFileConfigEntryReader.h"

namespace analysis {

namespace sample_merging{

NJets_HT_BinFileConfigEntryReader::NJets_HT_BinFileConfigEntryReader(NJets_HT_BinDescriptorCollection& _descriptors)
    : descriptors(&_descriptors)
{
}

void NJets_HT_BinFileConfigEntryReader::StartEntry(const std::string& name, const std::string& reference_name)
{
    ConfigEntryReader::StartEntry(name, reference_name);
    current = reference_name.size() ? descriptors->at(reference_name) : NJets_HT_BinFileDescriptor();
    current.name = name;
}

void NJets_HT_BinFileConfigEntryReader::EndEntry()
{
    CheckReadParamCounts("file_path", 0, Condition::greater_equal);
    CheckReadParamCounts("n_jet", 1, Condition::equal_to);
    CheckReadParamCounts("n_bjet", 1, Condition::equal_to);
    CheckReadParamCounts("n_ht", 1, Condition::equal_to);
    CheckReadParamCounts("fileType", 1, Condition::equal_to);

    (*descriptors)[current.name] = current;
}

void NJets_HT_BinFileConfigEntryReader::ReadParameter(const std::string& /*param_name*/,
                                                      const std::string& /*param_value*/,
                                                      std::istringstream& /*ss*/)
{
    ParseEntry("file_path", current.file_paths);
    ParseEntry("n_jet", current.n_jet);
    ParseEntry("n_bjet", current.n_bjet);
    ParseEntry("n_ht", current.n_ht);
    ParseEntry("fileType", current.fileType);
}

} //namespace sample_merging
} // namespace analysis
