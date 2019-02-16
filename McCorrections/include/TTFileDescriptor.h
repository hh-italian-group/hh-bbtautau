/*! Definition of the file descriptor for TTbar sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "SampleDescriptor.h"

namespace analysis {
namespace sample_merging{

struct TTBinDescriptor {
    std::string name;
    std::vector<std::string> file_paths;
    FileType fileType;
    Range<int> genType;

    PhysicalValue nu;
    PhysicalValue weight;
    size_t inclusive_integral;

    TTBinDescriptor();

    bool Contains(const analysis::GenEventType& genEventType) const;

    static std::vector<TTBinDescriptor> LoadConfig(const std::string& config_file);
    static std::ofstream SaveCfg(const std::string& output_file, const std::vector<TTBinDescriptor>& output_bins);
};

using TTBinDescriptorCollection = std::unordered_map<std::string, TTBinDescriptor>;

} //namespace sample_merging
} // namespace analysis
