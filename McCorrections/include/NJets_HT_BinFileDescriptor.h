/*! Definition of the file descriptor for DY and Wjets sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/PhysicalValue.h"
#include "SampleDescriptor.h"

namespace analysis {
namespace sample_merging{

struct NJets_HT_BinFileDescriptor {
    std::string name;
    std::vector<std::string> file_paths;
    FileType fileType;
    Range<size_t> n_jet;
    Range<size_t> n_bjet;
    Range<size_t> n_ht;

    PhysicalValue nu;
    PhysicalValue weight;
    std::string ref_sample;
    FileType ref_fileType;
    size_t inclusive_integral;

    NJets_HT_BinFileDescriptor();
    bool Contains(const ntuple::GenId& genId) const;

    static std::vector<NJets_HT_BinFileDescriptor> LoadConfig(const std::string& config_file);
    static std::ofstream SaveCfg(const std::string& output_file,
                                 const std::vector<NJets_HT_BinFileDescriptor>& output_bins);
};

using NJets_HT_BinDescriptorCollection = std::unordered_map<std::string, NJets_HT_BinFileDescriptor>;

} //namespace sample_merging

} // namespace analysis
