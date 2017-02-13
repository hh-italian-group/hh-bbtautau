/*! Definition of the file descriptor.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <list>
#include <istream>
#include <ostream>
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/PhysicalValue.h"

namespace analysis {

enum class FileType { inclusive, exclusive };
ENUM_NAMES(FileType) = {
    { FileType::inclusive, "inclusive" },
    { FileType::exclusive, "exclusive" }
};

struct DYBinDescriptor {
    std::string name;
    std::string file_path;
    FileType fileType;
    analysis::Range<int> n_jet;
    analysis::Range<int> n_bjet;
    analysis::Range<int> n_ht;

    PhysicalValue nu;
    double weight;
    std::string ref_sample;

    DYBinDescriptor()
        : n_jet(0,0), n_bjet(0,0),n_ht(0,0),nu(0.0, std::numeric_limits<double>::infinity()),
          weight(std::numeric_limits<double>::quiet_NaN()) {}

    static std::vector<DYBinDescriptor> LoadConfig(const std::string& config_file)
    {
        static std::mutex m;
        std::lock_guard<std::mutex> lock(m);
        std::vector<DYBinDescriptor> dyBinDescriptors;
        dyBinDescriptors.clear();
        std::ifstream cfg(config_file);
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            if (!cfgLine.size() || cfgLine.at(0) == '#') continue;
            std::istringstream ss(cfgLine);
            DYBinDescriptor descriptor;
            ss >> descriptor.n_jet >> descriptor.n_bjet >> descriptor.n_ht;
            dyBinDescriptors.push_back(descriptor);
        }
        return dyBinDescriptors;
    }
};

typedef std::unordered_map<std::string, DYBinDescriptor> DYBinDescriptorCollection;

} // namespace analysis
