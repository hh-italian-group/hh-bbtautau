/*! Definition of the file descriptor.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <list>
#include <istream>
#include <ostream>
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"

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

    PhysicalValue nu{0, inf};
    double weight;
    std::string ref_sample;

    FileDescriptor()
        : n_jet(0,0), n_bjet(0,0),n_ht(0,0) {}

    static std::vector<DYBinDescriptor> LoadConfig(const std::string& config_file)
    {
        static std::mutex m;
        std::lock_guard<std::mutex> lock(m);
        JetParametersMap jetParameters;
        jetParameters.clear();
        size_t line_number = 0;
        std::ifstream cfg(config_file);
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            if (!cfgLine.size() || cfgLine.at(0) == '#') continue;
            ++line_number;
            std::istringstream ss(cfgLine);
            JetParameters jetParam;
            ss >> jetParam.n_jet >> jetParam.n_bjet >> jetParam.n_ht;
            jetParameters[line_number] = jetParam;
        }
        return jetParameters;
    }
};

typedef std::unordered_map<std::string, FileDescriptor> FileDescriptorCollection;

} // namespace analysis
