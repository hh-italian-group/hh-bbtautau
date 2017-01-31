/*! Definition of the file descriptor.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <list>
#include <istream>
#include <ostream>
#include "AnalysisTools/Core/include/EnumNameMap.h"

namespace analysis {

enum class FileType { inclusive, exclusive };
ENUM_NAMES(FileType) = {
    { FileType::inclusive, "inclusive" },
    { FileType::exclusive, "exclusive" }
};

struct FileDescriptor {
    std::string name;
    std::string file_path;
    FileType fileType;
    double n_jet_min;
    double n_jet_max;
    double n_bjet_min;
    double n_bjet_max;
    double n_ht_min;
    double n_ht_max;

    FileDescriptor()
        : n_jet_min(0.0), n_jet_max(0.0), n_bjet_min(0.0), n_bjet_max(0.0), n_ht_min(0.0),
            n_ht_max(0.0){}
};

typedef std::unordered_map<std::string, FileDescriptor> FileDescriptorCollection;

} // namespace analysis
