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

struct FileDescriptor {
    std::string name;
    std::string file_path;
    FileType fileType;
    analysis::Range<int> n_jet;
    analysis::Range<int> n_bjet;
    analysis::Range<int> n_ht;

    FileDescriptor()
        : n_jet(0,0), n_bjet(0,0),n_ht(0,0) {}
};

typedef std::unordered_map<std::string, FileDescriptor> FileDescriptorCollection;

} // namespace analysis
