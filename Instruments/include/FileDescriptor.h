/*! Definition of the file descriptor.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <list>
#include <istream>
#include <ostream>
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/PhysicalValue.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"

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
    PhysicalValue weight;
    std::string ref_sample;
    FileType ref_fileType;

    DYBinDescriptor()
        : n_jet(0,0), n_bjet(0,0),n_ht(0,0),nu(0.0, std::numeric_limits<double>::infinity()),
          weight(std::numeric_limits<double>::quiet_NaN()) {}

    bool Contains(const ntuple::GenId& genId) const
    {
        return n_jet.Contains(genId.n_partons) &&
              n_bjet.Contains(genId.n_b_partons) &&
              n_ht.Contains(genId.ht10_bin);
    }

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

    static std::ofstream SaveCfg(const std::string& output_file, const std::vector<DYBinDescriptor>& output_bins)
    {
        std::ofstream cfg(output_file);
        if(cfg.fail())
            throw analysis::exception("Unable to create outputfile'");
        cfg.exceptions(std::ofstream::failbit | std::ofstream::badbit);

        cfg << "#n_jet_min n_jet_max n_bjet_min n_bjet_max ht_bin_min ht_bin_max weight rel_err_w "
               "nu rel_err_nu ref_sample\n";

        for(unsigned i = 0; i < output_bins.size(); ++i)
        {
            DYBinDescriptor output_bin = output_bins.at(i);
            cfg << output_bin.n_jet.min() << " " <<output_bin.n_jet.max()  << " " <<
                   output_bin.n_bjet.min() << " " << output_bin.n_bjet.max() << " " <<
                   output_bin.n_ht.min() << " " << output_bin.n_ht.max() << " " <<
                   output_bin.weight.GetValue() << " " << output_bin.weight.GetRelativeStatisticalError() << " " <<
                   output_bin.nu.GetValue() << " " << output_bin.nu.GetRelativeStatisticalError() <<
                " " << output_bin.ref_sample <<  "\n";
        }
        return cfg;
    }
};

typedef std::unordered_map<std::string, DYBinDescriptor> DYBinDescriptorCollection;

} // namespace analysis
