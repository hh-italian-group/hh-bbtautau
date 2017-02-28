/*! Definition of the file descriptor for DY sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <list>
#include <istream>
#include <ostream>
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/PhysicalValue.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"
#include "AnalysisTools/Core/include/ConfigReader.h"

namespace analysis {

namespace sample_merging{

enum class FileType { inclusive, exclusive };
ENUM_NAMES(FileType) = {
    { FileType::inclusive, "inclusive" },
    { FileType::exclusive, "exclusive" }
};

struct DYBinDescriptor {
    std::string name;
    std::vector<std::string> file_paths;
    FileType fileType;
    Range<int> n_jet;
    Range<int> n_bjet;
    Range<int> n_ht;

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
        std::vector<DYBinDescriptor> dyBinDescriptors;
        dyBinDescriptors.clear();
        std::ifstream cfg(config_file);
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            if (!cfgLine.size() || cfgLine.at(0) == '#') continue;
            auto columns = ConfigEntryReader::ParseOrderedParameterList(cfgLine, true);
            std::istringstream ss(cfgLine);
            DYBinDescriptor descriptor;
            if(columns.size() >= 6)
                ss >> descriptor.n_jet >> descriptor.n_bjet >> descriptor.n_ht;
            if (columns.size() >= 11){
                double col_weight = analysis::Parse<double>(columns.at(6));
                double col_weight_err = analysis::Parse<double>(columns.at(7))*col_weight;
                descriptor.weight = PhysicalValue(col_weight,col_weight_err);
                double col_nu = analysis::Parse<double>(columns.at(8));
                double col_nu_err = analysis::Parse<double>(columns.at(9))*col_nu;
                descriptor.nu = PhysicalValue(col_nu,col_nu_err);
                descriptor.ref_sample = columns.at(10);
            }
            if(columns.size() != 6 && columns.size() != 11)
                throw exception("Bad configuration file.");
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

        static const std::vector<int> column_widths = { 9,9,9,9,9,11,15,15,15,15,20 };

        cfg << std::setw(column_widths.at(0)) << "#n_jet_min " <<
               std::setw(column_widths.at(1)) << "n_jet_max " <<
               std::setw(column_widths.at(2)) << "n_bjet_min " <<
               std::setw(column_widths.at(3)) << "n_bjet_max " <<
               std::setw(column_widths.at(4)) << "ht_bin_min " <<
               std::setw(column_widths.at(5)) << "ht_bin_max " <<
               std::setw(column_widths.at(6)) << "weight " <<
               std::setw(column_widths.at(7)) << " rel_err_w " <<
               std::setw(column_widths.at(8)) << " nu " <<
               std::setw(column_widths.at(9)) << " rel_err_nu  " <<
               std::setw(column_widths.at(10)) << " ref_sample \n";

        for(auto& output_bin : output_bins)
        {
            cfg << std::setw(column_widths.at(0)) << output_bin.n_jet.min() << " " <<
                   std::setw(column_widths.at(1)) << output_bin.n_jet.max()  << " " <<
                   std::setw(column_widths.at(2)) << output_bin.n_bjet.min() << " " <<
                   std::setw(column_widths.at(3)) << output_bin.n_bjet.max() << " " <<
                   std::setw(column_widths.at(4)) << output_bin.n_ht.min() << " " <<
                   std::setw(column_widths.at(5)) << output_bin.n_ht.max() << " " <<
                   std::setw(column_widths.at(6)) << output_bin.weight.GetValue() << " " <<
                   std::setw(column_widths.at(7)) << output_bin.weight.GetRelativeStatisticalError() << " " <<
                   std::setw(column_widths.at(8)) << output_bin.nu.GetValue() << " " <<
                   std::setw(column_widths.at(9)) << output_bin.nu.GetRelativeStatisticalError() << " " <<
                   std::setw(column_widths.at(10)) << output_bin.ref_sample <<  "\n";
        }
        return cfg;
    }
};

using DYBinDescriptorCollection = std::unordered_map<std::string, DYBinDescriptor>;

} //namespace sample_merging

} // namespace analysis
