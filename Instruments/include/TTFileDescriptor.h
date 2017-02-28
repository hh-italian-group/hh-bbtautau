/*! Definition of the file descriptor for TTbar sample merging.
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

struct TTBinDescriptor {
    std::string name;
    std::vector<std::string> file_paths;
    FileType fileType;
    Range<int> genType;

    PhysicalValue nu;
    PhysicalValue weight;
    std::string ref_sample;
    FileType ref_fileType;

    TTBinDescriptor()
        : genType(0,0),nu(0.0, std::numeric_limits<double>::infinity()),
          weight(std::numeric_limits<double>::quiet_NaN()) {}

    bool Contains(const ntuple::GenId& genId) const  //?????
    {
        return n_jet.Contains(genId.n_partons) &&
              n_bjet.Contains(genId.n_b_partons) &&
              n_ht.Contains(genId.ht10_bin);
    }

    static std::vector<TTBinDescriptor> LoadConfig(const std::string& config_file)
    {
        std::vector<TTBinDescriptor> ttBinDescriptors;
        ttBinDescriptors.clear();
        std::ifstream cfg(config_file);
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            if (!cfgLine.size() || cfgLine.at(0) == '#') continue;
            auto columns = ConfigEntryReader::ParseOrderedParameterList(cfgLine, true);
            std::istringstream ss(cfgLine);
            TTBinDescriptor descriptor;
            if(columns.size() >= 2)
                ss >> descriptor.genType;
            if (columns.size() >= 7){
                double col_weight = analysis::Parse<double>(columns.at(2));
                double col_weight_err = analysis::Parse<double>(columns.at(3))*col_weight;
                descriptor.weight = PhysicalValue(col_weight,col_weight_err);
                double col_nu = analysis::Parse<double>(columns.at(4));
                double col_nu_err = analysis::Parse<double>(columns.at(5))*col_nu;
                descriptor.nu = PhysicalValue(col_nu,col_nu_err);
                descriptor.ref_sample = columns.at(6);
            }
            if(columns.size() != 2 && columns.size() != 7)
                throw exception("Bad configuration file.");
            ttBinDescriptors.push_back(descriptor);
        }
        return ttBinDescriptors;
    }


    static std::ofstream SaveCfg(const std::string& output_file, const std::vector<TTBinDescriptor>& output_bins)
    {
        std::ofstream cfg(output_file);
        if(cfg.fail())
            throw analysis::exception("Unable to create outputfile'");
        cfg.exceptions(std::ofstream::failbit | std::ofstream::badbit);

        static const std::vector<int> column_widths = { 9,9,15,15,15,15,20 };

        cfg << std::setw(column_widths.at(0)) << "#genType_min " <<
               std::setw(column_widths.at(1)) << "genType_max " <<
               std::setw(column_widths.at(2)) << "weight " <<
               std::setw(column_widths.at(3)) << " rel_err_w " <<
               std::setw(column_widths.at(4)) << " nu " <<
               std::setw(column_widths.at(5)) << " rel_err_nu  " <<
               std::setw(column_widths.at(6)) << " ref_sample \n";

        for(auto& output_bin : output_bins)
        {
            cfg << std::setw(column_widths.at(0)) << output_bin.genType.min() << " " <<
                   std::setw(column_widths.at(1)) << output_bin.genType.max()  << " " <<
                   std::setw(column_widths.at(6)) << output_bin.weight.GetValue() << " " <<
                   std::setw(column_widths.at(7)) << output_bin.weight.GetRelativeStatisticalError() << " " <<
                   std::setw(column_widths.at(8)) << output_bin.nu.GetValue() << " " <<
                   std::setw(column_widths.at(9)) << output_bin.nu.GetRelativeStatisticalError() << " " <<
                   std::setw(column_widths.at(10)) << output_bin.ref_sample <<  "\n";
        }
        return cfg;
    }
};

using TTBinDescriptorCollection = std::unordered_map<std::string, TTBinDescriptor>;

} //namespace sample_merging

} // namespace analysis
