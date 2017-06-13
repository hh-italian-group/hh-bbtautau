/*! Definition of the file descriptor for DY and Wjets sample merging.
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
#include "Instruments/include/SampleDescriptor.h"

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

    NJets_HT_BinFileDescriptor()
        : n_jet(0,0), n_bjet(0,0),n_ht(0,0),nu(0.0, std::numeric_limits<double>::infinity()),
          weight(std::numeric_limits<double>::quiet_NaN()) {}

    bool Contains(const ntuple::GenId& genId) const
    {
        return n_jet.Contains(genId.n_partons) &&
              n_bjet.Contains(genId.n_b_partons) &&
              n_ht.Contains(genId.ht10_bin);
    }

    static std::vector<NJets_HT_BinFileDescriptor> LoadConfig(const std::string& config_file)
    {
        std::vector<NJets_HT_BinFileDescriptor> jetBinDescriptors;
        std::ifstream cfg(config_file);
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            if (!cfgLine.size() || cfgLine.at(0) == '#') continue;
            auto columns = ConfigEntryReader::ParseOrderedParameterList(cfgLine, true);
            std::istringstream ss(cfgLine);
            NJets_HT_BinFileDescriptor descriptor;
            if(columns.size() >= 6)
                ss >> descriptor.n_jet >> descriptor.n_bjet >> descriptor.n_ht;
            if (columns.size() >= 12){
                double col_weight = analysis::Parse<double>(columns.at(6));
                double col_weight_err = analysis::Parse<double>(columns.at(7))*col_weight;
                descriptor.weight = PhysicalValue(col_weight,col_weight_err);
                double col_nu = analysis::Parse<double>(columns.at(8));
                double col_nu_err = analysis::Parse<double>(columns.at(9))*col_nu;
                descriptor.nu = PhysicalValue(col_nu,col_nu_err);
                descriptor.inclusive_integral = static_cast<size_t>(analysis::Parse<double>(columns.at(10)));
                descriptor.ref_sample = columns.at(11);

            }
            if(columns.size() != 6 && columns.size() != 12)
                throw exception("Bad configuration file.");
            jetBinDescriptors.push_back(descriptor);
        }
        return jetBinDescriptors;
    }


    static std::ofstream SaveCfg(const std::string& output_file, const std::vector<NJets_HT_BinFileDescriptor>& output_bins)
    {
        std::ofstream cfg(output_file);
        if(cfg.fail())
            throw analysis::exception("Unable to create outputfile'");
        cfg.exceptions(std::ofstream::failbit | std::ofstream::badbit);

        std::vector<std::string> headers = { "#n_jet_min","n_jet_max","n_bjet_min","n_bjet_max",
                                           "ht_bin_min","ht_bin_max","weight","rel_err_w","nu",
                                           "rel_err_nu","incl_integral","ref_sample"};
        std::cout << "headers size:" << headers.size() << std::endl;

        std::vector<int> column_widths_2(headers.size());

        std::cout << "column_widths size:" << column_widths_2.size() << std::endl;

        column_widths_2[0] = static_cast<int>(headers[0].size());
        for(unsigned n = 1; n <= 5; ++n) {
            column_widths_2[n] = static_cast<int>(headers[n].size()+1);
            std::cout << "column_widths n:" << column_widths_2[n] << std::endl;
        }

//        for(unsigned h = 6; h <= 9; ++h){
        for(unsigned h = 6; h < headers.size(); ++h){
          column_widths_2[h] = std::numeric_limits<double>::digits10 + 6;
          std::cout << "column_widths n:" << column_widths_2[h] << std::endl;
        }

//        for(unsigned h = 10; h < headers.size(); ++h){
//          column_widths_2[h] = headers[h].size()+1;
//          std::cout << "column_widths n:" << column_widths_2[h] << std::endl;
//        }

//        static const std::vector<int> column_widths = { 10,9,10,10,10,10,21,21,21,21,21,26 };

        cfg << std::setw(column_widths_2.at(0)) << "#n_jet_min" <<
               std::setw(column_widths_2.at(1)) << "n_jet_max" <<
               std::setw(column_widths_2.at(2)) << "n_bjet_min" <<
               std::setw(column_widths_2.at(3)) << "n_bjet_max" <<
               std::setw(column_widths_2.at(4)) << "ht_bin_min" <<
               std::setw(column_widths_2.at(5)) << "ht_bin_max" <<
               std::setw(column_widths_2.at(6)) << "weight" <<
               std::setw(column_widths_2.at(7)) << "rel_err_w" <<
               std::setw(column_widths_2.at(8)) << "nu" <<
               std::setw(column_widths_2.at(9)) << "rel_err_nu" <<
               std::setw(column_widths_2.at(10)) << "incl_integral" <<
               std::setw(column_widths_2.at(11)) << "ref_sample\n";

        cfg << std::setprecision(std::numeric_limits<double>::digits10);
        for(auto& output_bin : output_bins)
        {
            cfg << std::setw(column_widths_2.at(0)) << output_bin.n_jet.min() <<
                   std::setw(column_widths_2.at(1)) << output_bin.n_jet.max()  <<
                   std::setw(column_widths_2.at(2)) << output_bin.n_bjet.min() <<
                   std::setw(column_widths_2.at(3)) << output_bin.n_bjet.max() <<
                   std::setw(column_widths_2.at(4)) << output_bin.n_ht.min() <<
                   std::setw(column_widths_2.at(5)) << output_bin.n_ht.max() <<
                   std::setw(column_widths_2.at(6)) << output_bin.weight.GetValue() <<
                   std::setw(column_widths_2.at(7)) << output_bin.weight.GetRelativeStatisticalError() <<
                   std::setw(column_widths_2.at(8)) << output_bin.nu.GetValue() <<
                   std::setw(column_widths_2.at(9)) << output_bin.nu.GetRelativeStatisticalError() <<
                   std::setw(column_widths_2.at(10)) << output_bin.inclusive_integral <<
                   std::setw(column_widths_2.at(11)) << output_bin.ref_sample <<  "\n";
        }
        return cfg;
    }
};

using NJets_HT_BinDescriptorCollection = std::unordered_map<std::string, NJets_HT_BinFileDescriptor>;

} //namespace sample_merging

} // namespace analysis
