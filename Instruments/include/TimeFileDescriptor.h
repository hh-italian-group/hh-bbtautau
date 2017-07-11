/*! Definition of the file descriptor for TTbar sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <list>
#include <istream>
#include <ostream>
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/PhysicalValue.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/ConfigReader.h"

namespace analysis {

struct TimeFileDescriptor {

    std::string file;
    double scale_factor;
    size_t n_evt_per_job_prod_v2;
    size_t n_evt_per_job_prod_v3;

    TimeFileDescriptor() {}

    static std::vector<TimeFileDescriptor> LoadConfig(const std::string& config_file)
    {
        std::vector<TimeFileDescriptor> descriptors;
        std::ifstream cfg(config_file);
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            if (!cfgLine.size() || cfgLine.at(0) == '#') continue;
            auto columns = ConfigEntryReader::ParseOrderedParameterList(cfgLine, true);
            std::istringstream ss(cfgLine);
            TimeFileDescriptor descriptor;
            if(columns.size() >= 1)
                ss >> descriptor.file;
            if (columns.size() >= 4){
                descriptor.scale_factor = analysis::Parse<double>(columns.at(1));
                descriptor.n_evt_per_job_prod_v2 = analysis::Parse<size_t>(columns.at(2));
                descriptor.n_evt_per_job_prod_v3 = analysis::Parse<size_t>(columns.at(3));
            }
            if(columns.size() != 1 && columns.size() != 4)
                throw exception("Bad configuration file.");
            descriptors.push_back(descriptor);
        }
        return descriptors;
    }


    static std::ofstream SaveCfg(const std::string& output_file, const std::vector<TimeFileDescriptor>& outputs)
    {
        std::ofstream cfg(output_file);
        if(cfg.fail())
            throw analysis::exception("Unable to create outputfile'");
        cfg.exceptions(std::ofstream::failbit | std::ofstream::badbit);

        std::vector<std::string> headers = { "#file","scale_factor","n_evt_per_job_prod_v2","n_evt_per_job_prod_v3"};
        std::cout << "headers size:" << headers.size() << std::endl;

        std::vector<int> column_widths_2(headers.size());

        std::cout << "column_widths size:" << column_widths_2.size() << std::endl;

        column_widths_2[0] = static_cast<int>(headers[0].size()+30);
        column_widths_2[1] = static_cast<int>(headers[1].size()+10);
        column_widths_2[2] = static_cast<int>(headers[2].size()+2);
        column_widths_2[3] = static_cast<int>(headers[3].size()+2);

        cfg << std::setw(column_widths_2.at(0)) << "#file" <<
               std::setw(column_widths_2.at(1)) << "scale_factor" <<
               std::setw(column_widths_2.at(2)) << "n_evt_per_job_prod_v2" <<
               std::setw(column_widths_2.at(3)) << "n_evt_per_job_prod_v3\n";

        cfg << std::setprecision(std::numeric_limits<double>::digits10);
        for(auto& output : outputs)
        {
            cfg << std::setw(column_widths_2.at(0)) << output.file <<
                   std::setw(column_widths_2.at(1)) << output.scale_factor  <<
                   std::setw(column_widths_2.at(2)) << output.n_evt_per_job_prod_v2 <<
                   std::setw(column_widths_2.at(3)) << output.n_evt_per_job_prod_v3 << "\n";
        }
        return cfg;
    }

};

using TimeFileDescriptorCollection = std::unordered_map<std::string, TimeFileDescriptor>;

} // namespace analysis
