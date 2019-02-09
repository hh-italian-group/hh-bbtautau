/*! Definition of the file descriptor for TTbar sample merging.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/McCorrections/include/TTFileDescriptor.h"

#include <fstream>

namespace analysis {
namespace sample_merging{

TTBinDescriptor::TTBinDescriptor()
    : genType(0,0),nu(0.0, std::numeric_limits<double>::infinity()),
      weight(std::numeric_limits<double>::quiet_NaN()) {}

bool TTBinDescriptor::Contains(const analysis::GenEventType& genEventType) const
{
    return genType.Contains(static_cast<int>(genEventType));
}

std::vector<TTBinDescriptor> TTBinDescriptor::LoadConfig(const std::string& config_file)
{
    std::vector<TTBinDescriptor> ttBinDescriptors;
    ttBinDescriptors.clear();
    std::ifstream cfg(config_file);
    while (cfg.good()) {
        std::string cfgLine;
        std::getline(cfg,cfgLine);
        if (!cfgLine.size() || cfgLine.at(0) == '#') continue;
        auto columns = SplitValueList(cfgLine, true);
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
            descriptor.inclusive_integral = static_cast<size_t>(analysis::Parse<double>(columns.at(6)));
        }
        if(columns.size() != 2 && columns.size() != 7)
            throw exception("Bad configuration file.");
        ttBinDescriptors.push_back(descriptor);
    }
    return ttBinDescriptors;
}

std::ofstream TTBinDescriptor::SaveCfg(const std::string& output_file, const std::vector<TTBinDescriptor>& output_bins)
{
    std::ofstream cfg(output_file);
    if(cfg.fail())
        throw analysis::exception("Unable to create outputfile'");
    cfg.exceptions(std::ofstream::failbit | std::ofstream::badbit);

    std::vector<std::string> headers = { "#genType_min","genType_max","weight","rel_err_w",
                                       "nu","rel_err_nu","inclusive_integral"};
    std::cout << "headers size:" << headers.size() << std::endl;

    std::vector<int> column_widths(headers.size());

    std::cout << "column_widths size:" << column_widths.size() << std::endl;

    column_widths[0] = static_cast<int>(headers[0].size());
    column_widths[1] = static_cast<int>(headers[1].size()+1);

    for(unsigned h = 2; h < headers.size(); ++h){
      column_widths[h] = std::numeric_limits<double>::digits10 + 6;
      std::cout << "column_widths n:" << column_widths[h] << std::endl;
    }

    //static const std::vector<int> column_widths = { 9,9,15,15,15,15,15 };

    cfg << std::setw(column_widths.at(0)) << "#genType_min" <<
           std::setw(column_widths.at(1)) << "genType_max" <<
           std::setw(column_widths.at(2)) << "weight" <<
           std::setw(column_widths.at(3)) << " rel_err_w" <<
           std::setw(column_widths.at(4)) << " nu" <<
           std::setw(column_widths.at(5)) << " rel_err_nu" <<
           std::setw(column_widths.at(6)) << " inclusive_integral\n";

    cfg << std::setprecision(std::numeric_limits<double>::digits10);

    for(auto& output_bin : output_bins)
    {
        cfg << std::setw(column_widths.at(0)) << output_bin.genType.min() <<
               std::setw(column_widths.at(1)) << output_bin.genType.max()  <<
               std::setw(column_widths.at(2)) << output_bin.weight.GetValue() <<
               std::setw(column_widths.at(3)) << output_bin.weight.GetRelativeStatisticalError() <<
               std::setw(column_widths.at(4)) << output_bin.nu.GetValue() <<
               std::setw(column_widths.at(5)) << output_bin.nu.GetRelativeStatisticalError() <<
               std::setw(column_widths.at(6)) << output_bin.inclusive_integral << "\n";
    }
    return cfg;
}

} //namespace sample_merging
} // namespace analysis
