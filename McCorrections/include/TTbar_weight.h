/*! The ttbar weight.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"


namespace analysis {
namespace mc_corrections {

class TTbar_weight {
public:
    using Event = ntuple::Event;

    TTbar_weight() : { }

    double GetWeight(const std::string& config_file, const Event& event)
    {
        Range<int> genType;
        double ttbar_weight = 1;
        std::ifstream cfg(config_file);
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            if (!cfgLine.size() || cfgLine.at(0) == '#') continue;
            auto columns = ConfigEntryReader::ParseOrderedParameterList(cfgLine, true);
            std::istringstream ss(cfgLine);
            if(columns.size() >= 2)
                ss >> genType;
            if (columns.size() >= 7){
                double col_weight = analysis::Parse<double>(columns.at(2));
                if (genType.Contains(static_cast<int>(event.genEventType)))
                    ttbar_weight = col_weight;
            }
            if(columns.size() != 2 && columns.size() != 7)
                throw exception("Bad configuration file.");
        }

        return ttbar_weight;
    }

};

} // namespace mc_corrections
} // namespace analysis
