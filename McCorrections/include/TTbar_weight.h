/*! The ttbar weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "hh-bbtautau/Instruments/include/TTFileDescriptor.h"


namespace analysis {
namespace mc_corrections {

class TTbar_weight {
public:
    using Event = ntuple::Event;

    TTbar_weight(const std::string& ttbar_weight_file_name)
    {
        std::vector<analysis::sample_merging::TTBinDescriptor> ttbar_descriptors =
                analysis::sample_merging::TTBinDescriptor::LoadConfig(ttbar_weight_file_name);
        for (unsigned n = 0; n < ttbar_descriptors.size(); ++n){
            const analysis::sample_merging::TTBinDescriptor ttbin_descriptor = ttbar_descriptors.at(n);
            genEventType_weight_map[ttbin_descriptor.genType.min()]
                    = ttbin_descriptor.weight.GetValue()/ttbin_descriptor.inclusive_integral;
        }
    }

    template<typename Event>
    double Get(const Event& event) const
    {
        auto iter = genEventType_weight_map.find(event.genEventType);
        if(iter != genEventType_weight_map.end())
            return iter->second;
        throw analysis::exception("ttbar merge weight not found.");
    }

private:
    std::map<int, double> genEventType_weight_map;
};

} // namespace mc_corrections
} // namespace analysis
