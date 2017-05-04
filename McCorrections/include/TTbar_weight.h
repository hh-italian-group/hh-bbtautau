/*! The ttbar weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "../../Instruments/include/TTFileDescriptor.h"


namespace analysis {
namespace mc_corrections {

class TTbar_weight {
public:
    using Event = ntuple::Event;

    TTbar_weight(const std::string& ttbar_weight_file_name) :
    {
        ttbar_descriptors = analysis::sample_merging::TTBinDescriptor::LoadConfig(ttbar_weight_file_name);
        for (unsigned n = 0; n < ttbar_descriptors.size(); ++n){
            const analysis::sample_merging::TTBinDescriptor ttbin_descriptor = ttbar_descriptors.at(n);
            genEventType_weight_map[ttbin_descriptor.genType] = ttbin_descriptor.weight.GetValue();
        }
    }

    double Get(const Event& event)
    {
        double weight = 1;
        for (auto iter : genEventType_weight_map){
            Range<int> genEventType = iter.first;
            double tt_weight = iter.second;
            if (genEventType.Contains(event.genEventType))
                weight = tt_weight;
        }

        return weight;
    }

private:
    std::vector<analysis::sample_merging::TTBinDescriptor> ttbar_descriptors;
    std::map<Range<int>, double> genEventType_weight_map;


};

} // namespace mc_corrections
} // namespace analysis
