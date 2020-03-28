/*! The ttbar weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/McCorrections/include/TTbar_weight.h"

namespace analysis {
namespace mc_corrections {

TTbar_weight::TTbar_weight(const std::string& ttbar_weight_file_name)
{
    std::vector<analysis::sample_merging::TTBinDescriptor> ttbar_descriptors =
            analysis::sample_merging::TTBinDescriptor::LoadConfig(ttbar_weight_file_name);
    for (unsigned n = 0; n < ttbar_descriptors.size(); ++n){
        const analysis::sample_merging::TTBinDescriptor ttbin_descriptor = ttbar_descriptors.at(n);
        genEventType_weight_map[ttbin_descriptor.genType.min()]
                = ttbin_descriptor.weight.GetValue()/ttbin_descriptor.inclusive_integral;
    }
}

double TTbar_weight::Get(EventInfo& eventInfo) const { return Get(eventInfo->genEventType); }
double TTbar_weight::Get(const ntuple::ExpressEvent& event) const { return Get(event.genEventType); }

double TTbar_weight::Get(int genEventType) const
{
    auto iter = genEventType_weight_map.find(genEventType);
    if(iter != genEventType_weight_map.end())
        return iter->second;
    throw exception("ttbar merge weight not found for genEventType = %1%.") % genEventType;
}

} // namespace mc_corrections
} // namespace analysis
