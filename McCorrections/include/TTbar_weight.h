/*! The ttbar weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "hh-bbtautau/Instruments/include/TTFileDescriptor.h"
#include "h-tautau/McCorrections/include/WeightProvider.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"

namespace analysis {
namespace mc_corrections {

class TTbar_weight : public IWeightProvider {
public:
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

    virtual double Get(const ntuple::Event& event) const override
    {
        auto iter = genEventType_weight_map.find(event.genEventType);
        if(iter != genEventType_weight_map.end())
            return iter->second;
        throw exception("ttbar merge weight not found for genEventType = %1%.") % event.genEventType;
    }

    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override
    {
        throw exception("ExpressEvent is not supported in TTbar_weight::Get.");
    }

private:
    std::map<int, double> genEventType_weight_map;
};

class TTbar_weight_multiChannel : public IWeightProvider {
public:
    TTbar_weight_multiChannel(std::initializer_list<std::pair<Channel, std::string>>&& channel_weight_files)
    {
        for(const auto& pair : channel_weight_files)
            weights[pair.first] = std::make_shared<TTbar_weight>(pair.second);
    }

    virtual double Get(const ntuple::Event& event) const override
    {
        const Channel ch = static_cast<Channel>(event.channelId);
        auto iter = weights.find(ch);
        if(iter == weights.end())
            throw exception("Can't find ttbar weight for %1% channel.") % ch;
        return iter->second->Get(event);
    }

    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override
    {
        throw exception("ExpressEvent is not supported in TTbar_weight_multiChannel::Get.");
    }

private:
    std::map<Channel, std::shared_ptr<TTbar_weight>> weights;
};

} // namespace mc_corrections
} // namespace analysis
