/*! Collection of histogram containers for event analyzers.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <tuple>
#include "EventAnalyzerData.h"

namespace analysis {

class EventAnalyzerDataCollection {
public:
    using Data = EventAnalyzerData;
    using DataPtr = std::shared_ptr<Data>;
    using DataId = EventAnalyzerDataId;
    using DataMap = std::map<DataId, DataPtr>;
    using Tuple = bbtautau::AnaTuple;
    using NameSet = std::set<std::string>;
    using DescSet = PropertyConfigReader::ItemCollection;
    using SampleUnc = ModellingUncertainty::SampleUnc;
    using MucPtr = std::shared_ptr<ModellingUncertaintyCollection>;

    explicit EventAnalyzerDataCollection(std::shared_ptr<TFile> _file, Channel _channel, const Tuple* _tuple,
            const NameSet& _histNames, const DescSet& _histDescs, const NameSet& _backgrounds = {},
            MucPtr _unc_collection = MucPtr()) :
        file(_file), channel(_channel), tuple(_tuple), histNames(_histNames), histDescs(_histDescs),
        backgrounds(_backgrounds), unc_collection(_unc_collection)
    {}

    Data& Get(const DataId& id)
    {
        auto& anaData = anaDataMap[id];
        if(!anaData)
            anaData = Make(id);
        return *anaData;
    }

    const DataMap& GetAll() const { return anaDataMap; }
    Channel ChannelId() const { return channel; }
    bool ReadMode() const { return tuple != nullptr; }

    void Fill(const DataId& id, double weight)
    {
        Get(id).Fill(weight);
    }

private:
    DataPtr Make(const DataId& id) const
    {
        if(!id.IsComplete())
            throw exception("EventAnalyzerDataId '%1%' is not complete.") % id;
        const std::string dir_name = id.GetName();
        const auto& sample_unc = GetModellingUncertainty(id);
        return std::make_shared<Data>(file, dir_name, channel, id, tuple, histNames, histDescs, sample_unc);
    }

    SampleUnc GetModellingUncertainty(const DataId& id) const
    {
        const std::string& sample = id.Get<std::string>();
        if(!unc_collection || !backgrounds.count(sample))
            return SampleUnc();
        const auto& unc_values = unc_collection->Get(channel, id.Get<EventCategory>());
        if(!unc_values.samples.count(sample))
            throw exception("Modelling uncertainty not found for '%1%'.") % sample;
        return unc_values.samples.at(sample);
    }

private:
    std::shared_ptr<TFile> file;
    Channel channel;
    const Tuple* tuple;
    NameSet histNames;
    DescSet histDescs;
    DataMap anaDataMap;
    NameSet backgrounds;
    MucPtr unc_collection;

};

} // namespace analysis
