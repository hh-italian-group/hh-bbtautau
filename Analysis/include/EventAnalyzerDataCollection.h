/*! Collection of histogram containers for event analyzers.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

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
            MucPtr _unc_collection = MucPtr());

    Data& Get(const DataId& id);
    const DataMap& GetAll() const;
    Channel ChannelId() const;
    bool ReadMode() const;
    void Fill(const DataId& id, double weight);

private:
    DataPtr Make(const DataId& id) const;
    SampleUnc GetModellingUncertainty(const DataId& id) const;

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
