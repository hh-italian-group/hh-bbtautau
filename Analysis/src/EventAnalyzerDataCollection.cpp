/*! Collection of histogram containers for event analyzers.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "hh-bbtautau/Analysis/include/EventAnalyzerDataCollection.h"

namespace analysis {


EventAnalyzerDataCollection::EventAnalyzerDataCollection(std::shared_ptr<TFile> _file, Channel _channel,
                                                         const NameSet& _histNames,
                                                         const DescSet& _histDescs, bool _readMode,
                                                         const NameSet& _backgrounds, MucPtr _unc_collection) :
    file(_file), channel(_channel), histNames(_histNames), histDescs(_histDescs), readMode(_readMode),
    backgrounds(_backgrounds), unc_collection(_unc_collection)
{
}

EventAnalyzerDataCollection::Data& EventAnalyzerDataCollection::Get(const DataId& id)
{
    std::lock_guard<Mutex> lock(mutex);
    auto& anaData = anaDataMap[id];
    if(!anaData)
        anaData = Make(id);
    return *anaData;
}

const EventAnalyzerDataCollection::DataMap& EventAnalyzerDataCollection::GetAll() const { return anaDataMap; }
Channel EventAnalyzerDataCollection::ChannelId() const { return channel; }
bool EventAnalyzerDataCollection::ReadMode() const { return readMode; }

EventAnalyzerDataCollection::DataPtr EventAnalyzerDataCollection::Make(const DataId& id) const
{
    if(!id.IsComplete())
        throw exception("EventAnalyzerDataId '%1%' is not complete.") % id;
    const std::string dir_name = id.GetName();
    const auto& sample_unc = GetModellingUncertainty(id);
    return std::make_shared<Data>(file, dir_name, channel, id, histNames, histDescs, sample_unc, readMode);
}

EventAnalyzerDataCollection::SampleUnc EventAnalyzerDataCollection::GetModellingUncertainty(const DataId& id) const
{
    const std::string& sample = id.Get<std::string>();
    if(!unc_collection || !backgrounds.count(sample))
        return SampleUnc();
    const auto& unc_values = unc_collection->Get(channel, id.Get<EventCategory>());
    if(!unc_values.samples.count(sample))
        throw exception("Modelling uncertainty not found for '%1%'.") % sample;
    return unc_values.samples.at(sample);
}

} // namespace analysis
