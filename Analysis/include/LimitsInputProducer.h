/*! Definition of the limits input producer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "EventAnalyzerDataCollection.h"
#include "SampleDescriptor.h"

namespace analysis {

template<typename _FirstLeg, typename _SecondLeg>
class LimitsInputProducer {
public:
    using FirstLeg = _FirstLeg;
    using SecondLeg = _SecondLeg;
    using EventInfo = ::analysis::EventInfo<FirstLeg, SecondLeg>;
    using AnaData = ::analysis::EventAnalyzerData<FirstLeg, SecondLeg>;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection<AnaData>;
    using Sample = ::analysis::SampleDescriptorBase;
    using SampleCollection = std::vector<const Sample*>;
    using Hist = TH1D;
    using HistPtr = std::shared_ptr<root_ext::SmartHistogram<Hist>>;

    static constexpr Channel ChannelId() { return ChannelInfo::IdentifyChannel<FirstLeg, SecondLeg>(); }

    static std::string FullDataCardName(const std::string& datacard_name, EventEnergyScale eventEnergyScale)
    {
        if(eventEnergyScale == EventEnergyScale::Central)
            return datacard_name;

        std::ostringstream full_name;
        full_name << datacard_name << "_CMS_scale_";
        if(eventEnergyScale == EventEnergyScale::TauUp || eventEnergyScale == EventEnergyScale::TauDown)
            full_name << "t_";
        else if(eventEnergyScale == EventEnergyScale::JetUp || eventEnergyScale == EventEnergyScale::JetDown)
            full_name << "j";
        else
            throw exception("Unsupported event energy scale %1%.") % eventEnergyScale;
        full_name << "_13TeV";
        if(eventEnergyScale == EventEnergyScale::TauUp || eventEnergyScale == EventEnergyScale::JetUp)
            full_name << "Up";
        else
            full_name << "Down";
        return full_name.str();
    }

    LimitsInputProducer(AnaDataCollection& _anaDataCollection, const SampleDescriptorCollection& sample_descs,
                        const CombineSampleDescriptorCollection& combined_sample_descs) :
        anaDataCollection(&_anaDataCollection)
    {
        for(const auto& sample : sample_descs) {
            if(sample.datacard_name.size() || sample.datacard_name_ex.size())
                samples.push_back(&sample);
        }
        for(const auto& sample : combined_sample_descs) {
            if(sample.datacard_name.size() || sample.datacard_name_ex.size())
                samples.push_back(&sample);
        }
    }

    void ProduceLimitsInput(const std::string& outputFileNamePrefix, const std::string& hist_name,
                            const EventCategorySet& eventCategories, EventSubCategory eventSubCategory,
                            const EventEnergyScaleSet& eventEnergyScales)
    {
        static const std::map<EventCategory, std::string> categoryToDirectoryNameSuffix = {
            { EventCategory::TwoJets_ZeroBtag, "res0b" }, { EventCategory::TwoJets_OneBtag, "res1b" },
            { EventCategory::TwoJets_TwoBtag, "res2b" }
        };

        static const double tiny_value = 1e-9;
        static const double tiny_value_error = tiny_value;

        std::ostringstream s_file_name;
        s_file_name << outputFileNamePrefix << "_" << hist_name;
        if(eventSubCategory != EventSubCategory::NoCuts())
            s_file_name << "_" << eventSubCategory;
        s_file_name << ".root";
        const std::string file_name = s_file_name.str();

        auto outputFile = root_ext::CreateRootFile(file_name);
        for(EventCategory eventCategory : eventCategories) {
            if(!categoryToDirectoryNameSuffix.count(eventCategory)) continue;
            const std::string directoryName = ToString(ChannelId()) + "_"
                                            + categoryToDirectoryNameSuffix.at(eventCategory);
            TDirectory* directory = root_ext::GetDirectory(*outputFile, directoryName, true);
            for(const Sample* sample : samples) {
                for(const EventEnergyScale& eventEnergyScale : eventEnergyScales) {
                    std::vector<std::string> fullNames;
                    for(size_t n = 0; n < sample.GetNSignalPoints(); ++n)
                        fileNames.emplace_back(sample.GetFullName(n), sample.GetFileName(n));
                    if(!fullNames.size())
                        fullNames.push_back(sample->name);
                    const EventAnalyzerDataId dataId(eventCategory, eventSubCategory, eventEnergyScale,
                                                     EventRegion::SignalRegion(), sample->name)
                    const EventAnalyzerDataMetaId_noRegion_noName meta_id(eventCategory, eventSubCategory,
                                                                         eventEnergyScale);
                    std::shared_ptr<TH1D> hist;
                    if(auto hist_orig = GetSignalHistogram(meta_id, dataCategory->name, hist_name))
                        hist = std::shared_ptr<TH1D>(new TH1D(*hist_orig));
                    else {
                        std::cout << "Warning - Datacard histogram '" << hist_name
                                  << "' not found for data category '" << dataCategory->name << "' in '"
                                  << eventCategory << "/" << eventSubCategory << "/" << eventEnergyScale
                                  << "'. Using histogram with a tiny yield in the central bin instead.\n";

                        EventAnalyzerData& anaData = GetAnaData(meta_id.MakeId(EventRegion::OS_Isolated,
                                                                              dataCategory->name));
                        auto& new_hist = (anaData.*histogramAccessor)();
                        hist = std::shared_ptr<TH1D>(new TH1D(new_hist));
                        const Int_t central_bin = hist->GetNbinsX() / 2;
                        hist->SetBinContent(central_bin, tiny_value);
                        hist->SetBinError(central_bin, tiny_value_error);
                    }
                    const std::string full_datacard_name = FullDataCardName(sample->datacard_name, eventEnergyScale);
                    root_ext::WriteObject(*hist, directory, full_datacard_name);
                }
            }
        }
    }

private:
    AnaDataCollection* anaDataCollection;
    SampleCollection samples;
};

} // namespace analysis
