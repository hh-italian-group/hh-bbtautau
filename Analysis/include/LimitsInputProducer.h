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
    using SampleWP = Sample::Point;
    using SampleWPCollection = std::map<std::string, SampleWP>;
    using Hist = TH1D;
    using HistPtr = std::shared_ptr<root_ext::SmartHistogram<Hist>>;

    static constexpr Channel ChannelId() { return ChannelInfo::IdentifyChannel<FirstLeg, SecondLeg>(); }

    static std::string FullDataCardName(const std::string& datacard_name, EventEnergyScale es)
    {
        if(es == EventEnergyScale::Central)
            return datacard_name;

        std::ostringstream full_name;
        full_name << datacard_name << "_CMS_";
        if(es == EventEnergyScale::TauUp || es == EventEnergyScale::TauDown)
            full_name << "scale_t";
        else if(es == EventEnergyScale::JetUp || es == EventEnergyScale::JetDown)
            full_name << "scale_j";
        else if(es == EventEnergyScale::TopPtUp || es == EventEnergyScale::TopPtDown)
            full_name << "topPt";
        else
            throw exception("Unsupported event energy scale %1%.") % es;
        full_name << "_13TeV";
        if(es == EventEnergyScale::TauUp || es == EventEnergyScale::JetUp || es == EventEnergyScale::TopPtUp)
            full_name << "Up";
        else
            full_name << "Down";
        return full_name.str();
    }

    template<typename ...Args>
    LimitsInputProducer(AnaDataCollection& _anaDataCollection, Args&&... sample_descriptors) :
        anaDataCollection(&_anaDataCollection)
    {
        CollectWorkingPoints(std::forward<Args>(sample_descriptors)...);
    }

    void Produce(const std::string& outputFileNamePrefix, const std::string& hist_name,
                 const std::map<EventCategory, std::string>& eventCategories, EventSubCategory eventSubCategory,
                 const EventEnergyScaleSet& eventEnergyScales, const EventRegionSet& eventRegions)
    {
        static constexpr double tiny_value = 1e-9;
        static constexpr double tiny_value_error = tiny_value;
        static const std::string dirNamePrefix = ToString(ChannelId()) + "_";

        std::ostringstream s_file_name;
        s_file_name << outputFileNamePrefix << "_" << hist_name;
        if(eventSubCategory != EventSubCategory::NoCuts())
            s_file_name << "_" << eventSubCategory;
        s_file_name << ".root";
        const std::string file_name = s_file_name.str();
        auto outputFile = root_ext::CreateRootFile(file_name);

        std::set<EventAnalyzerDataId> empty_histograms;

        for(const EventAnalyzerDataId& metaId : EventAnalyzerDataId::MetaLoop(eventCategories, eventEnergyScales,
                                                                              sampleWorkingPoints, eventRegions))
        {
            std::string directoryName = dirNamePrefix + eventCategories.at(metaId.Get<EventCategory>());
            if(metaId.Get<EventRegion>() != EventRegion::SignalRegion())
                directoryName += "_" + ToString(metaId.Get<EventRegion>());
            TDirectory* directory = root_ext::GetDirectory(*outputFile, directoryName, true);
            const SampleWP& sampleWP = sampleWorkingPoints.at(metaId.Get<std::string>());
            const auto anaDataId = metaId.Set(eventSubCategory).Set(metaId.Get<EventRegion>());
            const auto& anaData = anaDataCollection->Get(anaDataId);
            auto& hist_entry = anaData.template GetEntryEx<TH1D>(hist_name);
            std::shared_ptr<TH1D> hist;
            if(hist_entry.GetHistograms().count(""))
                hist = std::make_shared<TH1D>(hist_entry());
            if(hist)
                hist->Scale(sampleWP.datacard_sf);
            if(!hist || hist->Integral() == 0.) {
                bool print_warning;
                if(CanHaveEmptyHistogram(anaDataId, print_warning)) continue;
                if(print_warning)
                    empty_histograms.insert(anaDataId);
                if(!hist)
                    hist = std::make_shared<TH1D>(hist_entry());
                const Int_t central_bin = hist->GetNbinsX() / 2;
                hist->SetBinContent(central_bin, tiny_value);
                hist->SetBinError(central_bin, tiny_value_error);
            }
            const auto datacard_name = FullDataCardName(sampleWP.datacard_name, metaId.Get<EventEnergyScale>());
            root_ext::WriteObject(*hist, directory, datacard_name);
        }

        if(empty_histograms.size()) {
            std::cout << "\t\tWarning: following datacard histograms are empty:\n";
            for(const auto& id : empty_histograms)
                std::cout << "\t\t\t" << id << "\n";
            std::cout << "\t\tUsing histograms with a tiny yield in the central bin instead." << std::endl;
        }
    }

private:
    void CollectWorkingPoints() {}

    template<typename SampleCollection, typename ...Args>
    void CollectWorkingPoints(const SampleCollection& sample_descriptors, Args&&... other_sample_descriptors)
    {
        for(const auto& sample : sample_descriptors) {
            for(const SampleWP& wp : sample.second.working_points) {
                if(wp.datacard_name.size())
                    sampleWorkingPoints[wp.full_name] = wp;
            }
        }
        CollectWorkingPoints(std::forward<Args>(other_sample_descriptors)...);
    }

    bool CanHaveEmptyHistogram(const EventAnalyzerDataId& id, bool& print_warning) const
    {
        const auto& es = id.Get<EventEnergyScale>();
        const SampleWP& sampleWP = sampleWorkingPoints.at(id.Get<std::string>());
        print_warning = id.Get<EventRegion>() == EventRegion::SignalRegion();
        if(es == EventEnergyScale::Central)
            return false;
        if(sampleWP.sampleType == SampleType::Data || sampleWP.sampleType == SampleType::QCD)
            return true;
        if(es == EventEnergyScale::TopPtUp || es == EventEnergyScale::TopPtDown)
            return sampleWP.sampleType != SampleType::TT;
        return false;
    }

private:
    AnaDataCollection* anaDataCollection;
    SampleWPCollection sampleWorkingPoints;
};

} // namespace analysis
