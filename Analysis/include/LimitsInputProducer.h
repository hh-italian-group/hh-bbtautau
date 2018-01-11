/*! Definition of the limits input producer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "EventAnalyzerDataCollection.h"
#include "SampleDescriptor.h"

namespace analysis {

class LimitsInputProducer {
public:
    using AnaData = ::analysis::EventAnalyzerData;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection;
    using Sample = ::analysis::SampleDescriptorBase;
    using SampleWP = Sample::Point;
    using SampleWPCollection = std::map<std::string, SampleWP>;
    using Hist = TH1D;
    using HistPtr = std::shared_ptr<root_ext::SmartHistogram<Hist>>;

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

    static std::string EventRegionSuffix(EventRegion region)
    {
        static const std::map<EventRegion, std::string> regions {
            { EventRegion::OS_AntiIsolated(), "_OS_antiiso" },
            { EventRegion::SS_AntiIsolated(), "_SS_antiiso" },
            { EventRegion::SS_Isolated(), "_SS_iso" },
        };

        if(region == EventRegion::SignalRegion())
            return "";
        if(regions.count(region))
            return regions.at(region);
        return "_" + ToString(region);
    }

    template<typename ...Args>
    LimitsInputProducer(AnaDataCollection& _anaDataCollection, Args&&... sample_descriptors) :
        anaDataCollection(&_anaDataCollection)
    {
        CollectWorkingPoints(std::forward<Args>(sample_descriptors)...);
    }

    void Produce(const std::string& outputFileNamePrefix, const std::string& hist_name,
                 const std::map<EventCategory, std::string>& eventCategories, EventSubCategory eventSubCategory,
                 const EventEnergyScaleSet& eventEnergyScales, const EventRegionSet& eventRegions,
                 const std::map<SelectionCut, std::string>& sel_aliases)
    {
        static constexpr double tiny_value = 1e-9;
        static constexpr double tiny_value_error = tiny_value;
        static const std::string dirNamePrefix = ToString(anaDataCollection->ChannelId()) + "_";

        std::ostringstream s_file_name;
        s_file_name << outputFileNamePrefix << "_" << hist_name;
        if(eventSubCategory != EventSubCategory::NoCuts())
            s_file_name << "_" << eventSubCategory.ToString(sel_aliases);
        const std::string file_name = s_file_name.str();
        auto outputFile = root_ext::CreateRootFile(file_name + ".root");
        std::set<EventAnalyzerDataId> empty_histograms;

        for(const EventAnalyzerDataId& metaId : EventAnalyzerDataId::MetaLoop(eventCategories, eventEnergyScales,
                                                                              sampleWorkingPoints, eventRegions))
        {
            const std::string directoryName = dirNamePrefix + eventCategories.at(metaId.Get<EventCategory>())
                    + EventRegionSuffix(metaId.Get<EventRegion>());
            TDirectory* directory = root_ext::GetDirectory(*outputFile, directoryName, true);
            const SampleWP& sampleWP = sampleWorkingPoints.at(metaId.Get<std::string>());
            const auto anaDataId = metaId.Set(eventSubCategory).Set(metaId.Get<EventRegion>());
            const auto& anaData = anaDataCollection->Get(anaDataId);
            auto& hist_entry = anaData.GetEntryEx<TH1D>(hist_name);
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
            const std::string of_name = file_name + "_emptyShapes.txt";
            std::ofstream of(of_name);
            of.exceptions(std::ios::failbit);
            for(const auto& id : empty_histograms)
                of << id << "\n";
            std::cout << "\t\t\tWarning: some datacard histograms are empty.\n"
                      << "\t\t\tThey are replaced with histograms with a tiny yield in the central bin.\n"
                      << "\t\t\tSee '" << of_name << "' for details." << std::endl;
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
