/*! Definition of the stacked prots producer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "EventAnalyzerDataCollection.h"
#include "SampleDescriptor.h"
#include "AnalysisTools/Print/include/RootPrintToPdf.h"

namespace analysis {

class StackedPlotsProducer {
public:
    using AnaData = ::analysis::EventAnalyzerData;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection;
    using Sample = ::analysis::SampleDescriptorBase;
    using SampleCollection = std::vector<const Sample*>;
    using Hist = TH1D;
    using HistPtr = std::shared_ptr<root_ext::SmartHistogram<Hist>>;

    static SampleCollection CreateOrderedSampleCollection(const std::vector<std::string>& draw_sequence,
                                                          const SampleDescriptorCollection& samples,
                                                          const CombinedSampleDescriptorCollection& combined_samples,
                                                          const std::vector<std::string>& signals,
                                                          const std::vector<std::string>& data, Channel channel)
    {
        SampleCollection selected_samples;
        for(const auto& name : draw_sequence) {
            if(samples.count(name)) {
                selected_samples.push_back(&samples.at(name));
            } else if(combined_samples.count(name)) {
                selected_samples.push_back(&combined_samples.at(name));
            } else if(name == "signals") {
                for(const auto& signal_name : signals) {
                    if(!samples.count(signal_name))
                        throw exception("Signal sample '%1%' not found while creating drawing sequence.") % signal_name;
                    selected_samples.push_back(&samples.at(signal_name));
                }
            } else if(name == "data") {
                for(const auto& data_name : data) {
                    if(!samples.count(data_name))
                        throw exception("Data sample '%1%' not found while creating drawing sequence.") % data_name;
                    selected_samples.push_back(&samples.at(data_name));
                }
            } else
                throw exception("Sample '%1%' not found while creating drawing sequence.") % name;
        }
        SampleCollection result;
        for(const Sample* sample : selected_samples) {
            if(!sample->channels.size() || sample->channels.count(channel))
                result.push_back(sample);
        }
        return result;
    }


    StackedPlotsProducer(AnaDataCollection& _anaDataCollection, const SampleCollection& _samples,
                         bool _isBlind, bool _drawRatio, const std::set<std::string>& _histogramNames = {}) :
        anaDataCollection(&_anaDataCollection), samples(_samples), histogramNames(_histogramNames), isBlind(_isBlind),
        drawRatio(_drawRatio)
    {
        if(!histogramNames.size()) {
            for(const auto& anaData : anaDataCollection->GetAll()) {
                for(const auto& item : anaData.second->template GetHistogramsEx<TH1D>()) {
                    histogramNames.insert(item.first);
                }
            }
        }
    }

    Channel ChannelId() const { return anaDataCollection->ChannelId(); }
    const std::string& ChannelNameLatex() const { return __Channel_names_latex.EnumToString(ChannelId()); }

    void PrintStackedPlots(const std::string& outputFileNamePrefix, const EventRegion& eventRegion,
                           const EventCategorySet& eventCategories,
                           const EventSubCategorySet& eventSubCategories, const std::set<std::string>& signals) const
    {
        const std::string blindCondition = isBlind ? "_blind" : "_noBlind";
        const std::string ratioCondition = drawRatio ? "_ratio" : "_noRatio";
        std::ostringstream outputFileName;
        outputFileName << outputFileNamePrefix << blindCondition << ratioCondition << "_" << eventRegion << ".pdf";
        root_ext::PdfPrinter printer(outputFileName.str());

        for(const EventCategory& eventCategory : eventCategories) {
            for (const auto& hist_name : histogramNames) {
                for(const EventSubCategory& subCategory : eventSubCategories) {
                    const EventAnalyzerDataId anaDataMetaId(eventRegion, eventCategory, subCategory,
                                                            EventEnergyScale::Central);
                    std::ostringstream ss_title;
                    ss_title << eventCategory;
                    if(subCategory != EventSubCategory::NoCuts())
                        ss_title << " " << subCategory;
                    ss_title << ": " << hist_name;

                    StackedPlotDescriptor stackDescriptor(ss_title.str(), false, ChannelNameLatex(),
                                                          ToString(eventCategory), drawRatio, false);

                    for(const Sample* sample : samples) {
                        for(const Sample::Point& item : sample->working_points) {
                            if(!item.draw) continue;
                            const std::string& item_name = item.full_name;
                            const root_ext::Color& color = item.color;
                            const auto histogram = GetHistogram(anaDataMetaId, item_name, hist_name);
                            if(!histogram || histogram->Integral() == 0.) continue;

                            if(signals.count(sample->name))
                                stackDescriptor.AddSignalHistogram(*histogram, item.title, color, sample->draw_sf);
                            else if(sample->sampleType == SampleType::Data)
                                stackDescriptor.AddDataHistogram(*histogram, item.title, isBlind,
                                                                 GetBlindRegion(subCategory, hist_name));
                            else
                                stackDescriptor.AddBackgroundHistogram(*histogram, item.title, color);
                        }
                    }

                    printer.PrintStack(stackDescriptor);
                }
            }
        }
    }

private:
    HistPtr GetHistogram(const EventAnalyzerDataId& metaId, const std::string& sample_name,
                         const std::string& hist_name) const
    {
        const EventAnalyzerDataId dataId = metaId.Set(sample_name);
        const auto& anaData = anaDataCollection->Get(dataId);
        if(anaData.ReadMode()) {
            auto& entry = anaData.template GetEntryEx<Hist>(hist_name);
            entry.Read();
        }
        return anaData.template TryGetHistogramEx<Hist>(hist_name);
    }

    static const std::vector< std::pair<double, double> >& GetBlindRegion(const EventSubCategory& /*subCategory*/,
                                                                          const std::string& /*hist_name*/)
    {
        static const std::vector< std::vector< std::pair<double, double> > > blindingRegions = {
            { { std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest() } },
        };

        return blindingRegions.at(0);
    }

private:
    AnaDataCollection* anaDataCollection;
    SampleCollection samples;
    std::set<std::string> histogramNames;
    bool isBlind, drawRatio;
};

} // namespace analysis
