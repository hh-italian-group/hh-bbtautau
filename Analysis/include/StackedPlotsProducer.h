/*! Definition of the stacked prots producer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "EventAnalyzerDataCollection.h"
#include "SampleDescriptor.h"
#include "AnalysisTools/Print/include/PdfPrinter.h"
#include "AnalysisTools/Print/include/StackedPlotDescriptor.h"

namespace analysis {

class StackedPlotsProducer {
public:
    using AnaData = ::analysis::EventAnalyzerData;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection;
    using Sample = ::analysis::SampleDescriptorBase;
    using SampleCollection = std::vector<const Sample*>;
    using Hist = TH1D;
    using HistPtr = std::shared_ptr<root_ext::SmartHistogram<Hist>>;
    using PlotConfigReader = ::analysis::PropertyConfigReader;
    using PlotConfig = PlotConfigReader::ItemCollection;
    using PageOptions = ::root_ext::draw_options::Page;
    using StackedPlotDescriptor = ::root_ext::StackedPlotDescriptor;

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
                         const std::string& plot_cfg_name, const std::string& page_opt_name,
                         const std::set<std::string>& _histogramNames = {}) :
        anaDataCollection(&_anaDataCollection), samples(_samples), histogramNames(_histogramNames)
    {
        if(!histogramNames.size()) {
            for(const auto& anaData : anaDataCollection->GetAll()) {
                for(const auto& item : anaData.second->template GetHistogramsEx<TH1D>()) {
                    histogramNames.insert(item.first);
                }
            }
        }

        PlotConfigReader cfg_reader;
        cfg_reader.Parse(plot_cfg_name);
        plot_cfg = cfg_reader.GetItems();
        if(!plot_cfg.count(page_opt_name)) {
            throw exception("Page drawing options = '%1%' not found in config file '%2%'.")
                % page_opt_name % plot_cfg_name;
        }
        page_opt = PageOptions(plot_cfg.at(page_opt_name));
    }

    Channel ChannelId() const { return anaDataCollection->ChannelId(); }
    const std::string& ChannelNameLatex() const { return __Channel_names_latex.EnumToString(ChannelId()); }

    void PrintStackedPlots(const std::string& outputFileNamePrefix, const EventRegion& eventRegion,
                           const EventCategorySet& eventCategories,
                           const EventSubCategorySet& eventSubCategories, const std::set<std::string>& signals)
    {
        std::ostringstream outputFileName;
        outputFileName << outputFileNamePrefix << "_" << eventRegion << ".pdf";
        root_ext::PdfPrinter printer(outputFileName.str(), plot_cfg, page_opt);

        std::vector<std::pair<std::string, StackedPlotDescriptor>> descriptors;
        for(auto category_iter = eventCategories.begin(); category_iter != eventCategories.end(); ++category_iter) {
            const auto& eventCategory = *category_iter;
            for (auto hist_name_iter = histogramNames.begin(); hist_name_iter != histogramNames.end();
                 ++hist_name_iter) {
                const auto& hist_name = *hist_name_iter;
                for(auto sub_category_iter = eventSubCategories.begin(); sub_category_iter != eventSubCategories.end();
                    ++sub_category_iter) {
                    const auto& subCategory = *sub_category_iter;
                    const EventAnalyzerDataId anaDataMetaId(eventRegion, eventCategory, subCategory,
                                                            EventEnergyScale::Central);
                    std::ostringstream ss_title;
                    ss_title << eventCategory;
                    if(subCategory != EventSubCategory::NoCuts())
                        ss_title << " " << subCategory;
                    ss_title << ": " << hist_name;

                    if(plot_cfg.count("cat_text")) {
                        const auto cat_text = ChannelNameLatex() + ", " + ToString(eventCategory);
                        printer.GetLabelOptions("cat_text").SetText(cat_text);
                    }

                    StackedPlotDescriptor stackDescriptor(page_opt, plot_cfg);

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
                                stackDescriptor.AddDataHistogram(*histogram, item.title);
                            else
                                stackDescriptor.AddBackgroundHistogram(*histogram, item.title, color);
                        }
                    }

                    const bool is_last = std::next(category_iter) == eventCategories.end()
                            && std::next(hist_name_iter) == histogramNames.end()
                            && std::next(sub_category_iter) == eventSubCategories.end();
                    printer.Print(ss_title.str(), stackDescriptor, is_last);
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
    PlotConfig plot_cfg;
    PageOptions page_opt;
};

} // namespace analysis
