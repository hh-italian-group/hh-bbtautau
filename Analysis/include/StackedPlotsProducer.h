/*! Definition of BaseEventAnalyzer class, the base class for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

namespace analysis {

template<typename _FirstLeg, typename _SecondLeg>
class BaseEventAnalyzer {
public:
    using FirstLeg = _FirstLeg;
    using SecondLeg = _SecondLeg;
    using EventInfo = ::analysis::EventInfo<FirstLeg, SecondLeg>;
    using EventAnalyzerData = ::analysis::EventAnalyzerData<FirstLeg, SecondLeg>;
    using PhysicalValueMap = std::map<EventRegion, PhysicalValue>;

    void PrintStackedPlots(EventRegion eventRegion, bool isBlind, bool drawRatio)
    {
        const std::string blindCondition = isBlind ? "_blind" : "_noBlind";
        const std::string ratioCondition = drawRatio ? "_ratio" : "_noRatio";
        std::ostringstream eventRegionName;
        eventRegionName << args.outputFileName() << blindCondition << ratioCondition << "_" << eventRegion << ".pdf";
        root_ext::PdfPrinter printer(eventRegionName.str());

        for(EventCategory eventCategory : EventCategoriesToProcess()) {
            for (const auto& hist_name : EventAnalyzerData::template GetOriginalHistogramNames<TH1D>()) {
                for(EventSubCategory subCategory : EventSubCategoriesToProcess()) {
                    const EventAnalyzerDataMetaId_noRegion_noName anaDataMetaId(eventCategory, subCategory,
                                                                               EventEnergyScale::Central);
                    std::ostringstream ss_title;
                    ss_title << eventCategory;
                    if(subCategory != EventSubCategory::NoCuts)
                        ss_title << " " << subCategory;
                    ss_title << ": " << hist_name;

                    StackedPlotDescriptor stackDescriptor(ss_title.str(), false, ChannelNameLatex(),
                                                          __EventCategory_names<>::names.EnumToString(eventCategory), drawRatio,
                                                          false);

                    for(const DataCategory* category : dataCategoryCollection.GetAllCategories()) {
                        if(!category->draw) continue;

                        const auto histogram = GetHistogram(anaDataMetaId, eventRegion, category->name, hist_name);
                        if(!histogram) continue;

                        if(category->IsSignal() && eventCategory == EventCategory::TwoJets_Inclusive) continue;
                        else if(category->IsSignal())
                            stackDescriptor.AddSignalHistogram(*histogram, category->title, category->color,
                                                               category->draw_sf);
                        else if(category->IsBackground())
                            stackDescriptor.AddBackgroundHistogram(*histogram, category->title, category->color);
                        else if(category->IsData())
                            stackDescriptor.AddDataHistogram(*histogram, category->title, isBlind,
                                                             GetBlindRegion(subCategory, hist_name));
                    }

                    printer.PrintStack(stackDescriptor);
                }
            }
        }
    }

    static const std::vector< std::pair<double, double> >& GetBlindRegion(EventSubCategory subCategory,
                                                                          const std::string& hist_name)
    {
        static const std::vector< std::vector< std::pair<double, double> > > blindingRegions = {
            { { std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest() } },
            { { 100, 150 } },
            { { 200, 400 } },
            { { 100, 150 }, { 450, 500 }, { 800, 850 }, { 1150, 1200 }, { 1500, 1550 } }
        };

        static const std::map<std::string, size_t> histogramsToBlind = {
            { EventAnalyzerData::m_sv_Name(), 1 }, { EventAnalyzerData::m_vis_Name(), 1 },
            { EventAnalyzerData::m_bb_Name(), 1 }, { EventAnalyzerData::m_ttbb_Name(), 2 },
            { EventAnalyzerData::m_ttbb_kinfit_Name(), 2 }
        };

};

} // namespace analysis
