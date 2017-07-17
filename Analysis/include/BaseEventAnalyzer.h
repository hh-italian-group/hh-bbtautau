/*! Definition of BaseEventAnalyzer class, the base class for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <iostream>
#include <cmath>
#include <set>
#include <list>
#include <locale>

#include <TColor.h>
#include <TLorentzVector.h>

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Print/include/RootPrintToPdf.h"

#include "AnalysisCategories.h"
#include "EventAnalyzerDataCollection.h"

#include "h-tautau/Analysis/include/Htautau_2015.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/McCorrections/include/EventWeights.h"

namespace analysis {

struct AnalyzerArguments {
    REQ_ARG(std::string, source_cfg);
    REQ_ARG(std::string, inputPath);
    REQ_ARG(std::string, outputFileName);
    REQ_ARG(std::string, signal_list);
    OPT_ARG(bool, saveFullOutput, false);
};

template<typename _FirstLeg, typename _SecondLeg>
class BaseEventAnalyzer {
public:
    using FirstLeg = _FirstLeg;
    using SecondLeg = _SecondLeg;
    using EventInfo = ::analysis::EventInfo<FirstLeg, SecondLeg>;
    using EventAnalyzerData = ::analysis::EventAnalyzerData<FirstLeg, SecondLeg>;
    using PhysicalValueMap = std::map<EventRegion, PhysicalValue>;

    static constexpr Channel ChannelId() { return ChannelInfo::IdentifyChannel<FirstLeg, SecondLeg>(); }

    virtual const EventCategorySet& EventCategoriesToProcess() const
    {
        static const EventCategorySet categories = {
            EventCategory::TwoJets_Inclusive(), EventCategory::TwoJets_ZeroBtag(),
            EventCategory::TwoJets_OneBtag(), EventCategory::TwoJets_OneLooseBtag(),
            EventCategory::TwoJets_TwoBtag(), EventCategory::TwoJets_TwoLooseBtag()
        };
        return categories;
    }

    virtual const EventSubCategorySet& EventSubCategoriesToProcess() const
    {
        static const EventSubCategorySet sub_categories = {
            EventSubCategory().SetCutResult(SelectionCut::InsideMassWindow, true)
                              .SetCutResult(SelectionCut::MVA, true)
                              .SetCutResult(SelectionCut::KinematicFitConverged, true)
        };
        return sub_categories;
    }

    virtual const EventRegionSet& EventRegionsToProcess() const
    {
        static const EventRegionSet regions = {
            EventRegion::OS_Isolated(), EventRegion::OS_AntiIsolated(),
            EventRegion::SS_Isolated(), EventRegion::SS_AntiIsolated()
        };
        return regions;
    }

    virtual const EventEnergyScaleSet& EventEnergyScaleToProcess() const
    {
        static const EventEnergyScaleSet scales = {
            EventEnergyScale::Central, EventEnergyScale::TauUp, EventEnergyScale::TauDown,
            EventEnergyScale::JetUp, EventEnergyScale::JetDown
        };
        return scales;
    }

    EventCategorySet DetermineEventCategories(const  std::vector<float>& csv_Bjets,
                                              const EventInfoBase::JetPair& selected_bjets,
                                              double CSVL, double CSVM)
    {
        EventCategorySet categories;
        categories.insert(EventCategory::Inclusive());

        std::map<DiscriminatorWP, size_t> n_bjets;

        static const std::map< size_t, EventCategory> mediumCategories_map {
            {  0 , EventCategory::TwoJets_ZeroBtag }, { 1 , EventCategory::TwoJets_OneBtag },
            {  2 , EventCategory::TwoJets_TwoBtag }
        };

        static const std::map< size_t, EventCategory> looseCategories_map {
            { 0 , EventCategory::TwoJets_ZeroLooseBtag }, { 1, EventCategory::TwoJets_OneLooseBtag },
            { 2 , EventCategory::TwoJets_TwoLooseBtag }
        };

        if (selected_bjets.first < csv_Bjets.size() && selected_bjets.second < csv_Bjets.size()){
            categories.push_back(EventCategory::TwoJets_Inclusive);

            size_t n_mediumBtag = 0;
            if(doRetag) {
                n_mediumBtag = std::min<size_t>(nBjets_retagged, 2);
            } else {
                if(csv_Bjets.at(selected_bjets.first) > CSVM) ++n_mediumBtag;
                if(csv_Bjets.at(selected_bjets.second) > CSVM) ++n_mediumBtag;
            }

            if(mediumCategories_map.count(n_mediumBtag))
                categories.push_back(mediumCategories_map.at(n_mediumBtag));
            if(n_mediumBtag > 0)
                categories.push_back(EventCategory::TwoJets_AtLeastOneBtag);

            size_t n_looseBtag = 0;
            if(csv_Bjets.at(selected_bjets.first) > CSVL) ++n_looseBtag;
            if(csv_Bjets.at(selected_bjets.second) > CSVL) ++n_looseBtag;

            if(looseCategories_map.count(n_looseBtag))
                categories.push_back(looseCategories_map.at(n_looseBtag));
            if(n_looseBtag > 0)
                categories.push_back(EventCategory::TwoJets_AtLeastOneLooseBtag);
        }

        return categories;
    }


    BaseEventAnalyzer(const AnalyzerArguments& _args)
        : args(_args), dataCategoryCollection(args.source_cfg(), args.signal_list(), ChannelId()),
          anaDataCollection(args.outputFileName() + "_full.root", args.saveFullOutput()),
          weights(Period::Run2015, DiscriminatorWP::Medium)
    {
    }

    void Run()
    {
        std::cout << "Processing data categories... " << std::endl;
        for(const DataCategory* dataCategory : dataCategoryCollection.GetAllCategories()) {
            if(!dataCategory->sources_sf.size()) continue;
            std::cout << *dataCategory << "   isData: "<<dataCategory->IsData()<<std::endl;
            for(const auto& source_entry : dataCategory->sources_sf) {
                const std::string fullFileName = args.inputPath() + "/" + source_entry.first;
                auto file = root_ext::OpenRootFile(fullFileName);
                std::shared_ptr<ntuple::EventTuple> tree(new ntuple::EventTuple(TreeName(), file.get(), true,
                                                                            disabled_branches));
                ProcessDataSource(*dataCategory, tree, source_entry.second);
            }
        }

        std::set<std::string> histograms_to_report = { EventAnalyzerData::m_ttbb_kinfit_Name() };

        for (const auto& hist_name : EventAnalyzerData::template GetOriginalHistogramNames<TH1D>()) {
            for(auto subCategory : EventSubCategoriesToProcess()) {
                for(auto energyScale : EventEnergyScaleToProcess()) {
                    std::cout << "Processing '" << hist_name << "' in " << subCategory << "/" << energyScale
                              << "..." << std::endl;

                    std::ostringstream debug;
                    std::ostream& s_out = histograms_to_report.count(hist_name) ? std::cout : debug;

                    for (EventCategory eventCategory : EventCategoriesToProcess()) {
                        const EventAnalyzerDataMetaId_noRegion_noName anaDataMetaId(eventCategory, subCategory,
                                                                                   energyScale);
                        if(dataCategoryCollection.GetCategories(DataCategoryType::Data).size()) {
                            DataCategoryType dataCategoryType = DataCategoryType::QCD;
                            const auto qcd_yield = CalculateQCDYield(anaDataMetaId, hist_name, dataCategoryType, s_out);
                            s_out << eventCategory << ": QCD yield = " << qcd_yield << ".\n";
                            EstimateQCD(anaDataMetaId, hist_name, qcd_yield, dataCategoryType);
                        }
                        ProcessCompositDataCategories(anaDataMetaId, hist_name);
                    }
                }
            }
        }

        std::cout << "\nSaving tables... " << std::endl;
        PrintTables("comma", L",");

        std::cout << "Saving datacards... " << std::endl;
        ProduceFileForLimitsCalculation(EventAnalyzerData::m_ttbb_kinfit_Name(),
                                        EventSubCategory::KinematicFitConvergedWithMassWindow,
                                        &EventAnalyzerData::m_ttbb_kinfit);

        std::cout << "Printing stacked plots... " << std::endl;
        PrintStackedPlots(EventRegion::OS_Isolated, false, true);
        PrintStackedPlots(EventRegion::OS_Isolated, false, false);
        std::cout << "Saving output file..." << std::endl;
    }

protected:

    virtual std::string TreeName() const = 0;
    virtual EventRegion DetermineEventRegion(EventInfo& event, EventCategory eventCategory) = 0;

    static EventSubCategorySet DetermineEventSubCategories(EventInfo& event)
    {
        using namespace cuts::massWindow;

        EventSubCategorySet sub_categories;
        sub_categories.insert(EventSubCategory::NoCuts);
        if(event.HasBjetPair()) {
            const double mass_tautau = event.GetHiggsTTMomentum(true).M();
            const double mass_bb = event.GetHiggsBB().GetMomentum().M();

            if(event.GetKinFitResults().HasValidMass())
                sub_categories.insert(EventSubCategory::KinematicFitConverged);

            if(mass_tautau > m_tautau_low && mass_tautau < m_tautau_high
                    && mass_bb > m_bb_low && mass_bb < m_bb_high) {
                sub_categories.insert(EventSubCategory::MassWindow);
                if(event.GetKinFitResults().HasValidMass())
                    sub_categories.insert(EventSubCategory::KinematicFitConvergedWithMassWindow);
            } else {
                sub_categories.insert(EventSubCategory::OutsideMassWindow);
                if(event.GetKinFitResults().HasValidMass())
                    sub_categories.insert(EventSubCategory::KinematicFitConvergedOutsideMassWindow);
            }
        }
        return sub_categories;
    }

    virtual PhysicalValue CalculateQCDYield(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                            const std::string& hist_name, DataCategoryType dataCategoryType,
                                            std::ostream& s_out) = 0;
    virtual void EstimateQCD(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId, const std::string& hist_name,
                             const PhysicalValue& scale_factor, DataCategoryType dataCategoryType) = 0;
    virtual void CreateHistogramForZTT(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                       const std::string& hist_name, const PhysicalValueMap& ztt_yield,
                                       bool useEmbedded) = 0;
    virtual void CreateHistogramForVVcategory(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                              const std::string& hist_name) = 0;

    PhysicalValueMap CalculateZTTmatchedYield(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                              const std::string& hist_name, bool useEmbedded)
    {
        const DataCategory& ZTT_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_MC);

        PhysicalValueMap zttYield;
        if (useEmbedded){
            const DataCategory& embedded =  dataCategoryCollection.GetUniqueCategory(DataCategoryType::Embedded);
            const DataCategory& TTembedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::TT_Embedded);
            const EventAnalyzerDataMetaId_noRegion_noName anaDataMetaId_incl(EventCategory::Inclusive,
                                                                            anaDataMetaId.eventSubCategory,
                                                                            anaDataMetaId.eventEnergyScale);

            TH1D* hist_embedded_inclusive = GetSignalHistogram(anaDataMetaId_incl, embedded.name, hist_name);
            TH1D* hist_embedded_category = GetSignalHistogram(anaDataMetaId, embedded.name, hist_name);
            TH1D* hist_TTembedded_inclusive = GetSignalHistogram(anaDataMetaId_incl, TTembedded.name, hist_name);
            TH1D* hist_TTembedded_category = GetSignalHistogram(anaDataMetaId, TTembedded.name, hist_name);
            TH1D* hist_ztautau_inclusive = GetSignalHistogram(anaDataMetaId_incl, ZTT_MC.name, hist_name);
            if(!hist_embedded_inclusive)
                throw exception("embedded hist in inclusive category not found");
            if(!hist_embedded_category)
                throw exception("embedded hist not found in event category: %1%\n") % anaDataMetaId.eventCategory;
            if(!hist_TTembedded_inclusive)
                throw exception("TTembedded hist in inclusive category not found");
            if(!hist_TTembedded_category)
                throw exception("TTembedded hist not found in event category: %1%\n") % anaDataMetaId.eventCategory;
            if(!hist_ztautau_inclusive )
                throw exception("ztt hist in inclusive category not found");

            const PhysicalValue n_emb_inclusive =
                    Integral(*hist_embedded_inclusive, true) - Integral(*hist_TTembedded_inclusive, true);
            const PhysicalValue n_emb_category =
                    Integral(*hist_embedded_category, true) - Integral(*hist_TTembedded_category, true);
            const PhysicalValue n_ztautau_inclusive = Integral(*hist_ztautau_inclusive, true);
            const PhysicalValue embedded_eff = n_emb_category/n_emb_inclusive;
            zttYield[EventRegion::OS_Isolated] = n_ztautau_inclusive * embedded_eff;
        }

        for (EventRegion eventRegion : AllEventRegions){
            if (zttYield.count(eventRegion)) continue;
            TH1D* hist_ztautau = GetHistogram(anaDataMetaId, eventRegion, ZTT_MC.name, hist_name);

            if(!hist_ztautau ){
                if (eventRegion == EventRegion::OS_Isolated)
                    throw exception("ztt hist not found in event category %1%.") % anaDataMetaId.eventCategory;
                continue;
            }

            zttYield[eventRegion] = Integral(*hist_ztautau, true);
        }
        return zttYield;
    }

    virtual void CreateHistogramForZcategory(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                             const std::string& hist_name)
    {
        const std::map<DataCategoryType, DataCategoryType> z_type_category_map = {
            { DataCategoryType::ZL_MC, DataCategoryType::ZL }, { DataCategoryType::ZJ_MC, DataCategoryType::ZJ }
        };

        for (const auto& z_category : z_type_category_map){
            const DataCategory& originalZcategory = dataCategoryCollection.GetUniqueCategory(z_category.first);
            const DataCategory& newZcategory = dataCategoryCollection.GetUniqueCategory(z_category.second);

            PhysicalValueMap valueMap;

            for(EventRegion eventRegion : AllEventRegions) {
                auto z_hist_yield = GetHistogram(anaDataMetaId, eventRegion, originalZcategory.name, hist_name);
                if (z_hist_yield)
                    valueMap[eventRegion] = Integral(*z_hist_yield,true);
            }

            static const EventCategorySet categoriesToRelax = {
                EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_AtLeastOneBtag };

            for(const auto& yield_iter : valueMap) {
                const EventRegion eventRegion = yield_iter.first;
                const PhysicalValue& yield = yield_iter.second;
                const EventCategory shapeEventCategory = (categoriesToRelax.count(anaDataMetaId.eventCategory) &&
                        eventRegion == analysis::EventRegion::OS_Isolated)
                        ? MediumToLoose_EventCategoryMap.at(anaDataMetaId.eventCategory) : anaDataMetaId.eventCategory;
                const EventAnalyzerDataMetaId_noRegion_noName anaDataMetaId_shape(shapeEventCategory,
                                                                                 anaDataMetaId.eventSubCategory,
                                                                                 anaDataMetaId.eventEnergyScale);
                auto z_hist_shape = GetHistogram(anaDataMetaId_shape, eventRegion, originalZcategory.name, hist_name);
                if (z_hist_shape){
                    TH1D& z_hist = CloneHistogram(anaDataMetaId, eventRegion, newZcategory.name, *z_hist_shape);
                    RenormalizeHistogram(z_hist, yield, true);
                }
            }

        }
    }

    PhysicalValue CalculateYieldsForQCD(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                        EventRegion eventRegion, const std::string& hist_name, std::ostream& s_out)
    {
        const analysis::DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(analysis::DataCategoryType::QCD);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(analysis::DataCategoryType::Data);

        std::string bkg_yield_debug;
        const analysis::PhysicalValue bkg_yield =
                CalculateBackgroundIntegral(anaDataMetaId, eventRegion, hist_name, qcd.name, false, bkg_yield_debug);
        s_out << bkg_yield_debug;

        auto metaId_data = anaDataMetaId;
        metaId_data.eventEnergyScale = EventEnergyScale::Central;
        auto hist_data = GetHistogram(metaId_data, eventRegion, data.name, hist_name);
        if(!hist_data)
            throw exception("Unable to find data histograms for QCD yield estimation.");
        const auto data_yield = Integral(*hist_data, true);
        const PhysicalValue yield = data_yield - bkg_yield;
        s_out << "Data yield = " << data_yield << "\nData-MC yield = " << yield << std::endl;
        if(yield.GetValue() < 0) {
            if(&s_out != &std::cout)
                std::cout << bkg_yield_debug << "\nData yield = " << data_yield << std::endl;
            throw exception("Negative QCD yield for histogram '%1%' in %2% %3%.") % hist_name
                    % anaDataMetaId.eventCategory % eventRegion;
        }
        return yield;
    }

    const std::string& ChannelName() const { return __Channel_names<>::names.EnumToString(ChannelId()); }
    const std::string& ChannelNameLatex() const { return __Channel_names_latex.EnumToString(ChannelId()); }

    EventAnalyzerData& GetAnaData(const EventAnalyzerDataId& anaDataId)
    {
        return anaDataCollection.Get<FirstLeg>(anaDataId);
    }

    root_ext::SmartHistogram<TH1D>* GetHistogram(const EventAnalyzerDataId& anaDataId, const std::string& histogramName)
    {
        return GetAnaData(anaDataId).template GetPtr<TH1D>(histogramName);
    }

    root_ext::SmartHistogram<TH1D>* GetHistogram(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                                 EventRegion eventRegion, const std::string& dataCategoryName,
                                                 const std::string& histogramName)
    {
        return GetHistogram(anaDataMetaId.MakeId(eventRegion, dataCategoryName), histogramName);
    }

    root_ext::SmartHistogram<TH1D>* GetSignalHistogram(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                                       const std::string& dataCategoryName,
                                                       const std::string& histogramName)
    {
        return GetHistogram(anaDataMetaId, EventRegion::OS_Isolated, dataCategoryName, histogramName);
    }

    TH1D& CloneHistogram(const EventAnalyzerDataId& anaDataId, const root_ext::SmartHistogram<TH1D>& originalHistogram)
    {
        return GetAnaData(anaDataId).Clone(originalHistogram);
    }

    TH1D& CloneHistogram(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId, EventRegion eventRegion,
                         const std::string& dataCategoryName, const root_ext::SmartHistogram<TH1D>& originalHistogram)
    {
        return CloneHistogram(anaDataMetaId.MakeId(eventRegion, dataCategoryName), originalHistogram);
    }

    TH1D& CloneSignalHistogram(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                               const std::string& dataCategoryName,
                               const root_ext::SmartHistogram<TH1D>& originalHistogram)
    {
        return CloneHistogram(anaDataMetaId, EventRegion::OS_Isolated, dataCategoryName, originalHistogram);
    }

    static EventInfoBase::BjetPair SelectBjetPair(const ntuple::Event& event, bool order_bjet_by_csv)
    {
        if(!order_bjet_by_csv)
            return EventInfoBase::BjetPair(0, 1);

        std::vector<std::pair<size_t,double>> csvPosition;
        //std::vector<size_t> selected_jet_ids;
        for(size_t n = 0; n < event.jets_csv.size(); ++n) {
            const auto& p4 = event.jets_p4.at(n);
            if(p4.pt() < 30 || std::abs(p4.eta()) > 2.4) continue;
            std::pair<size_t,double> csvIndex(n, event.jets_csv.at(n));
            csvPosition.push_back(csvIndex);
        }

        auto pairCompare = [&](const std::pair<size_t,double> firstElem, const std::pair<size_t,double> secondElem) -> bool
        {  return firstElem.second > secondElem.second; } ;

        std::sort(csvPosition.begin(),csvPosition.end(),pairCompare);

        EventInfoBase::BjetPair selected_pair(event.jets_csv.size(), event.jets_csv.size() + 1);
        if(csvPosition.size() > 0)
            selected_pair.first = csvPosition.at(0).first;
        if(csvPosition.size() > 1)
            selected_pair.second = csvPosition.at(1).first;
        return selected_pair;
    }

    double ComputeWeight(const DataCategory& dataCategory, const ntuple::Event& event, double scale_factor)
    {
        static constexpr bool applybTagWeight = true;
        if(dataCategory.IsData()) return 1;
//        std::cout << event.p4_1.Pt() << " " << event.p4_2.Pt() << std::endl;
        return scale_factor * weights.GetTotalWeight(event, applybTagWeight, cuts::Htautau_2015::btag::CSVM);
    }

    void ProcessDataSource(const DataCategory& dataCategory, std::shared_ptr<ntuple::EventTuple> tree,
                           double scale_factor)
    {
        static constexpr bool order_bjet_by_csv = true;

//        const DataCategory& DYJets_incl = dataCategoryCollection.GetUniqueCategory(DataCategoryType::DYJets_incl);

        for(Long64_t current_entry = 0; current_entry < tree->GetEntries(); ++current_entry) {
            tree->GetEntry(current_entry);

            const EventInfoBase::BjetPair selected_bjet_pair = SelectBjetPair(tree->data(), order_bjet_by_csv);
            EventInfo event(tree->data(), selected_bjet_pair);

            // TODO
//            const int HTBin = 0;
//            if (dataCategory.name == DYJets_incl.name && HTBin != 0) continue;

            const EventCategoryVector eventCategories = DetermineEventCategories(event->jets_csv,
                                                                                 selected_bjet_pair,
                                                                                 0,
                                                                                 cuts::Htautau_2015::btag::CSVL,
                                                                                 cuts::Htautau_2015::btag::CSVM,
                                                                                 false);
            double weight = std::numeric_limits<double>::quiet_NaN();
            for(auto eventCategory : eventCategories) {
                if (!EventCategoriesToProcess().count(eventCategory)) continue;
                const EventRegion eventRegion = DetermineEventRegion(event, eventCategory);
                if(!EventRegionsToProcess().count(eventRegion)) continue;

                const EventSubCategorySet subCategories = DetermineEventSubCategories(event);
                for(auto subCategory : subCategories) {
                    if(!EventSubCategoriesToProcess().count(subCategory)) continue;
                    const EventAnalyzerDataId data_id(eventCategory, subCategory, eventRegion,
                                                     event.GetEnergyScale(), dataCategory.name);
                    if(std::isnan(weight))
                        weight = ComputeWeight(dataCategory, *event, scale_factor);
                    anaDataCollection.Fill(data_id, event, weight);
                }
            }
        }
    }

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

    std::string FullDataCardName(const std::string& datacard_name, EventEnergyScale eventEnergyScale) const
    {
        if(eventEnergyScale == EventEnergyScale::Central)
            return datacard_name;

        std::string channel_name = ChannelName();
        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);
        std::ostringstream full_name;
        full_name << datacard_name << "_CMS_scale_";
        if(eventEnergyScale == EventEnergyScale::TauUp || eventEnergyScale == EventEnergyScale::TauDown)
            full_name << "t_" << channel_name;
        else if(eventEnergyScale == EventEnergyScale::JetUp || eventEnergyScale == EventEnergyScale::JetDown)
            full_name << "j";
        else if(eventEnergyScale == EventEnergyScale::BtagEfficiencyUp
                || eventEnergyScale == EventEnergyScale::BtagEfficiencyDown)
            full_name << "btagEff";
        else if(eventEnergyScale == EventEnergyScale::BtagFakeUp || eventEnergyScale == EventEnergyScale::BtagFakeDown)
            full_name << "btagFake";
        else
            throw exception("Unsupported event energy scale %1%.") % eventEnergyScale;
        full_name << "_13TeV";
        if(eventEnergyScale == EventEnergyScale::TauUp || eventEnergyScale == EventEnergyScale::JetUp ||
                eventEnergyScale == EventEnergyScale::BtagEfficiencyUp ||
                eventEnergyScale == EventEnergyScale::BtagFakeUp )
            full_name << "Up";
        else
            full_name << "Down";
        return full_name.str();
    }

    void ProduceFileForLimitsCalculation(const std::string& hist_name, EventSubCategory eventSubCategory,
                                         typename EventAnalyzerData::HistogramAccessor histogramAccessor)
    {
        static const std::map<EventCategory, std::string> categoryToDirectoryNameSuffix = {
            { EventCategory::Inclusive, "inclusive" }, { EventCategory::TwoJets_ZeroBtag, "2jet0tag" },
            { EventCategory::TwoJets_OneBtag, "2jet1tag" }, { EventCategory::TwoJets_TwoBtag, "2jet2tag" }
        };

        static const std::map<std::string, std::string> channelNameForFolder = {
            { "eTau", "eleTau" }, { "muTau", "muTau" }, { "tauTau", "tauTau" }
        };

        static const double tiny_value = 1e-9;
        static const double tiny_value_error = tiny_value;

        std::ostringstream s_file_name;
        s_file_name << args.outputFileName() << "_" << hist_name;
        if(eventSubCategory != EventSubCategory::NoCuts)
            s_file_name << "_" << eventSubCategory;
        s_file_name << ".root";
        const std::string file_name = s_file_name.str();

        auto outputFile = root_ext::CreateRootFile(file_name);
        for(EventCategory eventCategory : EventCategoriesToProcess()) {
            if(!categoryToDirectoryNameSuffix.count(eventCategory)) continue;
            const std::string directoryName = channelNameForFolder.at(ChannelName()) + "_"
                    + categoryToDirectoryNameSuffix.at(eventCategory);
            outputFile->mkdir(directoryName.c_str());
            TDirectory* directory = outputFile->GetDirectory(directoryName.c_str());
            for(const DataCategory* dataCategory : dataCategoryCollection.GetCategories(DataCategoryType::Limits)) {
                if(!dataCategory->datacard.size())
                    throw exception("Empty datacard name for data category '%1%'.") % dataCategory->name;
                for(const EventEnergyScale& eventEnergyScale : EventEnergyScaleToProcess()) {
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
                    const std::string full_datacard_name = FullDataCardName(dataCategory->datacard, eventEnergyScale);
                    hist->Scale(dataCategory->limits_sf);
                    root_ext::WriteObject(*hist, directory, full_datacard_name);

                    if(eventEnergyScale == EventEnergyScale::Central && dataCategory->datacard == "ZL") {
                        std::string channel_name = ChannelName();
                        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);
                        const std::string name_syst_prefix = dataCategory->datacard + "_CMS_htt_"
                                + dataCategory->datacard + "Scale_" + channel_name + "_8TeV";
                        const std::string name_syst_up = name_syst_prefix + "Up";
                        const std::string name_syst_down = name_syst_prefix + "Down";
                        std::shared_ptr<TH1D> hist_syst_up(new TH1D(*hist));
                        hist_syst_up->Scale(1.02);
                        root_ext::WriteObject(*hist_syst_up, directory, name_syst_up);
                        std::shared_ptr<TH1D> hist_syst_down(new TH1D(*hist));
                        hist_syst_down->Scale(0.98);
                        root_ext::WriteObject(*hist_syst_down, directory, name_syst_down);
                    }
                    //added shape systematics for QCD
                    if(eventEnergyScale == EventEnergyScale::Central && dataCategory->datacard == "QCD_alternative") {
                        std::string channel_name = ChannelName();
                        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);
                        const std::string name_syst_prefix = "QCD_CMS_htt_" + dataCategory->datacard
                                + "Shape_" + channel_name + "_8TeV";
                        const std::string name_syst_up = name_syst_prefix + "Up";
                        const std::string name_syst_down = name_syst_prefix + "Down";
                        std::shared_ptr<TH1D> hist_syst_up(new TH1D(*hist));
                        root_ext::WriteObject(*hist_syst_up, directory, name_syst_up);
                        std::shared_ptr<TH1D> hist_syst_down(new TH1D(*hist));
                        root_ext::WriteObject(*hist_syst_down, directory, name_syst_down);
                    }
                }
            }
        }
    }

    void SubtractBackgroundHistograms(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                      EventRegion eventRegion, TH1D& histogram, const std::string& current_category,
                                      std::string& debug_info, std::string& negative_bins_info)
    {
        static const double correction_factor = 0.0000001;

        std::ostringstream ss_debug;

        ss_debug << "\nSubtracting background for '" << histogram.GetName() << "' in region " << eventRegion
                 << " for Event category '" << anaDataMetaId.eventCategory
                 << "' for data category '" << current_category
                 << "'.\nInitial integral: " << Integral(histogram, true) << ".\n";
        for (auto category : dataCategoryCollection.GetCategories(DataCategoryType::Background)) {
            if(category->IsComposit() || category->name == current_category || !category->isCategoryToSubtract)
                continue;

            ss_debug << "Sample '" << category->name << "': ";
            if(auto other_histogram = GetHistogram(anaDataMetaId, eventRegion, category->name, histogram.GetName())) {
                histogram.Add(other_histogram, -1);
                ss_debug << Integral(*other_histogram, true) << ".\n";
            } else
                ss_debug << "not found.\n";
        }

        const PhysicalValue original_Integral = Integral(histogram, true);
        ss_debug << "Integral after bkg subtraction: " << original_Integral << ".\n";
        debug_info = ss_debug.str();
        if (original_Integral.GetValue() < 0) {
            std::cout << debug_info << std::endl;
            throw exception("Integral after bkg subtraction is negative for histogram '%1%' in event category %2%"
                            " for event region %3%.") % histogram.GetName() % anaDataMetaId.eventCategory
                            % eventRegion;
        }

        std::ostringstream ss_negative;

        for (Int_t n = 1; n <= histogram.GetNbinsX(); ++n) {
            if (histogram.GetBinContent(n) >= 0) continue;
            const std::string prefix = histogram.GetBinContent(n) + histogram.GetBinError(n) >= 0 ? "WARNING" : "ERROR";

            ss_negative << prefix << ": " << histogram.GetName() << " Bin " << n << ", content = "
                        << histogram.GetBinContent(n) << ", error = " << histogram.GetBinError(n)
                        << ", bin limits=[" << histogram.GetBinLowEdge(n) << "," << histogram.GetBinLowEdge(n+1)
                        << "].\n";
            const double error = correction_factor - histogram.GetBinContent(n);
            const double new_error = std::sqrt(sqr(error) + sqr(histogram.GetBinError(n)));
            histogram.SetBinContent(n, correction_factor);
            histogram.SetBinError(n, new_error);
        }
        RenormalizeHistogram(histogram, original_Integral, true);
        negative_bins_info = ss_negative.str();
    }

    PhysicalValue CalculateBackgroundIntegral(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                              EventRegion eventRegion, const std::string& hist_name,
                                              const std::string& current_category,
                                              bool expect_at_least_one_contribution = false)
    {
        std::string debug_info;
        return CalculateBackgroundIntegral(anaDataMetaId, eventRegion, hist_name, current_category,
                                           expect_at_least_one_contribution, debug_info);
    }

    PhysicalValue CalculateBackgroundIntegral(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                              EventRegion eventRegion, const std::string& hist_name,
                                              const std::string& current_category,
                                              bool expect_at_least_one_contribution, std::string& debug_info)
    {
        DataCategoryPtrSet bkg_dataCategories;
        for (auto category : dataCategoryCollection.GetCategories(DataCategoryType::Background)) {
            if(category->IsComposit() || category->name == current_category || !category->isCategoryToSubtract )
                continue;
            bkg_dataCategories.insert(category);
        }

        return CalculateFullIntegral(anaDataMetaId, eventRegion, hist_name, bkg_dataCategories,
                                     expect_at_least_one_contribution, debug_info);
    }

    PhysicalValue CalculateFullIntegral(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                        EventRegion eventRegion, const std::string& hist_name,
                                        const DataCategoryPtrSet& dataCategories,
                                        bool expect_at_least_one_contribution = false)
    {
        std::string debug_info;
        return CalculateFullIntegral(anaDataMetaId, eventRegion, hist_name, dataCategories,
                                     expect_at_least_one_contribution, debug_info);
    }

    PhysicalValue CalculateFullIntegral(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                        EventRegion eventRegion, const std::string& hist_name,
                                        const DataCategoryPtrSet& dataCategories, bool expect_at_least_one_contribution,
                                        std::string& debug_info)
    {
        PhysicalValue integral;
        bool hist_found = false;

        std::ostringstream ss_debug;

        ss_debug << "\nCalculating full integral for '" << hist_name << "' in '"
                 << anaDataMetaId.MakeMetaId(eventRegion)
                 << "' considering data categories (" << dataCategories << ").\n";

        for(const auto& dataCategory : dataCategories) {
            ss_debug << "Sample '" << dataCategory->name << "': ";

            if(auto hist = GetHistogram(anaDataMetaId, eventRegion, dataCategory->name, hist_name)) {
                hist_found = true;
                const auto hist_integral = Integral(*hist, true);
                integral += hist_integral;
                ss_debug << hist_integral << ".\n";
            } else
                ss_debug << "not found.\n";
        }

        ss_debug << "Full integral: " << integral << ".\n";
        debug_info = ss_debug.str();

        if(!hist_found && expect_at_least_one_contribution) {
            std::cout << debug_info << std::endl;
            throw exception("No histogram with name '%1%' was found in the given data category set (%2%),"
                            " in eventCategory: '%3%' to calculate full integral.") % hist_name % dataCategories
                            % anaDataMetaId.MakeMetaId(eventRegion);
        }

        return integral;
    }

private:

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

        static const std::set<EventSubCategory> sidebandSubCategories = {
            EventSubCategory::OutsideMassWindow, EventSubCategory::KinematicFitConvergedOutsideMassWindow
        };

        const auto findRegionId = [&]() -> size_t {
            if(sidebandSubCategories.count(subCategory) || !histogramsToBlind.count(hist_name))
                return 0;
            return histogramsToBlind.at(hist_name);
        };

        const size_t regionId = findRegionId();
        if(regionId >= blindingRegions.size())
            throw exception("Bad blinding region index = %1%.") % regionId;
        return blindingRegions.at(regionId);
    }

    void PrintTables(const std::string& name_suffix, const std::wstring& sep)
    {
        std::wofstream of(args.outputFileName() + "_" + name_suffix + ".csv");

        static const std::set< std::pair<std::string, EventSubCategory> > interesting_histograms = {
            { EventAnalyzerData::m_sv_Name(), EventSubCategory::NoCuts },
            { EventAnalyzerData::m_sv_Name(), EventSubCategory::MassWindow },
            { EventAnalyzerData::m_ttbb_kinfit_Name(), EventSubCategory::KinematicFitConverged },
            { EventAnalyzerData::m_ttbb_kinfit_Name(), EventSubCategory::KinematicFitConvergedWithMassWindow }
        };

        for(const auto& hist_entry : interesting_histograms)
            PrintTables(of, sep, hist_entry.first, hist_entry.second, false, true);

        of.flush();
        of.close();
    }

    void PrintTables(std::wostream& of, const std::wstring& sep, const std::string& hist_name,
                     EventSubCategory subCategory, bool includeOverflow, bool includeError)
    {
        of << std::wstring(hist_name.begin(), hist_name.end());

        std::wstring table_name_suffix = L"";
        if(includeOverflow && includeError)
            table_name_suffix = L" with overflow and error";
        else if(includeOverflow && !includeError)
            table_name_suffix = L" with overflow";
        else if(!includeOverflow && includeError)
            table_name_suffix = L" with error";
        of << table_name_suffix << sep;

        for (EventCategory eventCategory : EventCategoriesToProcess())
            of << eventCategory << sep;
        of << std::endl;

        for (const DataCategory* dataCategory : dataCategoryCollection.GetAllCategories()) {
            of << std::wstring(dataCategory->title.begin(), dataCategory->title.end()) << sep;
            for (EventCategory eventCategory : EventCategoriesToProcess()) {
                const EventAnalyzerDataMetaId_noRegion_noName meta_id(eventCategory, subCategory,
                                                                     EventEnergyScale::Central);

                if( TH1D* histogram = GetSignalHistogram(meta_id, dataCategory->name, hist_name) ) {
                    const PhysicalValue integral = Integral(*histogram, includeOverflow);
                    of << integral.ToString<wchar_t>(includeError, false) << sep;
                }
                else
                    of << "not found" << sep;
            }
            of << std::endl;
        }
        of << std::endl << std::endl;
    }

    void ProcessCompositDataCategories(const EventAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                       const std::string& hist_name)
    {
        for (analysis::EventRegion eventRegion : analysis::AllEventRegions) {
            for(const DataCategory* composit : dataCategoryCollection.GetCategories(DataCategoryType::Composit)) {
                for(const std::string& sub_name : composit->sub_categories) {
                    const DataCategory& sub_category = dataCategoryCollection.FindCategory(sub_name);
                    auto sub_hist = GetHistogram(anaDataMetaId, eventRegion, sub_category.name, hist_name);
                    if(!sub_hist) continue;
                    if(auto composit_hist = GetHistogram(anaDataMetaId, eventRegion, composit->name, hist_name))
                        composit_hist->Add(sub_hist);
                    else
                        CloneHistogram(anaDataMetaId, eventRegion, composit->name, *sub_hist);
                }
            }
        }
    }

protected:
    AnalyzerArguments args;
    DataCategoryCollection dataCategoryCollection;
    EventAnalyzerDataCollection anaDataCollection;
    mc_corrections::EventWeights weights;
};

} // namespace analysis
