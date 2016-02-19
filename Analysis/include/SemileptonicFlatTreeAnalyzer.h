/*! Definition of the base class for semi-leptonic flat tree analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "BaseFlatTreeAnalyzer.h"

namespace analysis {

class SemileptonicFlatTreeAnalyzer : public BaseFlatTreeAnalyzer {
public:
    SemileptonicFlatTreeAnalyzer(const DataCategoryCollection& _dataCategoryCollection, const std::string& _inputPath,
                                 const std::string& _outputFileName, bool applyPostFitCorrections, bool saveFullOutput)
         : BaseFlatTreeAnalyzer(_dataCategoryCollection, _inputPath, _outputFileName, applyPostFitCorrections,
                                saveFullOutput)
    {
    }

protected:
    static bool IsHighMtRegion(const ntuple::Flat& event, analysis::EventCategory eventCategory)
    {
        using namespace cuts;
        if (eventCategory == analysis::EventCategory::TwoJets_TwoBtag)
            return event.mt_1 > WjetsBackgroundEstimation::HighMtRegion_low &&
                    event.mt_1 < WjetsBackgroundEstimation::HighMtRegion_high;
        else
            return event.mt_1 > WjetsBackgroundEstimation::HighMtRegion;
    }

    static bool IsAntiIsolatedRegion(const ntuple::Flat& event)
    {
        using namespace cuts;        
        return event.pfRelIso_1 > IsolationRegionForLeptonicChannel::isolation_low &&
                event.pfRelIso_1 < IsolationRegionForLeptonicChannel::isolation_high;
    }

    virtual void CreateHistogramForZTT(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                       const std::string& hist_name, const PhysicalValueMap& ztt_yield_map,
                                       bool useEmbedded) override
    {
        const DataCategory& embedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Embedded);
        const DataCategory& ZTT_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_MC);
        const DataCategory& ZTT = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT);
        const DataCategory& ZTT_L = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_L);
        const DataCategory& TTembedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::TT_Embedded);

        for(const auto& eventRegionKey : ztt_yield_map) {
            const EventRegion eventRegion = eventRegionKey.first;
            const PhysicalValue& ztt_yield = eventRegionKey.second;
            auto ztt_l_hist = GetHistogram(anaDataMetaId, eventRegion, ZTT_L.name, hist_name);
            const std::string embeddedName = useEmbedded && eventRegion == EventRegion::OS_Isolated ?
                                             embedded.name : ZTT_MC.name;
            auto embedded_hist = GetHistogram(anaDataMetaId, eventRegion, embeddedName, hist_name);
            auto TTembedded_hist = GetHistogram(anaDataMetaId, eventRegion, TTembedded.name, hist_name);

            if (embedded_hist){
                TH1D& ztt_hist = CloneHistogram(anaDataMetaId, eventRegion, ZTT.name, *embedded_hist);
                if (TTembedded_hist && useEmbedded)
                    ztt_hist.Add(TTembedded_hist, -1);
                RenormalizeHistogram(ztt_hist, ztt_yield, true);
                if (ztt_l_hist)
                    ztt_hist.Add(ztt_l_hist);
            }
            if (!embedded_hist && ztt_l_hist)
                CloneHistogram(anaDataMetaId, eventRegion, ZTT.name, *ztt_l_hist);
        }
    }

    virtual analysis::PhysicalValue CalculateQCDYield(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                                      const std::string& hist_name,
                                                      DataCategoryType /*dataCategoryType*/,
                                                      std::ostream& s_out) override
    {
        static const PhysicalValue sf(1.06, 0.001);
        static const EventCategorySet categories = { EventCategory::TwoJets_TwoBtag };

        analysis::EventCategory refEventCategory = anaDataMetaId.eventCategory;
        if(categories.count(anaDataMetaId.eventCategory))
            refEventCategory = EventCategory::TwoJets_Inclusive;

        const FlatAnalyzerDataMetaId_noRegion_noName anaDataMetaId_ref(refEventCategory,
                                                                       anaDataMetaId.eventSubCategory,
                                                                       anaDataMetaId.eventEnergyScale);

        const PhysicalValue yield_SSIso =
                CalculateYieldsForQCD(anaDataMetaId_ref, EventRegion::SS_Isolated, hist_name, s_out);

        s_out << "yield_ssIso: " << yield_SSIso << "\n";
        if(refEventCategory == anaDataMetaId.eventCategory)
            return sf * yield_SSIso;

        const DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        auto hist_data_EvtCategory = GetHistogram(anaDataMetaId, EventRegion::SS_AntiIsolated, data.name, hist_name);
        if(!hist_data_EvtCategory)
            throw exception("Unable to find hist_data_EvtCategory for QCD scale factors estimation - SS AntiIso");
        const PhysicalValue yield_Data_EvtCategory = Integral(*hist_data_EvtCategory, true);

        auto hist_data_RefCategory =
                GetHistogram(anaDataMetaId_ref, EventRegion::SS_AntiIsolated, data.name, hist_name);
        if(!hist_data_RefCategory)
            throw exception("Unable to find hist_data_RefCategory for QCD scale factors estimation - SS AntiIso");
        const PhysicalValue yield_Data_RefCategory = Integral(*hist_data_RefCategory, true);

        const auto evt_ToRef_category_sf = yield_Data_EvtCategory / yield_Data_RefCategory;
        s_out << "evt_ToRef_category_sf: " << evt_ToRef_category_sf << "\n";

        return sf * yield_SSIso * evt_ToRef_category_sf;
    }

    virtual void EstimateQCD(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId, const std::string& hist_name,
                             const analysis::PhysicalValue& scale_factor, DataCategoryType dataCategoryType) override
    {
        static const EventCategorySet categories=
            { EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_AtLeastOneBtag };
        static const analysis::EventCategorySet inclusive_categories =
            { EventCategory::Inclusive, EventCategory::TwoJets_Inclusive};

        EventCategory refEventCategory = anaDataMetaId.eventCategory;
        if(categories.count(anaDataMetaId.eventCategory))
            refEventCategory = MediumToLoose_EventCategoryMap.at(anaDataMetaId.eventCategory);

        EventRegion eventRegion = EventRegion::SS_AntiIsolated;
        bool subtractOtherBkg = false;
        if(inclusive_categories.count(anaDataMetaId.eventCategory)) {
            eventRegion = EventRegion::SS_Isolated;
            subtractOtherBkg = true;
        }

        return EstimateQCDEx(anaDataMetaId, refEventCategory, eventRegion, hist_name, scale_factor, subtractOtherBkg,
                             dataCategoryType);
    }

    void EstimateQCDEx(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId, EventCategory refEventCategory,
                       EventRegion eventRegion, const std::string& hist_name, const PhysicalValue& scale_factor,
                       bool subtractOtherBkg, DataCategoryType dataCategoryType)
    {
        const DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(dataCategoryType);
        const DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        const FlatAnalyzerDataMetaId_noRegion_noName anaDataMetaId_ref(refEventCategory,
                                                                       anaDataMetaId.eventSubCategory,
                                                                       anaDataMetaId.eventEnergyScale);

        auto hist_shape_data = GetHistogram(anaDataMetaId_ref, eventRegion, data.name, hist_name);
        if(!hist_shape_data) return;

        TH1D& histogram = CloneHistogram(anaDataMetaId, EventRegion::OS_Isolated, qcd.name, *hist_shape_data);
        if (subtractOtherBkg){
            std::string debug_info, negative_bins_info;
            SubtractBackgroundHistograms(anaDataMetaId_ref, eventRegion, histogram, qcd.name, debug_info,
                                         negative_bins_info);
        }
        RenormalizeHistogram(histogram, scale_factor, true);
    }

    virtual PhysicalValueMap CalculateWjetsYields(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                                  const std::string& hist_name,
                                                  bool /*fullEstimate*/) override
    {
        PhysicalValueMap valueMap;

        const DataCategory& wjets = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets);
        const DataCategoryPtrSet& wjets_mc_categories =
                dataCategoryCollection.GetCategories(DataCategoryType::WJets_MC);
        const DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        static const EventCategorySet categoriesToRelax = { EventCategory::TwoJets_TwoBtag };
        EventCategory refEventCategory = anaDataMetaId.eventCategory;
        if(categoriesToRelax.count(anaDataMetaId.eventCategory))
            refEventCategory = MediumToLoose_EventCategoryMap.at(anaDataMetaId.eventCategory);

        const FlatAnalyzerDataMetaId_noRegion_noName anaDataMetaId_ref(refEventCategory,
                                                                       anaDataMetaId.eventSubCategory,
                                                                       anaDataMetaId.eventEnergyScale);

        static const EventRegionMap HighMt_LowMt_RegionMap_forW =
                { { EventRegion::OS_Iso_HighMt, EventRegion::OS_Isolated },
                  { EventRegion::SS_Iso_HighMt, EventRegion::SS_Isolated }};

        for (const auto& eventRegion : HighMt_LowMt_RegionMap_forW) {
            std::string bkg_yield_debug;
            const PhysicalValue bkg_yield = CalculateBackgroundIntegral(anaDataMetaId, eventRegion.first, hist_name,
                                                                        wjets.name, false, bkg_yield_debug);

            auto hist_data = GetHistogram(anaDataMetaId, eventRegion.first, data.name, hist_name);
            if(!hist_data)
                throw exception("Unable to find data histograms for Wjet scale factors estimation");
            std::cout << "Data Integral in Wjets Yield: " << Integral(*hist_data, true) << std::endl;

            const auto data_yield = Integral(*hist_data, true);
            const PhysicalValue yield = data_yield - bkg_yield;
            if(yield.GetValue() < 0) {
                std::cout << bkg_yield_debug << "\nData yield = " << data_yield << std::endl;
                throw exception("Negative Wjets yield for histogram '")
                        << hist_name << "' in " << anaDataMetaId.eventCategory << " " << eventRegion.first << ".";
            }

            const PhysicalValue n_HighMt_mc =
                    CalculateFullIntegral(anaDataMetaId_ref, eventRegion.first, hist_name, wjets_mc_categories, true);
            const PhysicalValue n_LowMt_mc =
                    CalculateFullIntegral(anaDataMetaId_ref, eventRegion.second, hist_name, wjets_mc_categories, true);
            const PhysicalValue ratio_LowToHighMt = n_LowMt_mc / n_HighMt_mc;
            valueMap[eventRegion.second] = yield * ratio_LowToHighMt;
        }

        return valueMap;
    }

//    virtual PhysicalValueMap CalculateWjetsYields(const analysis::FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
//                                                  const std::string& hist_name, bool fullEstimate) override
//    {
//        using analysis::EventCategory;

//        const analysis::DataCategoryPtrSet& wjets_mc_categories =
//                dataCategoryCollection.GetCategories(analysis::DataCategoryType::WJets_MC);

//        PhysicalValueMap valueMap;
//        for (analysis::EventRegion eventRegion : analysis::QcdRegions)
//            valueMap[eventRegion] = CalculateFullIntegral(anaDataMetaId, eventRegion, hist_name,
//                                                          wjets_mc_categories);
//        return valueMap;
//    }

    virtual void CreateHistogramForVVcategory(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                              const std::string& hist_name) override
    {
        const std::map<DataCategoryType, DataCategoryType> diboson_category_map = {
            { DataCategoryType::DiBoson_MC, DataCategoryType::DiBoson }
        };

        for (const auto& diboson_category : diboson_category_map){
            const DataCategory& originalVVcategory = dataCategoryCollection.GetUniqueCategory(diboson_category.first);
            const DataCategory& newVVcategory = dataCategoryCollection.GetUniqueCategory(diboson_category.second);

            for(EventRegion eventRegion : AllEventRegions) {
                auto vv_hist_shape = GetHistogram(anaDataMetaId, eventRegion, originalVVcategory.name, hist_name);
                if (vv_hist_shape)
                    CloneHistogram(anaDataMetaId, eventRegion, newVVcategory.name, *vv_hist_shape);
            }
        }
    }
};

} // namespace analysis
