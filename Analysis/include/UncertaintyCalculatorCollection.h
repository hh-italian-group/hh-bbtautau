/*! Collection of methods to calculate uncertainties.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "UncertaintyConfiguration.h"
#include "h-tautau/Analysis/include/FlatAnalyzerDataCollection.h"

namespace analysis {
namespace limits {

namespace uncertainty_names {

const std::string TauES = "scale_t";
const std::string JetES = "scale_j";
const std::string Btag_efficiency = "eff_b";
const std::string Btag_fake = "fake_b";
const std::string TTbar_normalization = "ttbarNorm";
const std::string ZTT_extrapolation = "extrap_ztt";
const std::string ZLL_FakeTau = "ZLL_FakeTau";
const std::string JetFakeTau = "JetFakeTau";
const std::string LeptonFakeTau = "LeptonFakeTau";
const std::string QCD = "QCDSyst";

} // namespace uncertainty_names

class UncertaintyCalculatorCollection {
private:
    typedef std::function< UncertaintyInterval (EventCategory, const std::string&) > UncertaintyCalculator;
    typedef std::map<std::string, UncertaintyCalculator> UncertaintyCalculatorMap;

    template<typename MethodPtr>
    void Bind(const std::string& name, MethodPtr method_ptr)
    {
        using namespace std::placeholders;
        calculator_map[name] = std::bind(method_ptr, this, _1, _2);
    }

    static double DefaultPrecision() { return 0.001; }

public:
    UncertaintyCalculatorCollection(const UncertaintyDescriptorCollection& _uncertainties,
                                    const DataCategoryCollection& _dataCategories,
                                    const FlatAnalyzerDataCollectionReader& _reader,
                                    EventSubCategory _eventSubCategory, const std::string& _referenceHistName)
        : uncertainties(&_uncertainties), dataCategories(&_dataCategories), reader(&_reader),
          eventSubCategory(_eventSubCategory), referenceHistName(_referenceHistName)
    {
        using namespace uncertainty_names;

        Bind(Btag_efficiency, &UncertaintyCalculatorCollection::CalculateBtagEfficiencyUnc);
        Bind(Btag_fake, &UncertaintyCalculatorCollection::CalculateBtagFakeUnc);
        Bind(TTbar_normalization, &UncertaintyCalculatorCollection::CalculateTTUnc);
        Bind(ZTT_extrapolation, &UncertaintyCalculatorCollection::CalculateZTTextrapUnc);
        Bind(ZLL_FakeTau, &UncertaintyCalculatorCollection::CalculateZFakeTauUnc);
        Bind(QCD, &UncertaintyCalculatorCollection::CalculateQCDUnc);
    }

    UncertaintyInterval Calculate(const std::string& unc_name, EventCategory event_category,
                                  const std::string& sample_name)
    {
        if(!calculator_map.count(unc_name))
            throw exception("Calculator for uncertainty '") << unc_name << "' not found.";
        return calculator_map.at(unc_name)(event_category, sample_name);
    }

    static PhysicalValue CombineUpDownUncertainties(const UncertaintyInterval& unc)
    {
        return ( unc.up + (PhysicalValue::Two - unc.down) ) / PhysicalValue::Two;
    }

private:
    PhysicalValue GetYield(EventCategory eventCategory, const std::string& dataCategoryName,
                           EventRegion eventRegion = EventRegion::OS_Isolated,
                           EventEnergyScale eventEnergyScale = EventEnergyScale::Central, bool expect_hist = true) const
    {
        const FlatAnalyzerDataId id(eventCategory, eventSubCategory, eventRegion, eventEnergyScale, dataCategoryName);
        auto hist = reader->GetHistogram<TH1D>(id, referenceHistName);
        if(hist)
            return Integral(*hist, true);
        if(expect_hist)
            std::cout << "Histogram '" << referenceHistName << "' in " << id << " not found. Considering zero yield.\n";
        return PhysicalValue::Zero;
    }

    const std::string& GetDataCategoryName(const std::string& datacard) const
    {
        return dataCategories->FindCategoryForDatacard(datacard).name;
    }

    void AddUncertainty(PhysicalValue& physicalValue, const std::string& unc_name,
                        const std::string& sample_name = "") const
    {
        const auto& unc_desc = uncertainties->Get(unc_name);
        const double value = sample_name.size() && unc_desc.sample_values.count(sample_name)
                ? unc_desc.sample_values.at(sample_name) : unc_desc.value;
        physicalValue.AddSystematicUncertainty(unc_name, value - 1);
    }

    UncertaintyInterval CalculateEnergyScaleRelatedUncertainty(EventCategory eventCategory,
                                                               const std::string& dataCategoryName,
                                                               EventRegion eventRegion,
                                                               EventEnergyScale up_energy_scale,
                                                               EventEnergyScale down_energy_scale,
                                                               bool expect_hist = true) const
    {
        const PhysicalValue n_central = GetYield(eventCategory, dataCategoryName, eventRegion,
                                                 EventEnergyScale::Central, expect_hist);
        const PhysicalValue n_up = GetYield(eventCategory, dataCategoryName, eventRegion, up_energy_scale, expect_hist);
        const PhysicalValue n_down = GetYield(eventCategory, dataCategoryName, eventRegion, down_energy_scale,
                                              expect_hist);

        return UncertaintyInterval(n_down / n_central, n_up / n_central);
    }


    UncertaintyInterval CalculateBtagEfficiencyUnc(EventCategory eventCategory,
                                                   const std::string& sample_name) const
    {
        const auto& dataCategoryName = GetDataCategoryName(sample_name);
        return CalculateEnergyScaleRelatedUncertainty(eventCategory, dataCategoryName, EventRegion::OS_Isolated,
                                                      EventEnergyScale::BtagEfficiencyUp,
                                                      EventEnergyScale::BtagEfficiencyDown);
    }

    UncertaintyInterval CalculateBtagFakeUnc(EventCategory eventCategory,
                                             const std::string& sample_name) const
    {
        const auto& dataCategoryName = GetDataCategoryName(sample_name);
        return CalculateEnergyScaleRelatedUncertainty(eventCategory, dataCategoryName, EventRegion::OS_Isolated,
                                                      EventEnergyScale::BtagFakeUp,
                                                      EventEnergyScale::BtagFakeDown);
    }

    const PhysicalValue& GetZTT_SF(EventCategory eventCategory)
    {
        if(!ztt_sf_map.count(eventCategory)) {
            const std::string& DY_emb_name = dataCategories->GetUniqueCategory(DataCategoryType::Embedded).name;
            const std::string& TT_emb_name = dataCategories->GetUniqueCategory(DataCategoryType::TT_Embedded).name;

            PhysicalValue DY_cat = GetYield(eventCategory, DY_emb_name);
            PhysicalValue TT_cat = GetYield(eventCategory, TT_emb_name);
            PhysicalValue DY_incl = GetYield(EventCategory::Inclusive, DY_emb_name);
            PhysicalValue TT_incl = GetYield(EventCategory::Inclusive, TT_emb_name);

            AddUncertainty(TT_cat, uncertainty_names::TTbar_normalization);
            AddUncertainty(TT_incl, uncertainty_names::TTbar_normalization);

            const PhysicalValue sf = (DY_cat - TT_cat) / (DY_incl - TT_incl);
            std::cout << "ZTT SF = " << sf << std::endl;
            ztt_sf_map[eventCategory] = sf;
        }
        return ztt_sf_map.at(eventCategory);
    }

    UncertaintyInterval CalculateTTUnc(EventCategory eventCategory, const std::string& sample_name)
    {
        if(sample_name != "ZTT")
            throw exception("Sample '") << sample_name << "' not supported by "
                                        << uncertainty_names::TTbar_normalization << " uncertainty calculator.";

        const PhysicalValue& sf = GetZTT_SF(eventCategory);
        const double unc = sf.GetRelativeSystematicUncertainty(uncertainty_names::TTbar_normalization);
        return UncertaintyInterval(PhysicalValue(unc, DefaultPrecision()));
    }

    UncertaintyInterval CalculateZTTextrapUnc(EventCategory eventCategory, const std::string& sample_name)
    {
        if(sample_name != "ZTT")
            throw exception("Sample '") << sample_name << "' not supported by "
                                        << uncertainty_names::ZTT_extrapolation << " uncertainty calculator.";

        const PhysicalValue& sf = GetZTT_SF(eventCategory);
        const double unc = sf.GetRelativeStatisticalError();
        return UncertaintyInterval(PhysicalValue(unc, DefaultPrecision()));
    }

    UncertaintyInterval CalculateZFakeTauUnc(EventCategory eventCategory, const std::string& sample_name) const
    {
        if(sample_name != "ZLL")
            throw exception("Sample '") << sample_name << "' not supported by "
                                        << uncertainty_names::ZLL_FakeTau << " uncertainty calculator.";

        const std::string ZJ_name = dataCategories->GetUniqueCategory(DataCategoryType::ZJ_MC).name;
        const std::string ZL_name = dataCategories->GetUniqueCategory(DataCategoryType::ZL_MC).name;

        PhysicalValue n_ZJ = GetYield(eventCategory, ZJ_name);
        PhysicalValue n_ZL = GetYield(eventCategory, ZL_name);

        AddUncertainty(n_ZJ, uncertainty_names::JetFakeTau);
        AddUncertainty(n_ZL, uncertainty_names::LeptonFakeTau);

        const PhysicalValue n_ZLL = n_ZJ + n_ZL;
        std::cout << "ZLL yield = " << n_ZLL << std::endl;

        const double unc = n_ZLL.GetRelativeFullError();
        return UncertaintyInterval(PhysicalValue(unc, DefaultPrecision()));
    }

    const std::set<std::string>& GetBackgroundNames()
    {
        if(!bkg_names.size()) {
            const std::set<std::string> qcd_names = {
                dataCategories->GetUniqueCategory(DataCategoryType::QCD).name,
                dataCategories->GetUniqueCategory(DataCategoryType::QCD_alternative).name
            };

            for (auto category : dataCategories->GetCategories(DataCategoryType::Background)) {
                if(category->IsComposit() || qcd_names.count(category->name) || !category->isCategoryToSubtract )
                    continue;
                bkg_names.insert(category->name);
            }
        }
        return bkg_names;
    }

    void AddAllUncertainties(PhysicalValue& value, EventCategory eventCategory, const std::string& dataCategoryName,
                             EventRegion eventRegion) const
    {
        const auto apply_es_unc = [&] (const std::string& name, EventEnergyScale up, EventEnergyScale down) -> void {
            const auto es_interval = CalculateEnergyScaleRelatedUncertainty(eventCategory, dataCategoryName,
                                                                            eventRegion, up, down, false);
            const PhysicalValue es_unc = CombineUpDownUncertainties(es_interval);
            if(!std::isnan(es_unc.GetValue()))
                value.AddSystematicUncertainty(name, es_unc.GetValue() - 1);
        };

        apply_es_unc(uncertainty_names::TauES, EventEnergyScale::TauUp, EventEnergyScale::TauDown);
        apply_es_unc(uncertainty_names::JetES, EventEnergyScale::JetUp, EventEnergyScale::JetDown);
        apply_es_unc(uncertainty_names::Btag_efficiency, EventEnergyScale::BtagEfficiencyUp,
                     EventEnergyScale::BtagEfficiencyDown);
        apply_es_unc(uncertainty_names::Btag_fake, EventEnergyScale::BtagFakeUp, EventEnergyScale::BtagFakeDown);

        const DataCategory& dataCategory = dataCategories->FindCategory(dataCategoryName);
        for(const std::string& unc_name : dataCategory.uncertainties)
            AddUncertainty(value, unc_name, dataCategory.datacard);
    }

    PhysicalValue CalculateBackgroundYield(EventCategory eventCategory, EventRegion eventRegion)
    {
        PhysicalValue total_yield;

        for(const std::string& bkg_name : GetBackgroundNames()) {
            PhysicalValue bkg_yield = GetYield(eventCategory, bkg_name, eventRegion, EventEnergyScale::Central, false);
            AddAllUncertainties(bkg_yield, eventCategory, bkg_name, eventRegion);
            std::cout << "  " << bkg_name << ": " << bkg_yield << ".\n";
            total_yield += bkg_yield;
        }
        std::cout << "Total bkg yield: " << total_yield << ".\n";
        return total_yield;
    }

    PhysicalValue CalculateQcdYield(EventCategory eventCategory, EventRegion eventRegion, bool consider_data_stat)
    {

        std::cout << eventCategory << "/" << eventRegion << "\n";
        const std::string& data_name = dataCategories->GetUniqueCategory(DataCategoryType::Data).name;
        const PhysicalValue data_yield_stat = GetYield(eventCategory, data_name, eventRegion);
        const PhysicalValue data_yield = consider_data_stat
                ? data_yield_stat : PhysicalValue(data_yield_stat.GetValue());
        const PhysicalValue bkg_yield = CalculateBackgroundYield(eventCategory, eventRegion);
        const PhysicalValue qcd_yield = data_yield - bkg_yield;

        std::cout << "Data yield: " << data_yield_stat << ".\n"
                  << "QCD yield: " << qcd_yield << ".\n";

        return qcd_yield;
    }

    UncertaintyInterval CalculateQCDUnc(EventCategory eventCategory, const std::string& sample_name)
    {
        if(sample_name != "QCD")
            throw exception("Sample '") << sample_name << "' not supported by "
                                        << uncertainty_names::QCD << " uncertainty calculator.";

        static const analysis::EventCategorySet categories_to_loose = {
            EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag
        };

        EventCategory refEventCategory = eventCategory;
        if(categories_to_loose.count(eventCategory))
            refEventCategory = MediumToLoose_EventCategoryMap.at(eventCategory);


        const PhysicalValue yield_OSAntiIso = CalculateQcdYield(eventCategory, EventRegion::OS_AntiIsolated, true);
        const PhysicalValue yield_SSIso_ref = CalculateQcdYield(refEventCategory, EventRegion::SS_Isolated, true);
        const PhysicalValue yield_SSAntiIso = CalculateQcdYield(refEventCategory, EventRegion::SS_AntiIsolated, true);

        const PhysicalValue iso_antiIso_sf = yield_SSIso_ref / yield_SSAntiIso;
        const PhysicalValue qcd_yield = yield_OSAntiIso * iso_antiIso_sf;
        const double qcd_unc = qcd_yield.GetRelativeFullError();

        std::cout << "QCD iso/anti_iso SF: " << iso_antiIso_sf << ".\n"
                  << "QCD signal yiled: " << qcd_yield << ".\n"
                  << "QCD full uncertainty: " << qcd_unc << ".\n";

        const PhysicalValue yield_OSAntiIso_ref = CalculateQcdYield(refEventCategory, EventRegion::OS_AntiIsolated,
                                                                    true);
        const PhysicalValue yield_SSIso = CalculateQcdYield(eventCategory, EventRegion::SS_Isolated, true);

        const PhysicalValue os_ss_sf = yield_OSAntiIso_ref / yield_SSAntiIso;
        const PhysicalValue alt_qcd_yield = yield_SSIso * os_ss_sf;
        const double alt_qcd_unc = alt_qcd_yield.GetRelativeFullError();
        const PhysicalValue delta_qcd =
                (yield_OSAntiIso * yield_SSIso_ref - yield_OSAntiIso_ref * yield_SSIso) / yield_SSAntiIso;

        const PhysicalValue average_qcd = PhysicalValue::WeightedAverage({qcd_yield, alt_qcd_yield});
        const double qcd_cov = qcd_yield.Covariance(alt_qcd_yield);

        std::cout << "Alt QCD os/ss SF: " << os_ss_sf << ".\n"
                  << "Alt QCD signal yiled: " << alt_qcd_yield << ".\n"
                  << "Alt QCD full uncertainty: " << alt_qcd_unc << ".\n"
                  << "Delta QCD methods: " << delta_qcd << ".\n"
                  << "Average QCD signal value: " << average_qcd << ".\n"
                  << "QCD covariance: " << qcd_cov << ".\n";

        return UncertaintyInterval(PhysicalValue(qcd_unc, DefaultPrecision()));
    }

private:
    const UncertaintyDescriptorCollection* uncertainties;
    const DataCategoryCollection* dataCategories;
    const FlatAnalyzerDataCollectionReader* reader;
    EventSubCategory eventSubCategory;
    std::string referenceHistName;
    UncertaintyCalculatorMap calculator_map;
    std::map<EventCategory, PhysicalValue> ztt_sf_map;
    std::set<std::string> bkg_names;
};

} // namespace limits
} // namespace analysis
