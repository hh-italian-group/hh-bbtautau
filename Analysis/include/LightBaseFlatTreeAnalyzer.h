/*! Definition of LightBaseFlatTreeAnalyzer class, the base class for separate studies on flat trees.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <iostream>
#include <cmath>
#include <set>
#include <list>
#include <locale>
#include <numeric>
#include <algorithm>

#include <TColor.h>
#include <TLorentzVector.h>

#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "h-tautau/Analysis/include/FlatEventInfo.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTools/Core/include/exception.h"
#include "h-tautau/Analysis/include/Particles.h"
#include "AnalysisTools/Print/include/RootPrintToPdf.h"

#include "MVASelections/include/MvaReader.h"

#include "h-tautau/Analysis/include/FlatAnalyzerDataCollection.h"

namespace analysis {

class LightBaseFlatTreeAnalyzer {
public:
    typedef std::map<std::string, FlatEventInfo::BjetPair> PairSelectionMap;
    typedef FlatAnalyzerDataMetaId_noName MetaId;
    typedef std::set<MetaId> MetaIdSet;
    typedef std::map<FlatEventInfo::BjetPair, FlatEventInfoPtr> FlatEventInfoMap;

    LightBaseFlatTreeAnalyzer(const std::string& inputFileName, const std::string& outputFileName)
        : inputFile(root_ext::OpenRootFile(inputFileName)),
          outputFile(root_ext::CreateRootFile(outputFileName)),
          flatTree(new ntuple::FlatTree("flatTree", inputFile.get(), true)), recalc_kinfit(false), do_retag(true)
    {
        TH1::SetDefaultSumw2();
    }

    virtual ~LightBaseFlatTreeAnalyzer() {}

    void Run()
    {
        for(Long64_t current_entry = 0; current_entry < flatTree->GetEntries(); ++current_entry) {
            eventInfoMap.clear();
            flatTree->GetEntry(current_entry);
            const ntuple::Flat& event = flatTree->data;
            const auto& pairSelectionMap = SelectBjetPairs(event);
            for(const auto& selection_entry : pairSelectionMap) {
                const std::string& selection_label = selection_entry.first;
                const FlatEventInfo::BjetPair& bjet_pair = selection_entry.second;
                const MetaIdSet metaIds = CreateMetaIdSet(event, bjet_pair);
                for (const auto& metaId : metaIds) {
                    const FlatEventInfo& eventInfo = GetFlatEventInfo(event, bjet_pair);
                    AnalyzeEvent(eventInfo, metaId, selection_label);
                }
            }
        }
        EndOfRun();
    }

protected:
    virtual void AnalyzeEvent(const FlatEventInfo& eventInfo, const MetaId& metaId,
                              const std::string& selectionLabel) = 0;

    virtual void EndOfRun() {}

    virtual PairSelectionMap SelectBjetPairs(const ntuple::Flat& /*event*/)
    {
        PairSelectionMap pairMap;
        pairMap["CSV"] = FlatEventInfo::BjetPair(0, 1);
        return pairMap;
    }

    virtual const EventEnergyScaleSet GetEnergyScalesToProcess() const { return AllEventEnergyScales; }
    virtual const EventCategorySet& GetCategoriesToProcess() const { return AllEventCategories; }
    virtual const EventRegionSet& GetRegionsToProcess() const { return AllEventRegions; }
    virtual const EventSubCategorySet GetSubCategoriesToProcess() const { return AllEventSubCategories; }

    MetaIdSet CreateMetaIdSet(const ntuple::Flat& event, const FlatEventInfo::BjetPair& bjet_pair)
    {
        using namespace cuts::Htautau_Summer13::btag;

        MetaIdSet result;

        const EventEnergyScale energyScale = static_cast<EventEnergyScale>(event.eventEnergyScale);
        if(!GetEnergyScalesToProcess().count(energyScale)) return result;

        const EventCategoryVector categories =
                DetermineEventCategories(event.csv_Bjets, bjet_pair, event.nBjets_retagged, CSVL, CSVM, do_retag);
        for(EventCategory category : categories) {
            if(!GetCategoriesToProcess().count(category)) continue;
            const EventRegion region = DetermineEventRegion(event, category);
            if(!GetRegionsToProcess().count(region)) continue;

            const auto subCategories = DetermineEventSubCategories(GetFlatEventInfo(event, bjet_pair));
            for(EventSubCategory subCategory : subCategories) {
                if(!GetSubCategoriesToProcess().count(subCategory)) continue;
                result.insert(MetaId(category, subCategory, region, energyScale));
            }
        }

        return result;
    }

    static EventSubCategorySet DetermineEventSubCategories(const FlatEventInfo& eventInfo)
    {
        using namespace cuts::massWindow;

        EventSubCategorySet result;
        result.insert(EventSubCategory::NoCuts);

        const double mass_tautau = eventInfo.event->m_sv_MC;
        const bool inside_mass_window = mass_tautau > m_tautau_low && mass_tautau < m_tautau_high
                && eventInfo.Hbb.M() > m_bb_low && eventInfo.Hbb.M() < m_bb_high;

        if(inside_mass_window)
            result.insert(EventSubCategory::MassWindow);
        else
            result.insert(EventSubCategory::OutsideMassWindow);

        if(eventInfo.fitResults.has_valid_mass)
            result.insert(EventSubCategory::KinematicFitConverged);
        if(eventInfo.fitResults.has_valid_mass && inside_mass_window)
            result.insert(EventSubCategory::KinematicFitConvergedWithMassWindow);
        if(eventInfo.fitResults.has_valid_mass && !inside_mass_window)
            result.insert(EventSubCategory::KinematicFitConvergedOutsideMassWindow);

        return result;
    }

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

    std::shared_ptr<TFile> GetOutputFile() { return outputFile; }

    static EventRegion DetermineEventRegion(const ntuple::Flat& event, EventCategory category)
    {
        const Channel channel = static_cast<analysis::Channel>(event.channel);

        if (channel == Channel::MuTau){
            using namespace cuts::Htautau_Summer13::MuTau;

            if(!event.againstMuonTight_2
                    || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                    || (event.pfRelIso_1 >= muonID::pFRelIso && !IsAntiIsolatedRegion(event))
                    || (event.mt_1 >= muonID::mt && !IsHighMtRegion(event,category)))
                return EventRegion::Unknown;

            const bool os = event.q_1 * event.q_2 == -1;
            const bool iso = event.pfRelIso_1 < muonID::pFRelIso;
            const bool low_mt = event.mt_1 < muonID::mt;


            if(iso && os) return low_mt ? EventRegion::OS_Isolated : EventRegion::OS_Iso_HighMt;
            if(iso && !os) return low_mt ? EventRegion::SS_Isolated : EventRegion::SS_Iso_HighMt;
            if(os) return low_mt ? EventRegion::OS_AntiIsolated : EventRegion::OS_AntiIso_HighMt;
            return low_mt ? EventRegion::SS_AntiIsolated : EventRegion::SS_AntiIso_HighMt;
        }
        if (channel == Channel::ETau){

            using namespace cuts::Htautau_Summer13::ETau;

            if(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                    || (event.pfRelIso_1 >= electronID::pFRelIso && !IsAntiIsolatedRegion(event))
                    || (event.mt_1 >= electronID::mt && !IsHighMtRegion(event,category)) /*|| event.pt_2 <= 30*/ )
                return EventRegion::Unknown;

            const bool os = event.q_1 * event.q_2 == -1;
            const bool iso = event.pfRelIso_1 < electronID::pFRelIso;
            const bool low_mt = event.mt_1 < electronID::mt;


            if(iso && os) return low_mt ? EventRegion::OS_Isolated : EventRegion::OS_Iso_HighMt;
            if(iso && !os) return low_mt ? EventRegion::SS_Isolated : EventRegion::SS_Iso_HighMt;
            if(os) return low_mt ? EventRegion::OS_AntiIsolated : EventRegion::OS_AntiIso_HighMt;
            return low_mt ? EventRegion::SS_AntiIsolated : EventRegion::SS_AntiIso_HighMt;
        }
        if (channel == Channel::TauTau){

            using namespace cuts::Htautau_Summer13::TauTau::tauID;

            if(!event.againstElectronLooseMVA_2
                    || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= BackgroundEstimation::Isolation_upperLimit
                    || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= BackgroundEstimation::Isolation_upperLimit
                    || (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= byCombinedIsolationDeltaBetaCorrRaw3Hits
                        && event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= byCombinedIsolationDeltaBetaCorrRaw3Hits))
                return EventRegion::Unknown;

            const bool os = event.q_1 * event.q_2 == -1;
            const bool iso = event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                             event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < byCombinedIsolationDeltaBetaCorrRaw3Hits;

            if(iso) return os ? EventRegion::OS_Isolated : EventRegion::SS_Isolated;
            return os ? EventRegion::OS_AntiIsolated : EventRegion::SS_AntiIsolated;
        }
        throw exception("unsupported channel ") << channel;
    }

    const FlatEventInfo& GetFlatEventInfo(const ntuple::Flat& event, const FlatEventInfo::BjetPair& bjet_pair)
    {
        if(!eventInfoMap.count(bjet_pair))
            eventInfoMap[bjet_pair] = FlatEventInfoPtr(new FlatEventInfo(event, bjet_pair, recalc_kinfit));
        return *eventInfoMap.at(bjet_pair);
    }

private:
    std::shared_ptr<TFile> inputFile, outputFile;
    std::shared_ptr<ntuple::FlatTree> flatTree;
    FlatEventInfoMap eventInfoMap;

protected:
    bool recalc_kinfit;
    bool do_retag;
};

} // namespace analysis
