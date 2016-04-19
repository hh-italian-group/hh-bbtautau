/*! Study of different posibilities to select signal b-jets.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/LightBaseFlatTreeAnalyzer.h"
#include "hh-bbtautau/Analysis/include/FlatAnalyzerData.h"

class BjetSelectionStudyData : public root_ext::AnalyzerData {
public:
    BjetSelectionStudyData(std::shared_ptr<TFile> outputFile) : root_ext::AnalyzerData(outputFile) {}

    TH1D_ENTRY(MET, 35, 0, 350)
    TH1D_ENTRY(totalMatchedBjets, 3, -0.5, 2.5)
    TH1D_ENTRY(matchedBjets, 3, -0.5, 2.5)
    TH1D_ENTRY(matchedBjetsByPt, 3, -0.5, 2.5)
    TH1D_ENTRY(matchedBjetsByChi2, 3, -0.5, 2.5)
    TH1D_ENTRY(matchedBjets_byMassPair, 3, -0.5, 2.5)
    TH1D_ENTRY(Hbb_trueM, 30, 0, 300)
    TH1D_ENTRY(Hbb_truePt, 30, 0, 300)
    TH1D_ENTRY(Hbb_truePtOverM, 50, 0, 5)
    TH1D_ENTRY(Hbb_falseM, 30, 0, 300)
    TH1D_ENTRY(Hbb_falsePt, 30, 0, 300)
    TH1D_ENTRY(Hbb_falsePtOverM, 50, 0, 5)
    TH1D_ENTRY(Hbb_CSV_fail, 35, 0, 350)
    TH1D_ENTRY(Hbb_CSV_good, 35, 0, 350)
    TH1D_ENTRY(Chi2_CSVfail, 10, 0, 50)
    TH1D_ENTRY(Chi2_CSVgood, 10, 0, 50)
    TH1D_ENTRY(matchedBjets_combinedCSVandMASS, 3, -0.5, 2.5)
    TH1D_ENTRY(matchedBjets_CSVandBestMassPair, 3, -0.5, 2.5)
    TH1D_ENTRY(matchedBjets_combinedCSVandMASS_2, 3, -0.5, 2.5)
    TH1D_ENTRY(matchedBjets_CSVandBestMassPair_2, 3, -0.5, 2.5)
    TH2D_ENTRY(csv_b1_vs_ptb1, 20, 0, 200 ,25, 0.2, 1.2)
};

class BjetSelectionStudy : public analysis::LightBaseFlatTreeAnalyzer {
public:
    BjetSelectionStudy(const std::string& _inputFileName, const std::string& _outputFileName)
         : LightBaseFlatTreeAnalyzer(_inputFileName,_outputFileName), anaData(GetOutputFile())
    {
        recalc_kinfit = true;
        do_retag = false;
    }

protected:
    virtual PairSelectionMap SelectBjetPairs(const ntuple::Sync& event) override
    {
        PairSelectionMap pairMap;
        pairMap["CSV"] = SelectBestCsvPair(event);
        pairMap["Pt"] = SelectBestPtPair(event);
        pairMap["Chi2"] = SelectBestChi2Pair(event);
        pairMap["CSV_massWindow"] = SelectBestCsvPairWithMassWindow(event);
        pairMap["Pt_massWindow"] = SelectBestPtPairWithMassWindow(event);
        pairMap["Chi2_massWindow"] = SelectBestChi2PairWithMassWindow(event);
        return pairMap;
    }

    virtual const analysis::EventEnergyScaleSet GetEnergyScalesToProcess() const override
    {
        using analysis::EventEnergyScale;
        static const analysis::EventEnergyScaleSet energyScalesToProcess = { EventEnergyScale::Central };
        return energyScalesToProcess;
    }

    virtual const analysis::EventCategorySet& GetCategoriesToProcess() const override
    {
        using analysis::EventCategory;
        static const analysis::EventCategorySet categoriesToProcess = {
            EventCategory::TwoJets_Inclusive, EventCategory::TwoJets_ZeroBtag,
            EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag
        };

        return categoriesToProcess;
    }

    virtual const analysis::EventRegionSet& GetRegionsToProcess() const override
    {
        using analysis::EventRegion;
        static const analysis::EventRegionSet regionsToProcess = { EventRegion::OS_Isolated };
        return regionsToProcess;
    }

    virtual const analysis::EventSubCategorySet GetSubCategoriesToProcess() const override
    {
        using analysis::EventSubCategory;
        static const analysis::EventSubCategorySet subCategoriesToProcess = {
            EventSubCategory::NoCuts
            //EventSubCategory::KinematicFitConvergedWithMassWindow
        };
        return subCategoriesToProcess;
    }

    virtual void AnalyzeEvent(const analysis::SyncEventInfo& eventInfo, const MetaId& metaId,
                              const std::string& selectionLabel) override
    {
        std::ostringstream ss_label;
        ss_label << selectionLabel << "_" << metaId.eventCategory;
        const std::string label = ss_label.str();

        anaData.MET(label).Fill(eventInfo.MET.Pt());

        anaData.csv_b1_vs_ptb1(label).Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).Pt(),
                                           eventInfo.event->csv_bjets.at(eventInfo.selected_bjets.first));

        size_t totalMatchedBjets = 0;
        std::vector<size_t> indexes;
        TLorentzVector Hbb_true, Hbb_false;

        for (size_t k = 0; k < eventInfo.event->energy_bjets.size(); ++k) {
            if (eventInfo.event->isBjet_MC_Bjet.at(k)){
                Hbb_true += eventInfo.bjet_momentums.at(k);
                ++totalMatchedBjets;
                indexes.push_back(k);
            }
        }

        for (size_t h = 0; h < eventInfo.event->energy_Bjets.size(); ++h){
            for (size_t n = h+1; n < eventInfo.event->energy_Bjets.size(); ++n){
                if (eventInfo.event->isBjet_MC_Bjet.at(h) && eventInfo.event->isBjet_MC_Bjet.at(n)) continue;
                Hbb_false = eventInfo.bjet_momentums.at(h) + eventInfo.bjet_momentums.at(n);
                anaData.Hbb_falseM(label).Fill(Hbb_false.M());
                anaData.Hbb_falsePt(label).Fill(Hbb_false.Pt());
                anaData.Hbb_falsePtOverM(label).Fill(Hbb_false.Pt()/Hbb_false.M());
            }
        }


        if(totalMatchedBjets > 2)
            throw analysis::exception("Too many matched b-jets.");
        anaData.totalMatchedBjets(label).Fill(totalMatchedBjets);

        if (totalMatchedBjets < 2) return;

        anaData.Hbb_trueM(label).Fill(Hbb_true.M());
        anaData.Hbb_truePt(label).Fill(Hbb_true.Pt());
        anaData.Hbb_truePtOverM(label).Fill(Hbb_true.Pt()/Hbb_true.M());
        const size_t matchedBjets = N_SelectedMatchedBjets(*eventInfo.event, eventInfo.selected_bjets);
        anaData.matchedBjets(label).Fill(matchedBjets);

    }

private:
    static size_t N_SelectedMatchedBjets(const ntuple::Sync& event, const analysis::SyncEventInfo::BjetPair& bjet_pair)
    {
        size_t matchedBjets = 0;
        const std::vector<size_t> bjet_indexes = { bjet_pair.first, bjet_pair.second };
        for(size_t index : bjet_indexes) {
            if(index < event.isBjet_MC_Bjet.size() && event.isBjet_MC_Bjet.at(index))
                ++matchedBjets;
        }
        return matchedBjets;
    }

    template<typename Comparator>
    static std::vector<size_t> SortObjects(size_t n_objects, const Comparator& comparator)
    {
        std::vector<size_t> indexes(n_objects);
        std::iota(indexes.begin(), indexes.end(), 0);
        std::sort(indexes.begin(), indexes.end(), comparator);
        return indexes;
    }

    static analysis::SyncEventInfo::BjetPair DefaultPair() { return analysis::SyncEventInfo::BjetPair(0, 1); }

    template<typename Comparator>
    static analysis::SyncEventInfo::BjetPair SelectBestPair(const ntuple::Sync& event, const Comparator& comparator,
                                                            bool sort_pairs)
    {
        using analysis::FlatEventInfo;
        const size_t n_bjets = event.pt_Bjets.size();
        if(n_bjets < 2) return DefaultPair();
        const size_t n_objects = sort_pairs ? FlatEventInfo::NumberOfCombinationPairs(n_bjets) : n_bjets;
        const auto indexes = SortObjects(n_objects, comparator);
        return sort_pairs ? FlatEventInfo::CombinationIndexToPair(indexes.at(0), n_bjets)
                          : FlatEventInfo::BjetPair(indexes.at(0), indexes.at(1));
    }

    analysis::SyncEventInfo::BjetPair SelectBestCsvPair(const ntuple::Sync& event)
    {
        const auto comparator = [&] (size_t first, size_t second) -> bool
        {
            return event.csv_Bjets.at(first) > event.csv_Bjets.at(second);
        };

        return SelectBestPair(event, comparator, false);
    }

    analysis::SyncEventInfo::BjetPair SelectBestPtPair(const ntuple::Sync& event)
    {
        const auto comparator = [&] (size_t first, size_t second) -> bool
        {
            return event.pt_Bjets.at(first) > event.pt_Bjets.at(second);
        };

        return SelectBestPair(event, comparator, false);
    }

    analysis::SyncEventInfo::BjetPair SelectBestChi2Pair(const ntuple::Sync& event)
    {
        using analysis::FlatEventInfo;
        using analysis::kinematic_fit::four_body::FitResults;

        const size_t n_bjets = event.pt_Bjets.size();

        const auto comparator = [&] (size_t first, size_t second) -> bool
        {
            const auto first_pair = FlatEventInfo::CombinationIndexToPair(first, n_bjets);
            const auto second_pair = FlatEventInfo::CombinationIndexToPair(second, n_bjets);
            const FlatEventInfo& first_info = GetFlatEventInfo(event, first_pair);
            const FlatEventInfo& second_info = GetFlatEventInfo(event, second_pair);
            const FitResults& first_fit = first_info.fitResults;
            const FitResults& second_fit = second_info.fitResults;

            if(first_fit.has_valid_mass && !second_fit.has_valid_mass) return true;
            if(!first_fit.has_valid_mass) return false;
            return first_fit.chi2 < second_fit.chi2;
        };

        return SelectBestPair(event, comparator, true);
    }

    analysis::SyncEventInfo::BjetPair SelectBestCsvPairWithMassWindow(const ntuple::Sync& event)
    {
        using analysis::FlatEventInfo;
        using analysis::kinematic_fit::four_body::FitResults;
        using namespace cuts::massWindow;

        const size_t n_bjets = event.pt_Bjets.size();

        const auto comparator = [&] (size_t first, size_t second) -> bool
        {
            const auto first_pair = FlatEventInfo::CombinationIndexToPair(first, n_bjets);
            const auto second_pair = FlatEventInfo::CombinationIndexToPair(second, n_bjets);
            const FlatEventInfo& first_info = GetFlatEventInfo(event, first_pair);
            const FlatEventInfo& second_info = GetFlatEventInfo(event, second_pair);
            const FitResults& first_fit = first_info.fitResults;
            const FitResults& second_fit = second_info.fitResults;

            if(first_fit.has_valid_mass && !second_fit.has_valid_mass) return true;
            if(!first_fit.has_valid_mass) return false;

            const bool firstPair_inside_mass_window = first_info.Hbb.M() > m_bb_low && first_info.Hbb.M() < m_bb_high;
            const bool secondPair_inside_mass_window = second_info.Hbb.M() > m_bb_low && second_info.Hbb.M() < m_bb_high;

            if(firstPair_inside_mass_window && !secondPair_inside_mass_window) return true;
            if(!firstPair_inside_mass_window && secondPair_inside_mass_window) return false;
            const float_t csv_firstPair = event.csv_Bjets.at(first_pair.first) + event.csv_Bjets.at(first_pair.second);
            const float_t csv_secondPair = event.csv_Bjets.at(second_pair.first) + event.csv_Bjets.at(second_pair.second);
            return csv_firstPair > csv_secondPair;
        };

        return SelectBestPair(event, comparator, true);
    }

    analysis::SyncEventInfo::BjetPair SelectBestPtPairWithMassWindow(const ntuple::Sync& event)
    {

        using analysis::FlatEventInfo;
        using analysis::kinematic_fit::four_body::FitResults;
        using namespace cuts::massWindow;

        const size_t n_bjets = event.pt_Bjets.size();

        const auto comparator = [&] (size_t first, size_t second) -> bool
        {
            const auto first_pair = FlatEventInfo::CombinationIndexToPair(first, n_bjets);
            const auto second_pair = FlatEventInfo::CombinationIndexToPair(second, n_bjets);
            const FlatEventInfo& first_info = GetFlatEventInfo(event, first_pair);
            const FlatEventInfo& second_info = GetFlatEventInfo(event, second_pair);
            const FitResults& first_fit = first_info.fitResults;
            const FitResults& second_fit = second_info.fitResults;

            if(first_fit.has_valid_mass && !second_fit.has_valid_mass) return true;
            if(!first_fit.has_valid_mass) return false;

            const bool firstPair_inside_mass_window = first_info.Hbb.M() > m_bb_low && first_info.Hbb.M() < m_bb_high;
            const bool secondPair_inside_mass_window = second_info.Hbb.M() > m_bb_low && second_info.Hbb.M() < m_bb_high;

            if(firstPair_inside_mass_window && !secondPair_inside_mass_window) return true;
            if(!firstPair_inside_mass_window && secondPair_inside_mass_window) return false;
            const float_t pt_firstPair = event.pt_Bjets.at(first_pair.first) + event.pt_Bjets.at(first_pair.second);
            const float_t pt_secondPair = event.pt_Bjets.at(second_pair.first) + event.pt_Bjets.at(second_pair.second);
            return pt_firstPair > pt_secondPair;
        };

        return SelectBestPair(event, comparator, true);
    }

    analysis::SyncEventInfo::BjetPair SelectBestChi2PairWithMassWindow(const ntuple::Sync& event)
    {
        using analysis::FlatEventInfo;
        using analysis::kinematic_fit::four_body::FitResults;
        using namespace cuts::massWindow;

        const size_t n_bjets = event.pt_Bjets.size();

        const auto comparator = [&] (size_t first, size_t second) -> bool
        {
            const auto first_pair = FlatEventInfo::CombinationIndexToPair(first, n_bjets);
            const auto second_pair = FlatEventInfo::CombinationIndexToPair(second, n_bjets);
            const FlatEventInfo& first_info = GetFlatEventInfo(event, first_pair);
            const FlatEventInfo& second_info = GetFlatEventInfo(event, second_pair);
            const FitResults& first_fit = first_info.fitResults;
            const FitResults& second_fit = second_info.fitResults;

            if(first_fit.has_valid_mass && !second_fit.has_valid_mass) return true;
            if(!first_fit.has_valid_mass) return false;

            const bool firstPair_inside_mass_window = first_info.Hbb.M() > m_bb_low && first_info.Hbb.M() < m_bb_high;
            const bool secondPair_inside_mass_window = second_info.Hbb.M() > m_bb_low && second_info.Hbb.M() < m_bb_high;

            if(firstPair_inside_mass_window && !secondPair_inside_mass_window) return true;
            if(!firstPair_inside_mass_window && secondPair_inside_mass_window) return false;
            return first_fit.chi2 < second_fit.chi2;
        };

        return SelectBestPair(event, comparator, true);
    }

private:
    BjetSelectionStudyData anaData;
};
