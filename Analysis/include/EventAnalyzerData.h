/*! Definition of histogram containers for flat tree analyzers.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisCategories.h"

namespace analysis {

class BaseEventAnalyzerData : public root_ext::AnalyzerData {
public:    
    using BinVector = std::vector<double>;
    const bool SaveAll = true;

    const BinVector M_tt_Bins = { 10, 35, 60, 85, 110, 135, 160, 185, 210, 250, 300, 350, 400, 500 };
    const BinVector M_ttbb_Bins = { 250, 260, 270, 280, 290, 300, 325, 350, 375, 400, 425, 450, 475, 500, 550,
                                    600, 650, 700, 850, 1000 };
    const BinVector MT2_Bins = { 0, 25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 500 };

    TH1D_ENTRY_CUSTOM_EX(m_ttbb, M_ttbb_Bins, "M_{#tau#tau+jj} (GeV)", "dN/dm_{#tau#tau+jj} (1/GeV)", false, 1.5, true, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(m_ttbb_kinfit, M_ttbb_Bins, "M_{H}^{kinfit} (GeV)", "dN/dm_{H}^{kinfit} (1/GeV)", false, 1.4, true, true)
    TH1D_ENTRY_CUSTOM_EX(m_sv, M_tt_Bins, "M_{#tau#tau} (GeV)", "dN/dm_{#tau#tau} (1/GeV)", false, 1.3, true, true)
    TH1D_ENTRY_CUSTOM_EX(MT2, MT2_Bins, "MT2_{H} (GeV)", "dN/dm (1/GeV)", false, 1.4, true, true)
    TH1D_ENTRY_EX(mva_score, 40, -1, 1, "MVA score", "Events", true, 1.2, false, true)

    TH1D_ENTRY_CUSTOM_EX(m_tt_vis, M_tt_Bins, "M_{vis}(GeV)", "Events", false, 1.1, false, SaveAll)
    TH1D_ENTRY_EX(pt_H_tt, 20, 0, 300, "P_{T}(GeV)", "Events", false, 1.2, false, SaveAll)
    TH1D_ENTRY_EX(pt_H_tt_MET, 20, 0, 300, "P_{T}(GeV)", "Evnets", false, 1.1, false, SaveAll)
    TH1D_ENTRY_EX(pt_1, 20, 0, 400, "P_{T}(leading#tau_{h})(GeV)", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(eta_1, 25, -2.5, 2.5, "#eta(leading#tau_{h})", "Events", false, 2, false, SaveAll)
    TH1D_ENTRY_EX(iso_1, 100, 0, 10, "Iso#tau_{1}", "Events", false, 1, false, SaveAll)
    TH1D_ENTRY_EX(mt_1, 20, 0, 200, "M_{T}1(GeV)", "Events", false, 1.3, false, SaveAll)
    TH1D_ENTRY_EX(pt_2, 20, 0, 400, "P_{T}(subleading#tau_{h})(GeV)", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(eta_2, 25, -2.5, 2.5, "#eta(subleading#tau_{h})", "Events", false, 2, false, SaveAll)
    TH1D_ENTRY_EX(iso_2, 100, 0, 10, "Iso#tau_{2}", "Events", false, 1, false, SaveAll)
    TH1D_ENTRY_EX(mt_2, 20, 0, 200, "M_{T}2(GeV)", "Events", false, 1.3, false, SaveAll)

    TH1D_ENTRY_EX(npv, 40, 0, 40,  "Number of Primary Vertex", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(MET, 30, 0, 150, "E_{T}^{miss} (GeV)", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(phiMET, 30, 0, 3.2, "#phi_{E_{T}^{miss}}", "Events", false, 1.4, false, SaveAll)

    TH1D_ENTRY_EX(m_bb, 20, 0, 400, "M_{jj} (GeV)", "Events/ 20 GeV", false, 1.4, true, SaveAll)
    TH1D_ENTRY_EX(pt_H_bb, 20, 0, 300, "P_{T}(GeV)", "Events", false, 1.2, false, SaveAll)
    TH1D_ENTRY_EX(pt_b1, 25, 0, 500, "Leading selected jet p_{T} (GeV)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(eta_b1, 25, -2.5, 2.5, "Leading selected jet #eta", "Events", false, 1.8, false, SaveAll)
    TH1D_ENTRY_EX(csv_b1, 25, 0, 1, "Leading selected jet CSV", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(pt_b2, 25, 0, 500, "Subleading selected jet p_{T} (GeV)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(eta_b2, 25, -2.5, 2.5, "Subleading selected jet #eta", "Events", false, 1.8, false, SaveAll)
    TH1D_ENTRY_EX(csv_b2, 25, 0, 1, "Subleading selected jet CSV", "Events", false, 1.4, false, SaveAll)

    BaseEventAnalyzerData(const EventCategory& /*eventCategory*/, bool _fill_all) : fill_all(_fill_all) {}
    BaseEventAnalyzerData(std::shared_ptr<TFile> outputFile, const std::string& directoryName,
                          const EventCategory& /*eventCategory*/, bool _fill_all, bool readMode)
        : AnalyzerData(outputFile, directoryName, readMode), fill_all(_fill_all) {}

    void FillBase(EventInfoBase& event, double weight)
    {
        if(event.HasBjetPair()) {
            const double mX = event.GetResonanceMomentum(true, false).M();
            m_ttbb().Fill(mX, weight);
            const auto& kinfit = event.GetKinFitResults();
            if(kinfit.HasValidMass())
                m_ttbb_kinfit().Fill(kinfit.mass, weight);
            MT2().Fill(event.GetMT2(), weight);
            mva_score().Fill(event.GetMvaScore(), weight);
        }

        const double m_SVfit = event.GetHiggsTTMomentum(true).M();
        m_sv().Fill(m_SVfit, weight);

        if(!fill_all) return;

        const auto& Htt = event.GetHiggsTTMomentum(false);
        const auto& t1 = event.GetLeg(1);
        const auto& t2 = event.GetLeg(2);

        m_tt_vis().Fill(Htt.M(), weight);
        pt_H_tt().Fill(Htt.Pt(), weight);
        pt_H_tt_MET().Fill((Htt + event.GetMET().GetMomentum()).Pt(), weight);
        pt_1().Fill(t1.GetMomentum().pt(), weight);
        eta_1().Fill(t1.GetMomentum().eta(), weight);
        iso_1().Fill(t1.GetIsolation(), weight);
        mt_1().Fill(Calculate_MT(t1.GetMomentum(), event.GetMET().GetMomentum()), weight);
        pt_2().Fill(t2.GetMomentum().pt(), weight);
        eta_2().Fill(t2.GetMomentum().eta(), weight);
        iso_2().Fill(t2.GetIsolation(), weight);
        mt_2().Fill(Calculate_MT(t2.GetMomentum(), event.GetMET().GetMomentum()), weight);

        npv().Fill(event->npv,weight);
        MET().Fill(event.GetMET().GetMomentum().Pt(), weight);
        phiMET().Fill(event.GetMET().GetMomentum().Phi(),weight);

        if(!event.HasBjetPair()) return;

        const auto& Hbb = event.GetHiggsBB();
        const auto& b1 = Hbb.GetFirstDaughter();
        const auto& b2 = Hbb.GetSecondDaughter();
        m_bb().Fill(Hbb.GetMomentum().M(), weight);
        pt_H_bb().Fill(Hbb.GetMomentum().Pt(), weight);
        pt_b1().Fill(b1.GetMomentum().pt(), weight);
        eta_b1().Fill(b1.GetMomentum().Eta(), weight);
        csv_b1().Fill(b1->csv(), weight);
        pt_b2().Fill(b2.GetMomentum().Pt(), weight);
        eta_b2().Fill(b2.GetMomentum().Eta(), weight);
        csv_b2().Fill(b2->csv(), weight);
    }

protected:
    bool fill_all;
};

template<typename _FirstLeg, typename _SecondLeg>
class EventAnalyzerData : public BaseEventAnalyzerData {
public:
    using FirstLeg = _FirstLeg;
    using SecondLeg = _SecondLeg;
    using EventInfo = ::analysis::EventInfo<FirstLeg, SecondLeg>;

    using BaseEventAnalyzerData::BaseEventAnalyzerData;

    void Fill(EventInfo& event, double weight)
    {
        FillBase(event, weight);
    }
};

template<>
class EventAnalyzerData<TauCandidate, TauCandidate> : public BaseEventAnalyzerData {
public:
    using FirstLeg = TauCandidate;
    using SecondLeg = TauCandidate;
    using EventInfo = ::analysis::EventInfo<FirstLeg, SecondLeg>;

    EventAnalyzerData(const EventCategory& eventCategory, bool _fill_all) :
        BaseEventAnalyzerData(eventCategory, _fill_all)
    {
        Initialize(eventCategory);
    }

    EventAnalyzerData(std::shared_ptr<TFile> outputFile, const std::string& directoryName,
                          const EventCategory& eventCategory, bool _fill_all, bool readMode) :
        BaseEventAnalyzerData(outputFile, directoryName, eventCategory, _fill_all, readMode)
    {
        Initialize(eventCategory);
    }

    void Fill(EventInfo& event, double weight)
    {
        FillBase(event, weight);
    }

private:
    void Initialize(const EventCategory& eventCategory)
    {
        static const std::vector<double> res_mva_bins = { -1, 0.4, 0.55, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975, 1 };
        static const std::vector<double> boosted_mva_bins = { -1, 1 };
        const auto& mva_bins = eventCategory.HasBoostConstraint() && eventCategory.IsBoosted()
                             ? boosted_mva_bins : res_mva_bins;
        mva_score.SetMasterHist(mva_bins, "MVA score", "dN / bin width", true, 1.2, true, true);
    }
};


} // namespace analysis
