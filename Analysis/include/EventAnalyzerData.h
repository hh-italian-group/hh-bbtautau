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

    const BinVector M_tt_Bins = { 10, 30, 50, 75, 90, 105, 120, 135, 150, 175, 210, 250, 300, 350, 400, 500 };
    const BinVector M_ttbb_Bins = { 250, 260, 270, 280, 290, 300, 325, 350, 375, 400, 425, 450, 475, 500, 550,
                                    600, 650, 700, 850, 1000 };
    const BinVector MT2_Bins = { 0, 25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 500 };
    const BinVector Pt_Bins = { 0, 20, 40, 60, 80, 100, 125, 150, 175, 200, 250, 300};

    TH1D_ENTRY_CUSTOM_EX(m_ttbb, M_ttbb_Bins, "M_{#tau#tau+jj} (GeV)", "dN/dm_{#tau#tau+jj} (1/GeV)", false, 1.5, true, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(m_ttbb_kinfit, M_ttbb_Bins, "M_{H}^{kinfit} (GeV)", "dN/dm_{H}^{kinfit} (1/GeV)", false, 1.4, true, true)
    TH1D_ENTRY_CUSTOM_EX(m_sv, M_tt_Bins, "M_{#tau#tau} (GeV)", "dN/dm_{#tau#tau} (1/GeV)", false, 1.6, true, true)
    TH1D_ENTRY_CUSTOM_EX(MT2, MT2_Bins, "MT2_{H} (GeV)", "dN/dm (1/GeV)", false, 1.4, true, true)
    TH1D_ENTRY_EX(mva_score, 40, -1, 1, "MVA score", "Events", true, 1.8, false, true)
    TH1D_ENTRY_EX(mt_tot, 25, 0, 500, "M_{T}^{TOT}(GeV)", "Events", false, 1.5, false, SaveAll)
    TH1D_ENTRY_EX(deta_hbbhtautau, 25, -5, 5, "#Delta#eta_{#tau#tau, jj}", "Events", false, 1.8, false, SaveAll)
    TH1D_ENTRY_EX(dphi_hbbhtautau, 20, -4, 4, "#Delta#phi(#tau#tau, jj)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(MT_htautau, 20, 0, 250, "M_{T}^{#tau#tau}(GeV)", "Events", false, 1.4, false, SaveAll)

    TH1D_ENTRY_CUSTOM_EX(m_tt_vis, M_tt_Bins, "M_{vis}(GeV)", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(pt_H_tt, Pt_Bins, "P_{T}(GeV)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(pt_H_tt_MET, Pt_Bins, "P_{T}(GeV)", "Events", false, 1.5, false, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(pt_1, Pt_Bins, "P_{T}(leading#tau_{h})(GeV)", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(eta_1, 20, -2.5, 2.5, "#eta(leading#tau_{h})", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(iso_1, 100, 0, 10, "Iso#tau_{1}", "Events", false, 1, false, SaveAll)
    TH1D_ENTRY_EX(mt_1, 20, 0, 200, "M_{T}1(GeV)", "Events", false, 1.3, false, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(pt_2, Pt_Bins, "P_{T}(subleading#tau_{h})(GeV)", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(eta_2, 20, -2.5, 2.5, "#eta(subleading#tau_{h})", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(iso_2, 100, 0, 10, "Iso#tau_{2}", "Events", false, 1, false, SaveAll)
    TH1D_ENTRY_EX(mt_2, 20, 0, 200, "M_{T}2(GeV)", "Events", false, 1.3, false, SaveAll)
    TH1D_ENTRY_EX(dR_l1l2, 20, 0, 4, "#DeltaR(#tau#tau)", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(abs_dphi_l1MET, 20, 0, 3.2, "|#Delta#phi(leading#tau_{h}, miss)|", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(dphi_htautauMET, 25, -4, 4, "#Delta#phi(#tau#tau, miss)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(dR_l1l2MET, 20, 0, 6, "#Delta R(leading#tau_{h}+subleading#tau_{h}, miss)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(dR_l1l2Pt_htautau, 20, 0, 500, "#DeltaR (leading#tau_{h},subleading#tau_{h}) P_{T}(#tau#tau) (GeV)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(mass_l1l2MET, M_tt_Bins, "M_{leading#tau_{h}+subleading#tau_{h}+miss}", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(pt_l1l2MET, Pt_Bins, "P_{T}(leading#tau_{h}+subleading#tau_{h}+miss) (GeV)", "Events", false, 1.4, false, SaveAll)

    TH1D_ENTRY_EX(npv, 40, 0, 40,  "Number of Primary Vertex", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(MET, 25, 0, 200, "E_{T}^{miss}(GeV)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(phiMET, 20, 0, 3.2, "#phi_{E_{T}^{miss}}", "Events", false, 1.5, false, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(pt_MET, Pt_Bins, "P_{T}^{miss} (GeV)", "Events", false, 1.6, false, SaveAll)

    TH1D_ENTRY_CUSTOM_EX(m_bb, M_tt_Bins, "M_{jj} (GeV)", "Events/ 20 GeV", false, 1.6, true, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(pt_H_bb, Pt_Bins, "P_{T}(GeV)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(pt_b1, Pt_Bins, "Leading selected jet p_{T} (GeV)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(eta_b1, 20, -2.5, 2.5, "Leading selected jet #eta", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(csv_b1, 20, 0.8, 1, "Leading selected jet CSV", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(pt_b2, Pt_Bins, "Subleading selected jet p_{T} (GeV)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(eta_b2, 20, -2.5, 2.5, "Subleading selected jet #eta", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(csv_b2, 40, 0, 1, "Subleading selected jet CSV", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(costheta_METhbb, 20, -1, 1, "#cos #theta(jj, miss)", "Events", false, 1.4, false, SaveAll)

    TH1D_ENTRY_EX(dR_b1b2, 25, 0, 5, "#DeltaR_{jj}", "Events", false, 1.5, false, SaveAll)
    TH1D_ENTRY_EX(dR_b1b2_boosted, 25, 3, 8, "Boosted #Delta R_{jj} ", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(HT_otherjets, M_tt_Bins, "HT other jets (GeV)", "Events", false, 1.4, false, SaveAll)

    TH1D_ENTRY_CUSTOM_EX(mass_top1, M_tt_Bins, "mass top_1 (GeV)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(mass_top2, M_tt_Bins, "mass top_2 (GeV)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(p_zeta, 30,-50, 250, "p_{#zeta}", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(p_zetavisible, 25, 0, 200, "p_{#zeta}^{vis}", "Events", false, 1.4, false, SaveAll)

    BaseEventAnalyzerData(const EventCategory& eventCategory, bool _fill_all) :
        fill_all(_fill_all)
    {
        Initialize(eventCategory);
    }

    BaseEventAnalyzerData(std::shared_ptr<TFile> outputFile, const std::string& directoryName,
                          const EventCategory& eventCategory, bool _fill_all, bool readMode) :
        AnalyzerData(outputFile, directoryName, readMode), fill_all(_fill_all)
    {
        Initialize(eventCategory);
    }

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

        dR_l1l2().Fill(ROOT::Math::VectorUtil::DeltaR(t1.GetMomentum(),t2.GetMomentum()), weight);
        abs_dphi_l1MET().Fill(std::abs(ROOT::Math::VectorUtil::DeltaPhi(t1.GetMomentum(), event.GetMET().GetMomentum())), weight);
        dphi_htautauMET().Fill(ROOT::Math::VectorUtil::DeltaPhi(event.GetHiggsTTMomentum(true), event.GetMET().GetMomentum()),weight);
        dR_l1l2MET().Fill(ROOT::Math::VectorUtil::DeltaR(event.GetHiggsTTMomentum(false), event.GetMET().GetMomentum()),weight);
        dR_l1l2Pt_htautau().Fill(ROOT::Math::VectorUtil::DeltaR(t1.GetMomentum(),t2.GetMomentum())*event.GetHiggsTTMomentum(true).pt(), weight);
        mass_l1l2MET().Fill((event.GetHiggsTTMomentum(false)+event.GetMET().GetMomentum()).M(), weight);
        pt_l1l2MET().Fill((event.GetHiggsTTMomentum(false)+event.GetMET().GetMomentum()).pt(), weight);
        MT_htautau().Fill(Calculate_MT(event.GetHiggsTTMomentum(true), event.GetMET().GetMomentum()), weight);

        npv().Fill(event->npv,weight);
        MET().Fill(event.GetMET().GetMomentum().Pt(), weight);
        phiMET().Fill(event.GetMET().GetMomentum().Phi(),weight);

        pt_MET().Fill(event.GetMET().GetMomentum().pt(), weight);

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

        dphi_hbbhtautau().Fill(ROOT::Math::VectorUtil::DeltaPhi(Hbb.GetMomentum(), event.GetHiggsTTMomentum(true)), weight);
        deta_hbbhtautau().Fill((Hbb.GetMomentum()-event.GetHiggsTTMomentum(true)).Eta(), weight);
        mt_tot().Fill(Calculate_TotalMT(t1.GetMomentum(),t2.GetMomentum(),event.GetMET().GetMomentum()), weight);

        costheta_METhbb().Fill(four_bodies::Calculate_cosTheta_2bodies(event.GetMET().GetMomentum(), Hbb.GetMomentum()), weight);
        dR_b1b2().Fill(ROOT::Math::VectorUtil::DeltaR(b1.GetMomentum(), b2.GetMomentum()), weight);
        dR_b1b2_boosted().Fill(four_bodies::Calculate_dR_boosted(b1.GetMomentum(), b2.GetMomentum(), Hbb.GetMomentum()), weight);
        HT_otherjets().Fill(event->ht_other_jets, weight);

        mass_top1().Fill(four_bodies::Calculate_topPairMasses(t1.GetMomentum(), t2.GetMomentum(), b1.GetMomentum(), b2.GetMomentum(), event.GetMET().GetMomentum()).first, weight);
        mass_top2().Fill(four_bodies::Calculate_topPairMasses(t1.GetMomentum(), t2.GetMomentum(), b1.GetMomentum(), b2.GetMomentum(), event.GetMET().GetMomentum()).second, weight);

        p_zeta().Fill(Calculate_Pzeta(t1.GetMomentum(), t2.GetMomentum(),  event.GetMET().GetMomentum()), weight);
        p_zetavisible().Fill(Calculate_visiblePzeta(t1.GetMomentum(), t2.GetMomentum()), weight);

    }

private:
    void Initialize(const EventCategory& eventCategory)
    {
        static const BinVector boosted_MT2_bins = { 0, 150, 200, 250, 300, 500 };

        if(eventCategory.HasBoostConstraint() && eventCategory.IsBoosted()) {
            MT2.SetMasterHist(boosted_MT2_bins, "MT2_{H} (GeV)", "dN/dm (1/GeV)", false, 1.4, true, true);
        }
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
class EventAnalyzerData<ElectronCandidate, TauCandidate> : public BaseEventAnalyzerData {
public:
    using FirstLeg = ElectronCandidate;
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
        static const std::vector<double> res_mva_bins = { -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0,
                                                          0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1 };
        static const std::vector<double> boosted_mva_bins = { -1, 0.6, 0.8, 1 };

        auto mva_bins = &res_mva_bins;
        auto mva_draw_sf = 1.2;
        if(eventCategory.HasBoostConstraint() && eventCategory.IsBoosted()) {
            mva_bins = &boosted_mva_bins;
            mva_draw_sf = 2;
        }

        mva_score.SetMasterHist(*mva_bins, "MVA score", "dN / bin width", false, mva_draw_sf, true, true);
    }
};

template<>
class EventAnalyzerData<MuonCandidate, TauCandidate> : public BaseEventAnalyzerData {
public:
    using FirstLeg = MuonCandidate;
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
        static const std::vector<double> res_mva_bins = { -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0,
                                                          0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1 };
        static const std::vector<double> boosted_mva_bins = { -1, 0.6, 0.8, 0.9, 1 };

        auto mva_bins = &res_mva_bins;
        auto mva_draw_sf = 1.2;
        if(eventCategory.HasBoostConstraint() && eventCategory.IsBoosted()) {
            mva_bins = &boosted_mva_bins;
            mva_draw_sf = 2;
        }

        mva_score.SetMasterHist(*mva_bins, "MVA score", "dN / bin width", false, mva_draw_sf, true, true);
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
        static const BinVector res1b_mva_bins = { -1, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0,
                                                            0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1 };
        static const BinVector res2b_mva_bins = { -1, -0.8, -0.6, -0.3, 0.1, 0.4, 0.7, 0.85, 1 };
        static const BinVector boosted_mva_bins = { -1, 0, 0.6, 0.8, 1 };

        auto mva_bins = &res1b_mva_bins;
        auto mva_draw_sf = 1.4;
        if(eventCategory.HasBoostConstraint() && eventCategory.IsBoosted()) {
            mva_bins = &boosted_mva_bins;
            mva_draw_sf = 2;
        } else if(eventCategory.HasBtagConstraint() && eventCategory.N_btag() > 1) {
            mva_bins = &res2b_mva_bins;
        }

        mva_score.SetMasterHist(*mva_bins, "MVA score", "dN / bin width", false, mva_draw_sf, true, true);
    }
};

} // namespace analysis
