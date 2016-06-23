/*! Definition of histogram containers for flat tree analyzers.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisCategories.h"

namespace analysis {

class BaseEventAnalyzerData : public root_ext::AnalyzerData {
public:
    TH1D_ENTRY_CUSTOM_EX(m_vis, M_tt_Bins2(), "M_{vis}(GeV)", "Events", false, 1.1, false, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(m_ttbb, M_ttbb_Bins(), "M_{#tau#tau+jj} (GeV)", "dN/dm_{#tau#tau+jj} (1/GeV)", false, 1.5, true, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(m_ttbb_log, M_ttbb_Bins(), "M_{#tau#tau+jj} (GeV)", "dN/dm_{#tau#tau+jj} (1/GeV)", true, 1., true, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(m_ttbb_kinfit, M_ttbb_Bins(), "M_{H}^{kinfit} (GeV)", "dN/dm_{H}^{kinfit} (1/GeV)", false, 1.4, true, true)


    TH1D_ENTRY_EX(npv, 40, 0, 40,  "Number of Primary Vertex", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(pt_b1, 25, 0, 500, "Leading selected jet p_{T} (GeV)", "Events", false, 1.4, false, SaveAll)
//    TH1D_ENTRY_EX(pt_b1_etaRange, 20, 0, 200, "Leading selected jet p_{T} (GeV), |#eta| < 1.5 ", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(eta_b1, 25, -2.5, 2.5, "Leading selected jet #eta", "Events", false, 1.8, false, SaveAll)
    TH1D_ENTRY_EX(csv_b1, 25, 0, 1, "Leading selected jet CSV", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(pt_b2, 25, 0, 500, "Subleading selected jet p_{T} (GeV)", "Events", false, 1.4, false, SaveAll)
//    TH1D_ENTRY_EX(pt_b2_etaRange, 20, 0, 200, "Subleading selected jet p_{T} (GeV), |#eta| < 1.5 ", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(eta_b2, 25, -2.5, 2.5, "Subleading selected jet #eta", "Events", false, 1.8, false, SaveAll)
    TH1D_ENTRY_EX(csv_b2, 25, 0, 1, "Subleading selected jet CSV", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(pt_H_tt, 20, 0, 300, "P_{T}(GeV)", "Events", false, 1.2, false, SaveAll)
    TH1D_ENTRY_EX(pt_H_bb, 20, 0, 300, "P_{T}(GeV)", "Events", false, 1.2, false, SaveAll)
    TH1D_ENTRY_EX(pt_H_hh, 20, 0, 300, "P_{T}(GeV)", "Events", false, 1.2, false, SaveAll)
    TH1D_ENTRY_EX(m_bb, 20, 0, 400, "M_{jj} (GeV)", "dN/dm_{jj} (1/GeV)", false, 1.4, true, SaveAll)
    TH1D_ENTRY_EX(DeltaPhi_tt, 22, 0., 3.3, "#Delta#Phi_{#tau#tau}[rad]", "Events", false, 1.3, false, SaveAll)
    TH1D_ENTRY_EX(DeltaPhi_bb, 22, 0., 3.3, "#Delta#Phi_{bb}[rad]", "Events", false, 1.8, false, SaveAll)
    TH1D_ENTRY_EX(DeltaPhi_bb_MET, 22, 0., 3.3, "#Delta#Phi_{bb,MET}[rad]", "Events", false, 1.5, false, SaveAll)
    TH1D_ENTRY_EX(DeltaPhi_tt_MET, 22, 0., 3.3, "#Delta#Phi_{#tau#tau,MET}[rad]", "Events", false, 1.5, false, SaveAll)
    TH1D_ENTRY_EX(DeltaPhi_hh, 22, 0., 3.3, "#Delta#Phi_{#tau#taubb}[rad]", "Events", false, 1.5, false, SaveAll)
    TH1D_ENTRY_EX(DeltaR_tt, 40, 0, 6, "#DeltaR_{#tau#tau}", "Events", false, 1.1, false, SaveAll)
    TH1D_ENTRY_EX(DeltaR_bb, 40, 0, 6, "#DeltaR_{bb}[rad]", "Events", false, 1.7, false, SaveAll)
    TH1D_ENTRY_EX(DeltaR_hh, 40, 0, 6, "#DeltaR_{#tau#taubb}[rad]", "Events", false, 1.5, false, SaveAll)
    TH1D_ENTRY_EX(mt_2, 20, 0, 200, "M_{T}(GeV)", "Events", false, 1.3, false, SaveAll)
    TH1D_ENTRY_EX(pt_H_tt_MET, 20, 0, 300, "P_{T}(GeV)", "Evnets", false, 1.1, false, SaveAll)
    TH1D_ENTRY_EX(convergence, 10, -3.5, 6.5, "Fit_convergence", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(chi2, 20, 0, 100, "#chi^{2}", "Events", false, 1.3, false, SaveAll)
    TH1D_ENTRY_EX(fit_probability, 20, 0, 1, "Fit_probability", "Events", false, 1.1, false, SaveAll)
    TH1D_ENTRY_EX(pull_balance, 20, -10, 10, "pull_balance", "Events", false, 2, false, SaveAll)
    TH1D_ENTRY_EX(pull_balance_1, 100, -10, 10, "pull_balance_1", "Events", false, 1.1, false, SaveAll)
    TH1D_ENTRY_EX(pull_balance_2, 100, -10, 10, "pull_balance_1", "Events", false, 1.1, false, SaveAll)
    TH1D_ENTRY_EX(MET, 30, 0, 150, "E_{T}^{miss} (GeV)", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(MET_wide, 40, 0, 400, "E_{T}^{miss} (GeV)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(phiMET, 30, 0, 3.2, "#phi_{E_{T}^{miss}}", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(nJets_Pt30, 5, 0, 4, "nJets", "Events", false, 1.1, false, SaveAll)
    TH2D_ENTRY_EX(csv_b1_vs_ptb1, 20, 0, 200, 25, 0, 1, "P_{T}(GeV)(leading_jet)", "CSV(leading_jet)", false, 1, SaveAll)
    TH2D_ENTRY_EX(chi2_vs_ptb1, 20, 0, 200, 20, 0, 100, "P_{T}(GeV)(leading_jet)", "#chi^{2}", false, 1, SaveAll)
    TH2D_ENTRY_EX(mH_vs_chi2, 20, 0, 100, 50, 200, 700, "#chi^{2}", "M_{#tau#taubb}(GeV)", false, 1, SaveAll)

    static constexpr bool SaveAll = true;

    virtual root_ext::SmartHistogram<TH1D>& m_sv_base() = 0;

    explicit BaseEventAnalyzerData(bool _fill_all) : fill_all(_fill_all) {}

    BaseEventAnalyzerData(std::shared_ptr<TFile> outputFile, const std::string& directoryName, bool _fill_all)
        : AnalyzerData(outputFile, directoryName), fill_all(_fill_all) {}

    using HistogramAccessor = root_ext::SmartHistogram<TH1D>& (BaseEventAnalyzerData::*)();

    virtual const std::vector<double>& M_tt_Bins() const
    {
        static const std::vector<double> bins = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150,
                                                  160, 170, 180, 190, 200, 225, 250, 275, 300, 325, 350 };
        return bins;
    }

    virtual const std::vector<double>& M_tt_Bins2() const
    {
        // static const std::vector<double> bins = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 240, 280, 320, 360 };
        static const std::vector<double> bins = { 10, 35, 60, 85, 110, 135, 160, 185, 210, 250, 300, 350, 400, 500 };
        return bins;
    }

    virtual const std::vector<double>& M_ttbb_Bins() const
    {
        static const std::vector<double> bins = { 260, 300, 350, 400, 450, 500, 600, 700, 850, 1000 };
        return bins;
    }

    virtual void FillBase(EventInfoBase& event, double weight)
    {
        if(event.HasBjetPair()){
            const double mX = event.GetResonanceMomentum(true, false).M();
            m_ttbb().Fill(mX, weight);
            m_ttbb_log().Fill(mX, weight);
            const auto& kinfit = event.GetKinFitResults();
            if(kinfit.HasValidMass())
                m_ttbb_kinfit().Fill(kinfit.mass, weight);
        }
        if(!fill_all) return;

        npv().Fill(event->npv,weight);
        m_vis().Fill(event.GetHiggsTTMomentum(false).M(),weight);
        mt_2().Fill(event->pfmt_2, weight);
        MET().Fill(event.GetMET().GetMomentum().Pt(), weight);
        MET_wide().Fill(event.GetMET().GetMomentum().Pt(), weight);
        phiMET().Fill(event.GetMET().GetMomentum().Phi(),weight);
        nJets_Pt30().Fill(event.GetNJets(),weight);
        if(!event.HasBjetPair()) return;

        const auto& Hbb = event.GetHiggsBB();
        const auto& b1 = Hbb.GetFirstDaughter();
        const auto& b2 = Hbb.GetSecondDaughter();
        pt_b1().Fill(b1.GetMomentum().pt(), weight);
        eta_b1().Fill(b1.GetMomentum().Eta(), weight);
        csv_b1().Fill(b1->csv(), weight);
        pt_b2().Fill(b2.GetMomentum().Pt(), weight);
        eta_b2().Fill(b2.GetMomentum().Eta(), weight);
        csv_b2().Fill(b2->csv(), weight);
    }

    virtual void CreateAll()
    {
        m_vis(); m_ttbb(); m_ttbb_log(); m_ttbb_kinfit(); pt_b1(); eta_b1(); csv_b1(); pt_b2(); eta_b2();
        csv_b2(); pt_H_tt(); pt_H_bb(); pt_H_hh(); /*m_bb();*/ DeltaPhi_tt(); DeltaPhi_bb(); DeltaPhi_bb_MET();
        DeltaPhi_tt_MET(); DeltaPhi_hh(); DeltaR_tt(); DeltaR_bb(); DeltaR_hh(); mt_2(); pt_H_tt_MET(); convergence();
        chi2(); fit_probability(); pull_balance(); pull_balance_1(); pull_balance_2(); MET(); MET_wide(); phiMET(); nJets_Pt30();
        csv_b1_vs_ptb1(); chi2_vs_ptb1(); mH_vs_chi2(); npv();
    }

protected:
    bool fill_all;
};

template<typename _FirstLeg>
class EventAnalyzerData : public BaseEventAnalyzerData {
public:
    using FirstLeg = _FirstLeg;
    using EventInfo = analysis::EventInfo<FirstLeg>;

    TH1D_ENTRY_CUSTOM_EX(m_sv, M_tt_Bins(), "M_{#tau#tau} (GeV)", "dN/dm_{#tau#tau} (1/GeV)", false, 1.3, true, true)
    TH1D_ENTRY_CUSTOM_EX(m_sv_bin, M_tt_Bins2(), "M_{#tau#tau} (GeV)", "dN/dm_{#tau#tau} (1/GeV)", false, 1.3, true, true)
    TH1D_ENTRY_EX(pt_1, 20, 0, 400, "P_{T}(leading#tau_{h})(GeV)", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(pt_1_log, 20, 0, 400, "Leading tau p_{T} (GeV)", "Events", true, 1.0, false, SaveAll)
    TH1D_ENTRY_EX(eta_1, 25, -2.5, 2.5, "#eta(leading#tau_{h})", "Events", false, 2, false, SaveAll)
    TH1D_ENTRY_EX(pt_2, 20, 0, 400, "P_{T}(subleading#tau_{h})(GeV)", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(pt_2_log, 20, 0, 400, "Subleading tau p_{T} (GeV)", "Events", true, 1.0, false, SaveAll)
    TH1D_ENTRY_EX(eta_2, 25, -2.5, 2.5, "#eta(subleading#tau_{h})", "Events", false, 2, false, SaveAll)

    TH1D_ENTRY_EX(mt_1, 20, 0, 200, "M_{T}(GeV)", "Events", false, 1.3, false, SaveAll)
    TH1D_ENTRY_EX(m_bb, 20, 0, 400, "M_{jj} (GeV)", "Events/ 20 GeV", false, 1.4, true, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(m_bb_bin, M_tt_Bins2(), "M_{jj} (GeV)", "Events", false, 1.4, true, SaveAll)

    using BaseEventAnalyzerData::BaseEventAnalyzerData;

    virtual void Fill(EventInfo& event, double weight)
    {
        BaseEventAnalyzerData::FillBase(event, weight);
        const double m_SVfit = event.GetHiggsTTMomentum(true).M();
        m_sv().Fill(m_SVfit, weight);
        m_sv_bin().Fill(m_SVfit, weight);
        if(!fill_all) return;

        pt_1().Fill(event.GetLeg(1).GetMomentum().pt(), weight);
        pt_1_log().Fill(event.GetLeg(1).GetMomentum().pt(), weight);
        eta_1().Fill(event.GetLeg(1).GetMomentum().eta(), weight);
        pt_2().Fill(event.GetLeg(2).GetMomentum().pt(), weight);
        pt_2_log().Fill(event.GetLeg(2).GetMomentum().pt(), weight);
        eta_2().Fill(event.GetLeg(2).GetMomentum().eta(), weight);
        mt_1().Fill(event.GetFirstLeg()->mt(MetType::PF), weight);
        if(!event.HasBjetPair()) return;

        const auto& Hbb = event.GetHiggsBB();
        m_bb().Fill(Hbb.GetMomentum().M(), weight);
        m_bb_bin().Fill(Hbb.GetMomentum().M(), weight);
    }

    virtual root_ext::SmartHistogram<TH1D>& m_sv_base() override { return m_sv(); }

    virtual const std::vector<double>& M_ttbb_Bins() const override
    {
        static const std::vector<double> bins = { 260, 300, 350, 400, 450, 500, 600, 700, 850, 1000 };
        return bins;
    }

    virtual void CreateAll() override
    {
        BaseEventAnalyzerData::CreateAll();
        m_sv(); m_sv_bin(); pt_1(); pt_1_log(); eta_1(); pt_2(); pt_2_log(); eta_2(); mt_1();m_bb();m_bb_bin();
    }
};

template<typename FirstLeg>
class EventAnalyzerData_semileptonic_2tag : public EventAnalyzerData<FirstLeg> {
public:
    using EventAnalyzerData<FirstLeg>::EventAnalyzerData;

    virtual const std::vector<double>& M_tt_Bins() const override
    {
        static const std::vector<double> bins = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 250, 300, 350 };
        return bins;
    }
};

template<>
class EventAnalyzerData<TauCandidate> : public BaseEventAnalyzerData {
public:
    using FirstLeg = TauCandidate;
    using EventInfo = analysis::EventInfo<FirstLeg>;

    TH1D_ENTRY_EX(pt_1, 20, 0, 200, "Leading tau p_{T} (GeV)", "Events", false, 2.0, false, SaveAll)
    TH1D_ENTRY_EX(eta_1, 25, -2.5, 2.5, "Leading tau #eta", "Events", false, 1.8, false, SaveAll)
    TH1D_ENTRY_EX(pt_2, 20, 0, 200, "Subleading tau p_{T} (GeV)", "Events", false, 2.0, false, SaveAll)
    TH1D_ENTRY_EX(eta_2, 25, -2.5, 2.5, "Subleading tau #eta", "Events", false, 1.8, false, SaveAll)

    TH1D_ENTRY_EX(mt_1, 20, 0, 200, "M_{T}(GeV)", "Events", false, 1.1, false, SaveAll)
    TH1D_ENTRY_EX(iso_tau1, 100, 0, 10, "Iso#tau_{1}", "Events", false, 1, false, SaveAll)
    TH1D_ENTRY_EX(iso_tau2, 100, 0, 10, "Iso#tau_{2}", "Events", false, 1, false, SaveAll)

    using BaseEventAnalyzerData::BaseEventAnalyzerData;

    virtual void Fill(EventInfo& event, double weight)
    {
        BaseEventAnalyzerData::FillBase(event, weight);
        if(!fill_all) return;

        const auto& tau1 = event.GetFirstLeg();
        const auto& tau2 = event.GetSecondLeg();
        pt_1().Fill(tau1.GetMomentum().pt(), weight);
        eta_1().Fill(tau1.GetMomentum().eta(), weight);
        pt_2().Fill(tau2.GetMomentum().pt(), weight);
        eta_2().Fill(tau2.GetMomentum().eta(), weight);
        mt_1().Fill(tau1->mt(MetType::PF), weight);
        iso_tau1().Fill(tau1->byCombinedIsolationDeltaBetaCorrRaw3Hits(), weight);
        iso_tau2().Fill(tau1->byCombinedIsolationDeltaBetaCorrRaw3Hits(), weight);
    }

    virtual const std::vector<double>& M_ttbb_Bins() const override
    {
        static const std::vector<double> bins = { 260, 300, 350, 400, 450, 500, 600, 700, 850, 1000 };
        return bins;
    }

    virtual void CreateAll() override
    {
        BaseEventAnalyzerData::CreateAll();
        pt_1(); eta_1(); pt_2(); eta_2(); mt_1(); iso_tau1(); iso_tau2();
    }
};

class EventAnalyzerData_tautau_other_tag : public EventAnalyzerData<TauCandidate> {
public:
    TH1D_ENTRY_CUSTOM_EX(m_sv, M_tt_Bins(), "M_{#tau#tau} (GeV)", "dN/dm_{#tau#tau} (1/GeV)", false, 1.5, true, true)
    TH1D_ENTRY_EX(m_bb, 20, 0, 400, "M_{jj} (GeV)", "Events/ 20 GeV", false, 1.4, true, SaveAll)

    using EventAnalyzerData::EventAnalyzerData;

    virtual void Fill(EventInfo& event, double weight) override
    {
        EventAnalyzerData::Fill(event, weight);
        const double m_SVfit = event.GetHiggsTTMomentum(true).M();
        m_sv().Fill(m_SVfit, weight);
        if(!event.HasBjetPair()) return;
        m_bb().Fill(event.GetHiggsBB().GetMomentum().M(), weight);
    }

    virtual root_ext::SmartHistogram<TH1D>& m_sv_base() override { return m_sv(); }

    virtual void CreateAll() override
    {
        EventAnalyzerData::CreateAll();
        m_sv();
        m_bb();
    }
};


class EventAnalyzerData_tautau_2tag : public EventAnalyzerData<TauCandidate> {
public:
    TH1D_ENTRY_CUSTOM_EX(m_sv, M_tt_Bins(), "M_{#tau#tau} (GeV)", "dN/dm_{#tau#tau} (1/GeV)", false, 1.7, true, true)
    TH1D_ENTRY_EX(m_bb, 10, 0, 400, "M_{jj} (GeV)", "Events/ 40 GeV", false, 1.4, true, SaveAll)

    using EventAnalyzerData::EventAnalyzerData;

    virtual void Fill(EventInfo& event, double weight) override
    {
        EventAnalyzerData::Fill(event, weight);
        const double m_SVfit = event.GetHiggsTTMomentum(true).M();
        m_sv().Fill(m_SVfit, weight);
        if(!event.HasBjetPair()) return;
        m_bb().Fill(event.GetHiggsBB().GetMomentum().M(), weight);
    }

    virtual root_ext::SmartHistogram<TH1D>& m_sv_base() override { return m_sv(); }

    virtual const std::vector<double>& M_tt_Bins() const override
    {
        static const std::vector<double> bins = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 250, 300, 350 };
        return bins;
    }

    virtual void CreateAll() override
    {
        EventAnalyzerData::CreateAll();
        m_sv();
        m_bb();
    }
};

} // namespace analysis
