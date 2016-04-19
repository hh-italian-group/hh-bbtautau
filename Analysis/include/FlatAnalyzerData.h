/*! Definition of histogram containers for flat tree analyzers.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "h-tautau/Analysis/include/FlatEventInfo.h"
#include "AnalysisCategories.h"

namespace analysis {

class FlatAnalyzerData : public root_ext::AnalyzerData {
public:
    TH1D_ENTRY_CUSTOM_EX(m_vis, M_tt_Bins2(), "M_{vis}(GeV)", "Events", false, 1.1, false, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(m_bb_slice, M_tt_bbSlice_Bins(), "2DM_{sv}(GeV)", "Events", false, 1.1, false, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(m_ttbb, M_ttbb_Bins(), "M_{#tau#tau+jj} (GeV)", "dN/dm_{#tau#tau+jj} (1/GeV)", false, 1.5, true, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(m_ttbb_kinfit, M_ttbb_Bins(), "M_{H}^{kinfit} (GeV)", "dN/dm_{H}^{kinfit} (1/GeV)", false, 1.4, true, true)


    TH1D_ENTRY_EX(npv, 40, 0, 40,  "Number of Primary Vertex", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(pt_b1, 25, 0, 500, "Leading selected jet p_{T} (GeV)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(pt_b1_etaRange, 20, 0, 200, "Leading selected jet p_{T} (GeV), |#eta| < 1.5 ", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(eta_b1, 25, -2.5, 2.5, "Leading selected jet #eta", "Events", false, 1.8, false, SaveAll)
    TH1D_ENTRY_EX(csv_b1, 25, 0, 1, "Leading selected jet CSV", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(pt_b2, 25, 0, 500, "Subleading selected jet p_{T} (GeV)", "Events", false, 1.4, false, SaveAll)
    TH1D_ENTRY_EX(pt_b2_etaRange, 20, 0, 200, "Subleading selected jet p_{T} (GeV), |#eta| < 1.5 ", "Events", false, 1.4, false, SaveAll)
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

    explicit FlatAnalyzerData(bool _fill_all) : fill_all(_fill_all) {}

    FlatAnalyzerData(std::shared_ptr<TFile> outputFile, const std::string& directoryName, bool _fill_all)
        : AnalyzerData(outputFile, directoryName), fill_all(_fill_all) {}

    typedef root_ext::SmartHistogram<TH1D>& (FlatAnalyzerData::*HistogramAccessor)();

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

    virtual const std::vector<double>& M_tt_bbSlice_Bins() const
    {
        static const size_t number_of_slices = 5;
        const static std::vector<double> bins = CreateMttSlicedBins(number_of_slices, M_tt_Bins());
        return bins;
    }

    virtual const std::vector<double>& M_ttbb_Bins() const
    {
        static const std::vector<double> bins = { 200, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 500, 550,
                                                  600, 650, 700 };
        return bins;
    }

    virtual void Fill(const SyncEventInfo& eventInfo, double weight){

        const ntuple::Sync& event = *eventInfo.event;

        if(!fill_all) return;

        npv().Fill(eventInfo.npv,weight);
      //  std::cout<<"  Filling number of primary vertex:   "<< eventInfo.npv<<std::endl;
        // DeltaPhi_tt().Fill(std::abs(eventInfo.lepton_momentums.at(0).DeltaPhi(eventInfo.lepton_momentums.at(1))), weight);
        // DeltaR_tt().Fill(eventInfo.lepton_momentums.at(0).DeltaR(eventInfo.lepton_momentums.at(1)), weight);
        // pt_H_tt().Fill(eventInfo.Htt.Pt(),weight);

        m_vis().Fill(eventInfo.Htt.M(),weight);
        // pt_H_tt_MET().Fill(eventInfo.Htt_MET.Pt(), weight);
        // DeltaPhi_tt_MET().Fill(std::abs(eventInfo.Htt.DeltaPhi(eventInfo.MET)), weight);
        mt_2().Fill(event.pfmt_2, weight);
        MET().Fill(eventInfo.MET.Pt(),weight);
        MET_wide().Fill(eventInfo.MET.Pt(),weight);
        phiMET().Fill(eventInfo.MET.Phi(),weight);
        nJets_Pt30().Fill(event.njets,weight);
        if(!eventInfo.has_bjet_pair) return;
        pt_b1().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).Pt(), weight);
        eta_b1().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).Eta(), weight);
        const bool jet1_etaRange = std::abs(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).Eta()) < 1.6 ;
        if (jet1_etaRange) pt_b1_etaRange().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).Pt(), weight);
        csv_b1().Fill(eventInfo.event->csv_jets.at(eventInfo.selected_bjets.first), weight);

        pt_b2().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second).Pt(), weight);
        const bool jet2_etaRange = std::abs(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second).Eta()) < 1.6 ;
        if (jet2_etaRange) pt_b2_etaRange().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second).Pt(), weight);
        eta_b2().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second).Eta(), weight);
        csv_b2().Fill(eventInfo.event->csv_jets.at(eventInfo.selected_bjets.second), weight);
        // DeltaPhi_bb().Fill(std::abs(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).DeltaPhi(
        //                                eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second))), weight);
        // DeltaR_bb().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).DeltaR(
        //                              eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second)), weight);
        // pt_H_bb().Fill(eventInfo.Hbb.Pt(),weight);
        m_ttbb().Fill(eventInfo.resonance.M(), weight);


    }

    virtual void CreateAll()
    {
        m_vis(); m_bb_slice(); m_ttbb(); m_ttbb_kinfit(); pt_b1(); eta_b1(); csv_b1(); pt_b2(); eta_b2();
        csv_b2(); pt_H_tt(); pt_H_bb(); pt_H_hh(); /*m_bb();*/ DeltaPhi_tt(); DeltaPhi_bb(); DeltaPhi_bb_MET();
        DeltaPhi_tt_MET(); DeltaPhi_hh(); DeltaR_tt(); DeltaR_bb(); DeltaR_hh(); mt_2(); pt_H_tt_MET(); convergence();
        chi2(); fit_probability(); pull_balance(); pull_balance_1(); pull_balance_2(); MET(); MET_wide(); phiMET(); nJets_Pt30();
        csv_b1_vs_ptb1(); chi2_vs_ptb1(); mH_vs_chi2(); npv();
    }

protected:
    static void FillSlice(TH1D& hist, double m_sv, double m_Hbb, double weight)
    {
        static const std::vector<double> slice_regions = { 60, 100, 140, 200, 600 };
        static const double slice_size = 350;

        if(m_sv < 0 || m_sv >= slice_size) return;
        const auto slice_region = std::find_if(slice_regions.begin(), slice_regions.end(),
                                               [&](double x) { return m_Hbb < x; });
        if(slice_region == slice_regions.end()) return;
        const ptrdiff_t slice_id = slice_region - slice_regions.begin();
        const double slice_shift = slice_size * slice_id;
        hist.Fill(m_sv + slice_shift, weight);
    }

    static std::vector<double> CreateMttSlicedBins(size_t number_of_slices, const std::vector<double>& mtt_bins)
    {
        if(!mtt_bins.size())
            throw exception("Invalid mtt bins.");
        if(!number_of_slices)
            throw exception("Invalid number of slices.");

        std::vector<double> sliced_bins;
        for(size_t n = 0; n < number_of_slices; ++n) {
            size_t k = n == 0 ? 0 : 1;
            const double shift = n * mtt_bins.back();
            for(; k < mtt_bins.size(); ++k) {
                const double value = mtt_bins.at(k) + shift;
                sliced_bins.push_back(value);
            }
        }
        return sliced_bins;
    }

    bool fill_all;
};

class FlatAnalyzerData_semileptonic : public FlatAnalyzerData {
public:
    TH1D_ENTRY_CUSTOM_EX(m_sv, M_tt_Bins(), "M_{#tau#tau} (GeV)", "dN/dm_{#tau#tau} (1/GeV)", false, 1.3, true, true)
    TH1D_ENTRY_CUSTOM_EX(m_sv_bin, M_tt_Bins2(), "M_{#tau#tau} (GeV)", "dN/dm_{#tau#tau} (1/GeV)", false, 1.3, true, true)
    TH1D_ENTRY_EX(pt_1, 20, 0, 400, "P_{T}(leading#tau_{h})(GeV)", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(pt_1_log, 20, 0, 400, "Leading tau p_{T} (GeV)", "Events", true, 2.0, false, SaveAll)
    TH1D_ENTRY_EX(eta_1, 25, -2.5, 2.5, "#eta(leading#tau_{h})", "Events", false, 2, false, SaveAll)
    TH1D_ENTRY_EX(pt_2, 20, 0, 400, "P_{T}(subleading#tau_{h})(GeV)", "Events", false, 1.6, false, SaveAll)
    TH1D_ENTRY_EX(pt_2_log, 20, 0, 400, "Subleading tau p_{T} (GeV)", "Events", true, 2.0, false, SaveAll)
    TH1D_ENTRY_EX(eta_2, 25, -2.5, 2.5, "#eta(subleading#tau_{h})", "Events", false, 2, false, SaveAll)

    TH1D_ENTRY_EX(mt_1, 20, 0, 200, "M_{T}(GeV)", "Events", false, 1.3, false, SaveAll)
    TH1D_ENTRY_EX(m_bb, 20, 0, 400, "M_{jj} (GeV)", "Events/ 20 GeV", false, 1.4, true, SaveAll)
    TH1D_ENTRY_CUSTOM_EX(m_bb_bin, M_tt_Bins2(), "M_{jj} (GeV)", "Events", false, 1.4, true, SaveAll)

    explicit FlatAnalyzerData_semileptonic(bool _fill_all) : FlatAnalyzerData(_fill_all) {}

    FlatAnalyzerData_semileptonic(std::shared_ptr<TFile> outputFile, const std::string& directoryName, bool fill_all)
        : FlatAnalyzerData(outputFile, directoryName, fill_all) {}

    virtual void Fill(const SyncEventInfo& eventInfo, double weight) override
    {
        FlatAnalyzerData::Fill(eventInfo, weight);
        m_sv().Fill(eventInfo.event->m_sv, weight);
        m_sv_bin().Fill(eventInfo.event->m_sv, weight);
        if(!fill_all) return;
        const ntuple::Sync& event = *eventInfo.event;
        pt_1().Fill(event.pt_1, weight);
        pt_1_log().Fill(event.pt_1, weight);
        eta_1().Fill(event.eta_1, weight);
        pt_2().Fill(event.pt_2, weight);
        pt_2_log().Fill(event.pt_2, weight);
        eta_2().Fill(event.eta_2, weight);
        mt_1().Fill(event.pfmt_1, weight);
        if(!eventInfo.has_bjet_pair) return;
        m_bb().Fill(eventInfo.Hbb.M(), weight);
        m_bb_bin().Fill(eventInfo.Hbb.M(), weight);
    }

    virtual root_ext::SmartHistogram<TH1D>& m_sv_base() override { return m_sv(); }

    virtual const std::vector<double>& M_ttbb_Bins() const override
    {
        static const std::vector<double> bins = { 100,200, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 500, 550,
                                                  600, 650, 700 };
        return bins;
    }

    virtual void CreateAll() override
    {
        FlatAnalyzerData::CreateAll();
        m_sv(); m_sv_bin(); pt_1(); pt_1_log(); eta_1(); pt_2(); pt_2_log(); eta_2(); mt_1();m_bb();m_bb_bin();
    }
};

class FlatAnalyzerData_semileptonic_2tag : public FlatAnalyzerData_semileptonic {
public:
    explicit FlatAnalyzerData_semileptonic_2tag(bool _fill_all) : FlatAnalyzerData_semileptonic(_fill_all) {}

    FlatAnalyzerData_semileptonic_2tag(std::shared_ptr<TFile> outputFile, const std::string& directoryName,
                                       bool fill_all)
        : FlatAnalyzerData_semileptonic(outputFile, directoryName, fill_all) {}

    virtual const std::vector<double>& M_tt_Bins() const override
    {
        static const std::vector<double> bins = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 250, 300, 350 };
        return bins;
    }
};

class FlatAnalyzerData_tautau : public FlatAnalyzerData {
public:
    TH1D_ENTRY_EX(pt_1, 20, 0, 200, "Leading tau p_{T} (GeV)", "Events", false, 2.0, false, SaveAll)
    TH1D_ENTRY_EX(eta_1, 25, -2.5, 2.5, "Leading tau #eta", "Events", false, 1.8, false, SaveAll)
    TH1D_ENTRY_EX(pt_2, 20, 0, 200, "Subleading tau p_{T} (GeV)", "Events", false, 2.0, false, SaveAll)
    TH1D_ENTRY_EX(eta_2, 25, -2.5, 2.5, "Subleading tau #eta", "Events", false, 1.8, false, SaveAll)

    TH1D_ENTRY_EX(mt_1, 20, 0, 200, "M_{T}(GeV)", "Events", false, 1.1, false, SaveAll)
    TH1D_ENTRY_EX(iso_tau1, 100, 0, 10, "Iso#tau_{1}", "Events", false, 1, false, SaveAll)
    TH1D_ENTRY_EX(iso_tau2, 100, 0, 10, "Iso#tau_{2}", "Events", false, 1, false, SaveAll)

    explicit FlatAnalyzerData_tautau(bool _fill_all) : FlatAnalyzerData(_fill_all) {}

    FlatAnalyzerData_tautau(std::shared_ptr<TFile> outputFile, const std::string& directoryName, bool fill_all)
        : FlatAnalyzerData(outputFile, directoryName, fill_all) {}

    virtual void Fill(const SyncEventInfo& eventInfo, double weight) override
    {
        FlatAnalyzerData::Fill(eventInfo, weight);
        if(!fill_all) return;
        const ntuple::Sync& event = *eventInfo.event;
        pt_1().Fill(event.pt_1, weight);
        eta_1().Fill(event.eta_1, weight);
        pt_2().Fill(event.pt_2, weight);
        eta_2().Fill(event.eta_2, weight);
        mt_1().Fill(event.mt_1, weight);
        iso_tau1().Fill(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1,weight);
        iso_tau2().Fill(event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2,weight);
    }

    virtual const std::vector<double>& M_ttbb_Bins() const override
    {
        static const std::vector<double> bins = { 200, 250, 280, 310, 340, 370, 400, 500, 600, 700 };
        return bins;
    }

    virtual void CreateAll() override
    {
        FlatAnalyzerData::CreateAll();
        pt_1(); eta_1(); pt_2(); eta_2(); mt_1(); iso_tau1(); iso_tau2();
    }
};

class FlatAnalyzerData_tautau_other_tag : public FlatAnalyzerData_tautau {
public:
    TH1D_ENTRY_CUSTOM_EX(m_sv, M_tt_Bins(), "M_{#tau#tau} (GeV)", "dN/dm_{#tau#tau} (1/GeV)", false, 1.5, true, true)
    TH1D_ENTRY_EX(m_bb, 20, 0, 400, "M_{jj} (GeV)", "Events/ 20 GeV", false, 1.4, true, SaveAll)

    explicit FlatAnalyzerData_tautau_other_tag(bool _fill_all) : FlatAnalyzerData_tautau(_fill_all) {}

    FlatAnalyzerData_tautau_other_tag(std::shared_ptr<TFile> outputFile, const std::string& directoryName, bool fill_all)
        : FlatAnalyzerData_tautau(outputFile, directoryName, fill_all) {}

    virtual void Fill(const SyncEventInfo& eventInfo, double weight) override
    {
        FlatAnalyzerData_tautau::Fill(eventInfo, weight);
        m_sv().Fill(eventInfo.event->m_sv, weight);
        if(!eventInfo.has_bjet_pair) return;
        m_bb().Fill(eventInfo.Hbb.M(), weight);
    }

    virtual root_ext::SmartHistogram<TH1D>& m_sv_base() override { return m_sv(); }

    virtual void CreateAll() override
    {
        FlatAnalyzerData_tautau::CreateAll();
        m_sv();
        m_bb();
    }
};


class FlatAnalyzerData_tautau_2tag : public FlatAnalyzerData_tautau {
public:
    TH1D_ENTRY_CUSTOM_EX(m_sv, M_tt_Bins(), "M_{#tau#tau} (GeV)", "dN/dm_{#tau#tau} (1/GeV)", false, 1.7, true, true)
    TH1D_ENTRY_EX(m_bb, 10, 0, 400, "M_{jj} (GeV)", "Events/ 40 GeV", false, 1.4, true, SaveAll)

    explicit FlatAnalyzerData_tautau_2tag(bool _fill_all) : FlatAnalyzerData_tautau(_fill_all) {}

    FlatAnalyzerData_tautau_2tag(std::shared_ptr<TFile> outputFile, const std::string& directoryName, bool fill_all)
        : FlatAnalyzerData_tautau(outputFile, directoryName, fill_all) {}

    virtual void Fill(const SyncEventInfo& eventInfo, double weight) override
    {
        FlatAnalyzerData_tautau::Fill(eventInfo, weight);
        m_sv().Fill(eventInfo.event->m_sv, weight);
        if(!eventInfo.has_bjet_pair) return;
        m_bb().Fill(eventInfo.Hbb.M(), weight);
    }

    virtual root_ext::SmartHistogram<TH1D>& m_sv_base() override { return m_sv(); }

    virtual const std::vector<double>& M_tt_Bins() const override
    {
        static const std::vector<double> bins = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 250, 300, 350 };
        return bins;
    }

    virtual void CreateAll() override
    {
        FlatAnalyzerData_tautau::CreateAll();
        m_sv();
        m_bb();
    }
};

} // namespace analysis
