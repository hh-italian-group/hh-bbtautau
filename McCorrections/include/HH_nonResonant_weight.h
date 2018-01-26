/*! The sm weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/McCorrections/include/WeightProvider.h"
#include "AnalysisTools/Core/include/SmartHistogram.h"

namespace analysis {
namespace NonResHH_EFT {

class WeightProvider : public ::analysis::mc_corrections::IWeightProvider {
public:
    using Hist = TH2D;
    using SmartHist = ::root_ext::SmartHistogram<Hist>;
    using HistPtr = std::shared_ptr<SmartHist>;

    static const std::string& PangeaName() { static const std::string name = "pangea"; return name; }

    WeightProvider()
    {
        static const std::vector<double> mhh_Bins = { 250, 260, 270, 280, 290, 300, 310, 320, 330, 340,
                                                      350, 360, 370, 380, 390, 400, 410, 420, 430, 440,
                                                      450, 460, 470, 480, 490, 500, 510, 520, 530, 540,
                                                      550, 600, 610, 620, 630, 640, 650, 660, 670, 680,
                                                      690, 700, 750, 800, 850, 900, 950, 1000, 1100, 1200,
                                                      1300, 1400, 1500, 1750, 2000, 50000 };
        static const std::vector<double> cosTheta_Bins = { 0.0, 0.4, 0.6, 0.8, 1.0 };
        pangea = std::make_shared<SmartHist>(PangeaName(), mhh_Bins, cosTheta_Bins);
    }

    void AddFile(TFile& file)
    {
        auto f_pangea = root_ext::ReadObject<Hist>(file, "pangea");
    }

    void WritePangea(TFile& file)
    {
        root_ext::WriteObject(*pangea, &file, "pangea");
    }

    virtual double Get(const ntuple::Event& event) const override { return GetT(event); }
    virtual double Get(const ntuple::ExpressEvent& event) const override { return GetT(event); }

private:
    template<typename Event>
    double GetT(const Event& event) const
    {
        double m_hh = event.lhe_hh_m;
        double cos_Theta = event.lhe_hh_cosTheta;
        const Int_t bin_x = sm_weight->GetXaxis()->FindBin(m_hh);
        const Int_t bin_y = sm_weight->GetYaxis()->FindBin(std::abs(cos_Theta));
        if(bin_x < 1 || bin_x > sm_weight->GetNbinsX() || bin_y < 1 || bin_y > sm_weight->GetNbinsY())
            throw exception("Unable to estimate HH BSM to SM weight for the event with m_hh = %1%"
                            " and cos(theta) = %2%.") % m_hh % cos_Theta;

        return sm_weight->GetBinContent(bin_x,bin_y);
    }

private:
    static HistPtr LoadSMweight(const std::string& sm_weight_file_name, const std::string& hist_name)
    {
        auto inputFile_weight = root_ext::OpenRootFile(sm_weight_file_name);
        return HistPtr(root_ext::ReadCloneObject<Hist>(*inputFile_weight, hist_name, "", true));
    }

private:
    HistPtr pangea_pdf;
};

} // namespace mc_corrections
} // namespace analysis
