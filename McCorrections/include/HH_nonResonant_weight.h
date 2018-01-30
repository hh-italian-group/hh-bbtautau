/*! Nonresonant EFT reweighing.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/McCorrections/include/WeightProvider.h"
#include "AnalysisTools/Core/include/SmartHistogram.h"
#include "hh-bbtautau/Analysis/include/NonResHH_EFT.h"

namespace analysis {
namespace NonResHH_EFT {

class WeightProvider : public ::analysis::mc_corrections::IWeightProvider {
public:
    using Hist = TH2D;
    using SmartHist = ::root_ext::SmartHistogram<Hist>;
    using HistPtr = std::shared_ptr<SmartHist>;
    using ParamMap = std::map<int, std::map<int, GF_Parameterization>>;

    static constexpr bool enable_negative_weight_warning = false;

    static const std::string& PangeaName() { static const std::string name = "pangea"; return name; }
    static const std::string& AllEventTupleName() { static const std::string name = "pangea"; return name; }

    WeightProvider(const std::string& param_cfg_name, TFile* file = nullptr)
    {
        static const std::vector<double> mhh_Bins = { 250, 260, 270, 280, 290, 300, 310, 320, 330, 340,
                                                      350, 360, 370, 380, 390, 400, 410, 420, 430, 440,
                                                      450, 460, 470, 480, 490, 500, 510, 520, 530, 540,
                                                      550, 600, 610, 620, 630, 640, 650, 660, 670, 680,
                                                      690, 700, 750, 800, 850, 900, 950, 1000, 1100, 1200,
                                                      1300, 1400, 1500, 1750, 2000, 50000 };
        static const std::vector<double> cosTheta_Bins = { 0.0, 0.4, 0.6, 0.8, 1.0 };
        pangea = std::make_shared<SmartHist>(PangeaName(), mhh_Bins, cosTheta_Bins);
        pangea_pdf = std::make_shared<SmartHist>(PangeaName() + "_pdf", mhh_Bins, cosTheta_Bins);
        sm_pdf = std::make_shared<SmartHist>("sm_pdf", mhh_Bins, cosTheta_Bins);
        sm_pangea_ratio = std::make_shared<SmartHist>("sm_" + PangeaName() + "_ratio", mhh_Bins, cosTheta_Bins);
        ReadParameterizationConfig(param_cfg_name);
        SetTargetPoint(Point::SM_Point());

        if(file) {
            AddFile(*file);
            CreatePdfs();
        }
    }

    void AddFile(TFile& file)
    {
        auto f_pangea = root_ext::TryReadObject<Hist>(file, PangeaName());
        if(f_pangea) {
            pangea->Add(f_pangea);
        } else {
            ntuple::ExpressTuple all_events("all_events", &file, true);
            for(const auto& event : all_events)
                pangea->Fill(event.lhe_hh_m, std::abs(event.lhe_hh_cosTheta));
        }
        has_pdf = false;
    }

    void CreatePdfs(TFile* file = nullptr)
    {
        CheckOverflows();
        pangea_pdf->CopyContent(*pangea);
        RenormalizeHistogram(*pangea_pdf, 1, false);
        sm_pangea_ratio->CopyContent(*sm_pdf);
        sm_pangea_ratio->Divide(pangea_pdf.get());
        has_pdf = true;

        if(file) {
            root_ext::WriteObject(*pangea, file);
            root_ext::WriteObject(*pangea_pdf, file);
            root_ext::WriteObject(*sm_pdf, file);
            root_ext::WriteObject(*sm_pangea_ratio, file);
        }
    }

    void SetTargetPoint(const Point& _point)
    {
        point = _point;
        param_denom = GF_Parameterization::Denominator_13TeV().Evaluate(point);
    }

    virtual double Get(const ntuple::Event& event) const override { return GetT(event); }
    virtual double Get(const ntuple::ExpressEvent& event) const override { return GetT(event); }

    template<typename Event>
    double Get(const Event& event, const Point& _point)
    {
        SetTargetPoint(_point);
        return Get(event);
    }

private:
    template<typename Event>
    double GetT(const Event& event) const
    {
        if(!has_pdf)
            throw exception("Pangea pdf is not created.");

        double m_hh = event.lhe_hh_m;
        double cos_theta = std::abs(event.lhe_hh_cosTheta);
        const Int_t bin_x = pangea_pdf->GetXaxis()->FindBin(m_hh);
        const Int_t bin_y = pangea_pdf->GetYaxis()->FindBin(cos_theta);
        if(bin_x < 1 || bin_x > pangea_pdf->GetNbinsX() || bin_y < 1 || bin_y > pangea_pdf->GetNbinsY())
            return 0;

        const double pdf_ratio = sm_pangea_ratio->GetBinContent(bin_x, bin_y);
        const double param_ratio = parameterization.at(bin_x).at(bin_y).Evaluate(point) / param_denom;
        const double weight = pdf_ratio * param_ratio;
        if(weight < 0) {
            if(enable_negative_weight_warning) {
                static std::set<std::tuple<int, int, int, int, int, int, int>> reported_bins;
                const auto bin_tuple = std::make_tuple(int(point.kt * 1000), int(point.kl * 1000), int(point.c2 * 1000),
                                                       int(point.cg * 1000), int(point.c2g * 1000), bin_x, bin_y);
                if(!reported_bins.count(bin_tuple)) {
                    std::cerr << boost::format("\nWarning: invalid EFT weight for m_hh = %1%, cos_theta = %2% in bin "
                                 "(%3%, %4%): n_pangea = %5%, pdf_pangea = %6% +/- %7%, pdf_sm = %8% +/- %9%, "
                                 "pdf_ratio = %10%, param_num = %11%, param_denom = %12%, param_ratio = %13%, "
                                 "weight = %14%.\n(kl kt c2 cg c2g) = (%15%).\nA[n] = [ %16% ]")
                                 % m_hh % cos_theta % bin_x % bin_y % pangea->GetBinContent(bin_x, bin_y)
                                 % pangea_pdf->GetBinContent(bin_x, bin_y) % pangea_pdf->GetBinError(bin_x, bin_y)
                                 % sm_pdf->GetBinContent(bin_x, bin_y) % sm_pdf->GetBinError(bin_x, bin_y) % pdf_ratio
                                 % parameterization.at(bin_x).at(bin_y).Evaluate(point) % param_denom % param_ratio
                                 % weight % point % parameterization.at(bin_x).at(bin_y) << std::endl;
                    reported_bins.insert(bin_tuple);
                }
            }

            return 0;
        }
        return weight;
    }

    void CheckOverflows() const
    {
        const Int_t N = pangea->GetNbinsX() + 1;
        const Int_t H = pangea->GetNbinsY() + 1;
        for (Int_t n = 1; n <= N; ++n){
            for (Int_t h = 0; h <= H; ++h){
                const bool zero_content = pangea->GetBinContent(n,h) == 0;
                const bool overflow_bin = n == 0 || n == N || h == 0 || h == H;
                if((zero_content && !overflow_bin) || (!zero_content && overflow_bin)) {
                    std::ostringstream ss;
                    const std::string prefix = overflow_bin ? "Non empty" : "Empty";
                    ss << prefix << " bin in HH nonResonant pangea: (" << n << ", " << h << ") with center at ("
                       << pangea->GetXaxis()->GetBinCenter(n) << ", " << pangea->GetYaxis()->GetBinCenter(h) << ").";
                    throw exception(ss.str());
                }
            }
        }
    }

    void ReadParameterizationConfig(const std::string& param_cfg_name)
    {
        std::ifstream cfg(param_cfg_name);
        if(cfg.fail())
            throw exception("Failed to open config file '%1%'.") % param_cfg_name;

        size_t n = 0;
        while(cfg.good()) {
            std::string line;
            std::getline(cfg, line);
            ++n;
            if(line.empty()) continue;
            std::istringstream ss(line);

            try {
                int id;
                double m_hh, cos_theta, n_sm, n_bsm;
                GF_Parameterization params;
                ss >> id >> m_hh >> cos_theta >> n_sm >> n_bsm >> params;
                if(ss.fail())
                    throw exception("Failed to parse '%1%'.") % line;

                const int bin_x = sm_pdf->GetXaxis()->FindBin(m_hh);
                const int bin_y = sm_pdf->GetYaxis()->FindBin(cos_theta);
                if(bin_x < 1 || bin_y > sm_pdf->GetNbinsX() || bin_y < 1 || bin_y > sm_pdf->GetNbinsY())
                    throw exception("(m_hh, cos_theta) = (%1%, %2%) is out of range.") % m_hh % cos_theta;
                sm_pdf->SetBinContent(bin_x, bin_y, n_sm);
                sm_pdf->SetBinError(bin_x, bin_y, std::sqrt(n_sm));

                if(parameterization[bin_x].count(bin_y))
                    throw exception("Duplicated parametrization for bin (%1%, %2%) with center at (%3%, %4%).")
                        % bin_x % bin_y % sm_pdf->GetXaxis()->GetBinCenter(bin_x)
                        % sm_pdf->GetYaxis()->GetBinCenter(bin_y);
                parameterization[bin_x][bin_y] = params;
            } catch(std::exception& e) {
                throw exception("Error while paring line %1% of config '%2%'. %3%") % n % param_cfg_name % e.what();
            }
        }

        RenormalizeHistogram(*sm_pdf, 1, false);
        for(int bin_x = 1; bin_x <= sm_pdf->GetNbinsX(); ++bin_x) {
            for(int bin_y = 1; bin_y <= sm_pdf->GetNbinsY(); ++bin_y) {
                if(!parameterization[bin_x].count(bin_y))
                    throw exception("Missing parameterization for bin (%1%, %2%) with center at (%3%, %4%).")
                        % bin_x % bin_y % sm_pdf->GetXaxis()->GetBinCenter(bin_x)
                        % sm_pdf->GetYaxis()->GetBinCenter(bin_y);
            }
        }
    }

private:
    HistPtr pangea, pangea_pdf, sm_pdf, sm_pangea_ratio;
    bool has_pdf{false};
    Point point;
    ParamMap parameterization;
    double param_denom;
};

} // namespace NonResHH_EFT
} // namespace analysis
