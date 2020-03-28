/*! Nonresonant EFT reweighing.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "h-tautau/McCorrections/include/WeightProvider.h"
#include "AnalysisTools/Core/include/SmartHistogram.h"
#include "NonResHH_EFT.h"

namespace analysis {
namespace NonResHH_EFT {

class WeightProvider : public ::analysis::mc_corrections::IWeightProvider {
public:
    using Hist = TH2D;
    using SmartHist = ::root_ext::SmartHistogram<Hist>;
    using HistPtr = std::shared_ptr<SmartHist>;
    using ParamMap = std::map<int, std::map<int, GF_Parameterization>>;

    static constexpr bool enable_negative_weight_warning = false;

    static const std::string& PangeaName();
    static const std::string& AllEventTupleName();

    WeightProvider(const std::string& param_cfg_name, TFile* file = nullptr);
    void AddFile(TFile& file);
    void CreatePdfs(TFile* file = nullptr);
    void SetTargetPoint(const Point& _point);

    virtual double Get(EventInfo& eventInfo) const override;
    virtual double Get(const ntuple::ExpressEvent& event) const override;

    template<typename Event>
    double Get(const Event& event, const Point& _point)
    {
        SetTargetPoint(_point);
        return Get(event.lhe_hh_m, event.lhe_hh_cosTheta);
    }

private:
    double Get(double lhe_hh_m, double lhe_hh_cosTheta) const;
    void CheckOverflows() const;
    void ReadParameterizationConfig(const std::string& param_cfg_name);

private:
    HistPtr pangea, pangea_pdf, sm_pdf, sm_pangea_ratio;
    bool has_pdf{false};
    Point point;
    ParamMap parameterization;
    double param_denom;
};

} // namespace NonResHH_EFT
} // namespace analysis
