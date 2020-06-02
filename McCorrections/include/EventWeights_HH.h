/*! Various hh event weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/McCorrections/include/EventWeights.h"
#include "hh-bbtautau/McCorrections/include/NonResHH_EFT.h"

namespace analysis {
namespace mc_corrections {

class EventWeights_HH : public EventWeights {
public:
    EventWeights_HH(Period period, const BTagger& bTagger, const WeightingMode& mode = {});

    ntuple::ProdSummary GetSummaryWithWeights(const std::shared_ptr<TFile>& file, const WeightingMode& weighting_mode,
                                              const boost::optional<double>& max_gen_weight,
                                              bool control_duplicates = true) const;

    std::vector<double> GetTotalShapeWeights(const std::shared_ptr<TFile>& file, const WeightingMode& weighting_mode,
                                             const std::vector<NonResHH_EFT::Point>& eft_points, bool orthogonal);

private:
    static std::string FullBSMtoSM_Name(const std::string& fileName);
    static std::string Full_Cfg_Name(const std::string& fileName);
};

} // namespace mc_corrections
} // namespace analysis
