/*! Various hh event weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/McCorrections/include/EventWeights.h"

namespace analysis {
namespace mc_corrections {

class EventWeights_HH : public EventWeights {
public:
    EventWeights_HH(Period period, JetOrdering jet_ordering, DiscriminatorWP btag_wp, bool use_LLR_weights,
                    WeightingMode mode = {});

    ntuple::ProdSummary GetSummaryWithWeights(const std::shared_ptr<TFile>& file,
                                              const WeightingMode& weighting_mode) const;

private:
    static std::string FullBSMtoSM_Name(const std::string& fileName);
    static std::string Full_Cfg_Name(const std::string& fileName);
};

} // namespace mc_corrections
} // namespace analysis
