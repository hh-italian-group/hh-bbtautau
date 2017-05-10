/*! The ttbar weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "../../Instruments/include/TTFileDescriptor.h"


namespace analysis {
namespace mc_corrections {

class TTbar_pt_weight {
public:
    using Event = ntuple::ExpressEvent;

    TTbar_pt_weight() : {}

    double Get(const Event& event)
    {
        return std::sqrt(std::exp(0.0615 - 0.0005 * event.gen_top_pt) * std::exp(0.0615 - 0.0005*event.gen_topBar_pt));
    }



};

} // namespace mc_corrections
} // namespace analysis
