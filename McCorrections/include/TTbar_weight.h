/*! The ttbar weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "../../Instruments/include/TTFileDescriptor.h"


namespace analysis {
namespace mc_corrections {

class TTbar_weight {
public:
    using Event = ntuple::Event;

    TTbar_weight(const std::string& sm_weight_file_name) :
    {
        ttbar_descriptors = analysis::sample_merging::TTBinDescriptor::LoadConfig(sm_weight_file_name);
        GenEventType genEventType = static_cast<GenEventType>(event.genEventType);
    }

    double Get(const Event& event)
    {
        double m_hh = event.lhe_hh_m;
        double cos_Theta = event.lhe_hh_cosTheta;
        const Int_t bin_x = sm_weight->GetXaxis()->FindBin(m_hh);
        const Int_t bin_y = sm_weight->GetYaxis()->FindBin(cos_Theta);

        return sm_weight->GetBinContent(bin_x,bin_y);
    }

private:
    std::vector<analysis::sample_merging::TTBinDescriptor> ttbar_descriptors;
    //std::map<> mappa di genEvent o del range??


};

} // namespace mc_corrections
} // namespace analysis
