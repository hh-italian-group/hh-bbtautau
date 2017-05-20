/*! The sm weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventTuple.h"


namespace analysis {
namespace mc_corrections {

class HH_BMStoSM_weight {
public:
    using Event = ntuple::Event;
    using Hist = TH2;
    using HistPtr = std::shared_ptr<Hist>;

    HH_BMStoSM_weight(const std::string& sm_weight_file_name, const std::string& hist_name) :
        sm_weight(LoadSMweight(sm_weight_file_name, hist_name)) { }

    template<typename Event>
    double Get(const Event& event)
    {
        double m_hh = event.lhe_hh_m;
        double cos_Theta = event.lhe_hh_cosTheta;
        const Int_t bin_x = sm_weight->GetXaxis()->FindBin(m_hh);
        const Int_t bin_y = sm_weight->GetYaxis()->FindBin(cos_Theta);

        return sm_weight->GetBinContent(bin_x,bin_y);
    }

private:
    static HistPtr LoadSMweight(const std::string& sm_weight_file_name, const std::string& hist_name)
    {
        auto inputFile_weight = root_ext::OpenRootFile(sm_weight_file_name);
        return HistPtr(root_ext::ReadCloneObject<Hist>(*inputFile_weight, hist_name, "", true));
    }

private:
    HistPtr sm_weight;
};

} // namespace mc_corrections
} // namespace analysis
