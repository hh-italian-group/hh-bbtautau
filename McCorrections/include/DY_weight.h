/*! The DrellYan weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "../../Instruments/include/DYFileDescriptor.h"


namespace analysis {
namespace mc_corrections {

class DY_weight {
public:
    using Event = ntuple::Event;
    using Range_weight_map = std::map<int, double>;
    using DoubleRange_map = std::map<int, Range_weight_map>;

    DY_weight(const std::string& dy_weight_file_name) :
        dy_descriptors(analysis::sample_merging::DYBinDescriptor::LoadConfig(dy_weight_file_name))
    {
        //std::vector<analysis::sample_merging::DYBinDescriptor> dy_descriptors = analysis::sample_merging::DYBinDescriptor::LoadConfig(dy_weight_file_name);

        for (unsigned n = 0; n < dy_descriptors.size(); ++n){
            const analysis::sample_merging::DYBinDescriptor dybin_descriptor = dy_descriptors.at(n);
            const double weight = dybin_descriptor.weight.GetValue()/dybin_descriptor.inclusive_integral;
            for(int n_jet = dybin_descriptor.n_jet.min(); n_jet <= dybin_descriptor.n_jet.max(); ++n_jet) {
                if(!dy_weight_map.count(n_jet)) {
                    dy_weight_map[n_jet] = DoubleRange_map();
                }
                DoubleRange_map& njet_map = dy_weight_map.at(n_jet);
                for(int n_bjet = dybin_descriptor.n_bjet.min(); n_bjet <= dybin_descriptor.n_bjet.max(); ++n_bjet) {
                    if(!njet_map.count(n_bjet))
                        njet_map[n_bjet] = Range_weight_map();
                    Range_weight_map& nbjet_map = njet_map.at(n_bjet);
                    for(int ht = dybin_descriptor.n_ht.min(); ht <= dybin_descriptor.n_ht.max(); ++ht) {
                        if(nbjet_map.count(ht))
                            throw exception("Repeated bin");
                        nbjet_map[ht] = weight;
                    }
                }
            }
        }
    }

    template<typename Event>
    double Get(const Event& event)
    {
        auto njet_iter = dy_weight_map.find(event.lhe_n_partons);
        if(njet_iter != dy_weight_map.end()) {
            const auto& nbjet_map = njet_iter->second;
            auto nbjet_iter = nbjet_map.find(event.lhe_n_b_partons);
            if(nbjet_iter != nbjet_map.end()){
                const auto& nht_map = nbjet_iter->second;
                auto nht_iter = nht_map.find(event.lhe_HT);
                if(nht_iter != nht_map.end())
                    return nht_iter->second;
            }
        }
        throw exception("weight not found.");
    }

private:
    std::map<int, DoubleRange_map> dy_weight_map;
    std::vector<analysis::sample_merging::DYBinDescriptor> dy_descriptors;


};

} // namespace mc_corrections
} // namespace analysis
