/*! The DrellYan weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "hh-bbtautau/Instruments/include/DYFileDescriptor.h"


namespace analysis {
namespace mc_corrections {

class DY_weight {
public:
    using Event = ntuple::Event;
    using Range_weight_map = std::map<int, double>;
    using DoubleRange_map = std::map<int, Range_weight_map>;

    DY_weight(const std::string& dy_weight_file_name)
    {
        std::vector<analysis::sample_merging::DYBinDescriptor> dy_descriptors =
                analysis::sample_merging::DYBinDescriptor::LoadConfig(dy_weight_file_name);

        for (unsigned n = 0; n < dy_descriptors.size(); ++n){
            const analysis::sample_merging::DYBinDescriptor dybin_descriptor = dy_descriptors.at(n);
            const double weight = dybin_descriptor.weight.GetValue()/dybin_descriptor.inclusive_integral;
            for(int n_jet = dybin_descriptor.n_jet.min(); n_jet <= dybin_descriptor.n_jet.max(); ++n_jet) {
                DoubleRange_map& njet_map = dy_weight_map[n_jet];
                for(int n_bjet = dybin_descriptor.n_bjet.min(); n_bjet <= dybin_descriptor.n_bjet.max(); ++n_bjet) {
                    Range_weight_map& nbjet_map = njet_map[n_bjet];
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
    double Get(const Event& event) const
    {
        return GetWeight(event.lhe_n_partons,event.lhe_n_b_partons,event.lhe_HT);
    }

    double GetWeight(Int_t n_partons, Int_t n_b_partons, Int_t ht) const
    {
        auto njet_iter = dy_weight_map.find(n_partons);
        if(njet_iter != dy_weight_map.end()) {
            const auto& nbjet_map = njet_iter->second;
            auto nbjet_iter = nbjet_map.find(n_b_partons);
            if(nbjet_iter != nbjet_map.end()){
                const auto& nht_map = nbjet_iter->second;
                auto nht_iter = nht_map.find(ht);
                if(nht_iter != nht_map.end())
                    return nht_iter->second;
            }
        }
        throw exception("weight not found.");
    }

private:
    std::map<int, DoubleRange_map> dy_weight_map;
};

} // namespace mc_corrections
} // namespace analysis
