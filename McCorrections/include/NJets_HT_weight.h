/*! The DrellYan weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "hh-bbtautau/Instruments/include/NJets_HT_BinFileDescriptor.h"


namespace analysis {
namespace mc_corrections {

class NJets_HT_weight {
public:
    using Event = ntuple::Event;
    using Range_weight_map = std::map<int, double>;
    using DoubleRange_map = std::map<int, Range_weight_map>;

    NJets_HT_weight(const std::string& weight_file_name)
    {
        std::vector<analysis::sample_merging::NJets_HT_BinFileDescriptor> descriptors =
                analysis::sample_merging::NJets_HT_BinFileDescriptor::LoadConfig(weight_file_name);

        for (unsigned n = 0; n < descriptors.size(); ++n){
            const analysis::sample_merging::NJets_HT_BinFileDescriptor bin_descriptor = descriptors.at(n);
            const double weight = bin_descriptor.weight.GetValue()/bin_descriptor.inclusive_integral;
            for(int n_jet = bin_descriptor.n_jet.min(); n_jet <= bin_descriptor.n_jet.max(); ++n_jet) {
                DoubleRange_map& njet_map = weight_map[n_jet];
                for(int n_bjet = bin_descriptor.n_bjet.min(); n_bjet <= bin_descriptor.n_bjet.max(); ++n_bjet) {
                    Range_weight_map& nbjet_map = njet_map[n_bjet];
                    for(int ht = bin_descriptor.n_ht.min(); ht <= bin_descriptor.n_ht.max(); ++ht) {
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
        auto njet_iter = weight_map.find(n_partons);
        if(njet_iter != weight_map.end()) {
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
    std::map<int, DoubleRange_map> weight_map;
};

} // namespace mc_corrections
} // namespace analysis
