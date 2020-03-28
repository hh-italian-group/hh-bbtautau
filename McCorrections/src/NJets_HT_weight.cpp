/*! The DrellYan weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/McCorrections/include/NJets_HT_weight.h"

namespace analysis {
namespace mc_corrections {

NJets_HT_weight::NJets_HT_weight(const std::string& _name, const std::string& weight_file_name) :
    name(_name)
{
    std::vector<analysis::sample_merging::NJets_HT_BinFileDescriptor> descriptors =
            analysis::sample_merging::NJets_HT_BinFileDescriptor::LoadConfig(weight_file_name);

    for (unsigned n = 0; n < descriptors.size(); ++n){
        const analysis::sample_merging::NJets_HT_BinFileDescriptor bin_descriptor = descriptors.at(n);
        const double weight = bin_descriptor.weight.GetValue()/bin_descriptor.inclusive_integral;
        for(size_t n_jet = bin_descriptor.n_jet.min(); n_jet <= bin_descriptor.n_jet.max(); ++n_jet) {
            DoubleRange_map& njet_map = weight_map[n_jet];
            for(size_t n_bjet = bin_descriptor.n_bjet.min(); n_bjet <= bin_descriptor.n_bjet.max(); ++n_bjet) {
                Range_weight_map& nbjet_map = njet_map[n_bjet];
                for(size_t ht = bin_descriptor.n_ht.min(); ht <= bin_descriptor.n_ht.max(); ++ht) {
                    if(nbjet_map.count(ht))
                        throw exception("Repeated bin");
                    nbjet_map[ht] = weight;
                }
            }
        }
    }
}

double NJets_HT_weight::Get(EventInfo& eventInfo) const
{
    static constexpr size_t ht_bin_size = 10;
    return GetWeight(eventInfo->lhe_n_partons, eventInfo->lhe_n_b_partons, static_cast<size_t>(eventInfo->lhe_HT / ht_bin_size));
}

double NJets_HT_weight::Get(const ntuple::ExpressEvent& event) const
{
    return GetWeight(event.lhe_n_partons, event.lhe_n_b_partons, event.lhe_ht10_bin);
}

double NJets_HT_weight::GetWeight(size_t n_partons, size_t n_b_partons, size_t ht_bin) const
{
    auto njet_iter = weight_map.find(n_partons);
    if(njet_iter != weight_map.end()) {
        const auto& nbjet_map = njet_iter->second;
        auto nbjet_iter = nbjet_map.find(n_b_partons);
        if(nbjet_iter != nbjet_map.end()){
            const auto& nht_map = nbjet_iter->second;
            auto nht_iter = nht_map.find(ht_bin);
            if(nht_iter != nht_map.end())
                return nht_iter->second;
        }
    }
    throw exception("%1% weight not found for n_partons = %2%, n_b_partons = %3%, ht_bin = %4%.")
            % name % n_partons % n_b_partons % ht_bin;
}

} // namespace mc_corrections
} // namespace analysis
