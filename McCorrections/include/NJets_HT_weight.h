/*! The DrellYan weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "h-tautau/McCorrections/include/WeightProvider.h"
#include "NJets_HT_BinFileDescriptor.h"

namespace analysis {
namespace mc_corrections {

class NJets_HT_weight : public IWeightProvider {
public:
    using Event = ntuple::Event;
    using Range_weight_map = std::map<size_t, double>;
    using DoubleRange_map = std::map<size_t, Range_weight_map>;

    NJets_HT_weight(const std::string& _name, const std::string& weight_file_name);

    virtual double Get(EventInfo& eventInfo) const override;
    virtual double Get(const ntuple::ExpressEvent& event) const override;

    double GetWeight(size_t n_partons, size_t n_b_partons, size_t ht_bin) const;

private:
    std::string name;
    std::map<size_t, DoubleRange_map> weight_map;
};

} // namespace mc_corrections
} // namespace analysis
