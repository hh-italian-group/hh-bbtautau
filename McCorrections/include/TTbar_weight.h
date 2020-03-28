/*! The ttbar weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "h-tautau/McCorrections/include/WeightProvider.h"
#include "TTFileDescriptor.h"

namespace analysis {
namespace mc_corrections {

class TTbar_weight : public IWeightProvider {
public:
    TTbar_weight(const std::string& ttbar_weight_file_name);
    virtual double Get(EventInfo& eventInfo) const override;
    virtual double Get(const ntuple::ExpressEvent& event) const override;

private:
    double Get(int genEventType) const;

private:
    std::map<int, double> genEventType_weight_map;
};

} // namespace mc_corrections
} // namespace analysis
