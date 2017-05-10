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
    using Range_weight_map = std::map<Range<int>, double>;
    using DoubleRange_map = std::map<Range<int>, Range_weight_map>;

    DY_weight(const std::string& dy_weight_file_name) :
    {
        std::vector<analysis::sample_merging::DYBinDescriptor> dy_descriptors = analysis::sample_merging::DYBinDescriptor::LoadConfig(dy_weight_file_name);
        for (unsigned n = 0; n < dy_descriptors.size(); ++n){
            const analysis::sample_merging::DYBinDescriptor dybin_descriptor = dy_descriptors.at(n);
            dy_weight_map[dybin_descriptor.n_jet][dybin_descriptor.n_bjet][dybin_descriptor.n_ht] = dybin_descriptor.weight.GetValue();
        }
    }

    double Get(const Event& event)
    {
        for (auto iter : dy_weight_map){
            Range<int> n_jet = iter.first;
            if (!(n_jet.Contains(event.lhe_n_partons))) continue;
            DoubleRange_map doubleRange_map = iter.second;
            for (auto iter_1 : doubleRange_map){
                Range<int> n_bjet = iter_1.first;
                if (!(n_bjet.Contains(event.lhe_n_b_partons))) continue;
                Range_weight_map range_weight_map = iter_1.second;
                for (auto iter_2 : range_weight_map){
                    Range<int> n_ht = iter_2.first;
                    if (n_ht.Contains(event.lhe_HT))
                        return iter_2.second;
                }
            }
        }
        throw analysis::exception("drell-yan merge weight not found.");
    }

private:
    std::map<Range<int>, DoubleRange_map> dy_weight_map;


};

} // namespace mc_corrections
} // namespace analysis
