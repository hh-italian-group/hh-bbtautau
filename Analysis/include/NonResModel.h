/*! Class than provides nonresonant EFT model for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "h-tautau/Analysis/include/EventInfo.h"
#include "hh-bbtautau/McCorrections/include/EventWeights_HH.h"
#include "hh-bbtautau/McCorrections/include/HH_nonResonant_weight.h"
#include "hh-bbtautau/Instruments/include/SkimmerConfig.h"

#include "AnaTuple.h"
#include "SampleDescriptor.h"

namespace analysis {

class NonResModel {
private:
    using WeightingMode = mc_corrections::WeightingMode;
    using WeightType = mc_corrections::WeightType;
    using Point = NonResHH_EFT::Point;

    struct ParamPositionDesc {
        using NameMap = std::map<std::string, size_t>;
        using ValueVec = std::vector<std::string>;
        using OptPos = boost::optional<size_t>;
        OptPos kl, kt, c2, cg, c2g;

        ParamPositionDesc(const NameMap& names);
        Point CreatePoint(const ValueVec& param_values) const;

        static void SetParamPosition(const NameMap& names, const std::string& name, OptPos& pos);
        static void SetValue(const ValueVec& values, const std::string& name, const OptPos& pos, double& value);
    };

public:
    NonResModel(Period period, const SampleDescriptor& sample, std::shared_ptr<TFile> file,
                tuple_skimmer::CrossSectionProvider& xs_provider);
    void ProcessEvent(const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight, double shape_weight,
                      bbtautau::AnaTupleWriter::DataIdMap& dataIds, double cross_section);

private:
    WeightingMode weighting_mode;
    mc_corrections::EventWeights_HH weights;
    std::shared_ptr<NonResHH_EFT::WeightProvider> eft_weights;
    std::vector<std::string> point_names;
    std::vector<Point> points;
    std::vector<double> total_shape_weights;
    std::vector<double> point_xs;
    bool points_are_orthogonal;
};

} // namespace analysis
