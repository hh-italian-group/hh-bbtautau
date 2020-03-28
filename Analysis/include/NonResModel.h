/*! Class than provides nonresonant EFT model for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "h-tautau/Analysis/include/EventInfo.h"
#include "hh-bbtautau/McCorrections/include/EventWeights_HH.h"
#include "hh-bbtautau/McCorrections/include/HH_nonResonant_weight.h"

#include "AnaTuple.h"
#include "SampleDescriptor.h"

namespace analysis {

class NonResModel {
private:
    using WeightingMode = mc_corrections::WeightingMode;
    using WeightType = mc_corrections::WeightType;
    using Point = NonResHH_EFT::Point;

    struct PointDesc {
        Point point;
        double total_shape_weight;

        PointDesc() {}
        PointDesc(const Point& _point, double _total_shape_weight);
    };

    struct ParamPositionDesc {
        using NameMap = std::map<std::string, size_t>;
        using ValueVec = std::vector<double>;
        using OptPos = boost::optional<size_t>;
        OptPos kl, kt, c2, cg, c2g;

        ParamPositionDesc(const NameMap& names);
        Point CreatePoint(const ValueVec& param_values) const;

        static void SetParamPosition(const NameMap& names, const std::string& name, OptPos& pos);
        static void SetValue(const ValueVec& values, const std::string& name, const OptPos& pos, double& value);
    };

public:
    NonResModel(Period period, const SampleDescriptor& sample, std::shared_ptr<TFile> file);
    void ProcessEvent(const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight, double shape_weight,
                      bbtautau::AnaTupleWriter::DataIdMap& dataIds);

private:
    WeightingMode weighting_mode;
    mc_corrections::EventWeights_HH weights;
    std::shared_ptr<NonResHH_EFT::WeightProvider> eft_weights;
    std::map<std::string, PointDesc> points;
};

} // namespace analysis
