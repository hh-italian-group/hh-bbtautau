/*! Class than provides nonresonant EFT model for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "h-tautau/Analysis/include/EventInfo.h"
#include "hh-bbtautau/McCorrections/include/EventWeights_HH.h"
#include "SampleDescriptor.h"
#include "AnaTuple.h"

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
        PointDesc(const Point& _point, double _total_shape_weight) :
            point(_point), total_shape_weight(_total_shape_weight) {}
    };

    struct ParamPositionDesc {
        using NameMap = std::map<std::string, size_t>;
        using ValueVec = std::vector<double>;
        using OptPos = boost::optional<size_t>;
        OptPos kl, kt, c2, cg, c2g;

        ParamPositionDesc(const NameMap& names)
        {
            SetParamPosition(names, "kl", kl);
            SetParamPosition(names, "kt", kt);
            SetParamPosition(names, "c2", c2);
            SetParamPosition(names, "cg", cg);
            SetParamPosition(names, "c2g", c2g);
        }

        Point CreatePoint(const ValueVec& param_values) const
        {
            Point point;
            SetValue(param_values, "kl", kl, point.kl);
            SetValue(param_values, "kt", kt, point.kt);
            SetValue(param_values, "c2", c2, point.c2);
            SetValue(param_values, "cg", cg, point.cg);
            SetValue(param_values, "c2g", c2g, point.c2g);
            return point;
        }

        static void SetParamPosition(const NameMap& names, const std::string& name, OptPos& pos)
        {
            auto iter = names.find(name);
            if(iter != names.end())
                pos = iter->second;
        }

        static void SetValue(const ValueVec& values, const std::string& name, const OptPos& pos, double& value)
        {
            if(pos) {
                if(values.size() <= *pos)
                    throw exception("Value not found for EFT parameter %1%.") % name;
                value = values.at(*pos);
            }
        }
    };

public:

    NonResModel(Period period, const SampleDescriptor& sample, std::shared_ptr<TFile> file)
        : weighting_mode({WeightType::PileUp, WeightType::BSM_to_SM}),
          weights(period, DiscriminatorWP::Medium, false, weighting_mode),
          eft_weights(weights.GetProviderT<NonResHH_EFT::WeightProvider>(WeightType::BSM_to_SM))
    {
        const ParamPositionDesc param_positions(sample.GetModelParameterNames());
        eft_weights->AddFile(*file);
        eft_weights->CreatePdfs();
        for(const auto& sample_wp : sample.working_points) {
            const Point point = param_positions.CreatePoint(sample_wp.param_values);
            eft_weights->SetTargetPoint(point);
            const auto summary = weights.GetSummaryWithWeights(file, weighting_mode);
            points[sample_wp.full_name] = PointDesc(point, summary.totalShapeWeight);
        }
    }

    void ProcessEvent(const EventAnalyzerDataId& anaDataId, EventInfoBase& event, double weight, double shape_weight,
                      bbtautau::AnaTupleWriter::DataIdMap& dataIds)
    {
        for(const auto& wp : points) {
            const auto final_id = anaDataId.Set(wp.first);
            eft_weights->SetTargetPoint(wp.second.point);
            const double eft_weight = eft_weights->Get(*event);
            const double final_weight = weight * (shape_weight / wp.second.total_shape_weight)
                    * (eft_weight / event->weight_bsm_to_sm);
            dataIds[final_id] = std::make_tuple(final_weight, event.GetMvaScore());
        }
    }

private:
    WeightingMode weighting_mode;
    mc_corrections::EventWeights_HH weights;
    std::shared_ptr<NonResHH_EFT::WeightProvider> eft_weights;
    std::map<std::string, PointDesc> points;
};

} // namespace analysis
