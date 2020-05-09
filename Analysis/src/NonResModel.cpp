/*! Class than provides nonresonant EFT model for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/NonResModel.h"

#include "h-tautau/McCorrections/include/PileUpWeight.h"

namespace analysis {

NonResModel::PointDesc::PointDesc(const Point& _point, double _total_shape_weight) :
    point(_point), total_shape_weight(_total_shape_weight) {}

NonResModel::ParamPositionDesc::ParamPositionDesc(const NameMap& names)
{
    SetParamPosition(names, "kl", kl);
    SetParamPosition(names, "kt", kt);
    SetParamPosition(names, "c2", c2);
    SetParamPosition(names, "cg", cg);
    SetParamPosition(names, "c2g", c2g);
}

NonResModel::Point NonResModel::ParamPositionDesc::CreatePoint(const ValueVec& param_values) const
{
    Point point;
    SetValue(param_values, "kl", kl, point.kl);
    SetValue(param_values, "kt", kt, point.kt);
    SetValue(param_values, "c2", c2, point.c2);
    SetValue(param_values, "cg", cg, point.cg);
    SetValue(param_values, "c2g", c2g, point.c2g);
    return point;
}

void NonResModel::ParamPositionDesc::SetParamPosition(const NameMap& names, const std::string& name, OptPos& pos)
{
    auto iter = names.find(name);
    if(iter != names.end())
        pos = iter->second;
}

void NonResModel::ParamPositionDesc::SetValue(const ValueVec& values, const std::string& name,
                                              const OptPos& pos, double& value)
{
    if(pos) {
        if(values.size() <= *pos)
            throw exception("Value not found for EFT parameter %1%.") % name;
        value = Parse<double>(values.at(*pos));
    }
}

NonResModel::NonResModel(Period period, const SampleDescriptor& sample, std::shared_ptr<TFile> file) :
        weighting_mode(WeightType::PileUp, WeightType::BSM_to_SM, WeightType::GenEventWeight),
        weights(period, BTagger(period, BTaggerKind::DeepFlavour), weighting_mode),
        eft_weights(weights.GetProviderT<NonResHH_EFT::WeightProvider>(WeightType::BSM_to_SM))
{
    const ParamPositionDesc param_positions(sample.GetModelParameterNames());
    eft_weights->AddFile(*file);
    eft_weights->CreatePdfs();

    if(weighting_mode.count(WeightType::PileUp) && period == Period::Run2017){
        auto pile_up_weight = weights.GetProviderT<mc_corrections::PileUpWeightEx>(mc_corrections::WeightType::PileUp);
        pile_up_weight->SetActiveDataset(sample.reference_pu_sample);
    }

    for(const auto& sample_wp : sample.working_points) {
        const Point point = param_positions.CreatePoint(sample_wp.param_values);
        eft_weights->SetTargetPoint(point);
        const auto summary = weights.GetSummaryWithWeights(file, weighting_mode);
        points[sample_wp.full_name] = PointDesc(point, summary.totalShapeWeight);
    }
}

void NonResModel::ProcessEvent(const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight,
                               double shape_weight, bbtautau::AnaTupleWriter::DataIdMap& dataIds)
{
    for(const auto& wp : points) {
        const auto final_id = anaDataId.Set(wp.first);
        eft_weights->SetTargetPoint(wp.second.point);
        const double eft_weight = eft_weights->Get(event);
        const double final_weight = weight * (shape_weight / wp.second.total_shape_weight)
                * (eft_weight / event->weight_bsm_to_sm);
        dataIds[final_id] = std::make_tuple(final_weight, event.GetMvaScore());
    }
}

} // namespace analysis
