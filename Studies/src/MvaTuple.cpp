/*! Definition of MvaVariablesStudy class, the main class for Mva studies.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Studies/include/MvaTuple.h"

namespace  analysis {
namespace mva_study {

GridPoint GetGridPoint(const MvaResults& results)
{
    GridPoint params;
    const size_t N = results.param_names.size();
    if(results.param_positions.size() != N || results.param_values.size() != N)
        throw exception("Incompatible grid point info in mva tuple.");
    for(size_t n = 0; n < N; ++n)
        params[results.param_names[n]] = GridParam{results.param_positions[n], results.param_values[n]};
    return params;
}

std::map<ChannelSampleIdSpin, PhysicalValue> GetRocIntegralMap(const std::string& name,
                                                               const std::vector<double>& vec_value,
                                                               const std::vector<double>& vec_err,
                                                               const std::vector<std::string>& vec_channel,
                                                               const std::vector<int>& vec_mass,
                                                               const std::vector<int>& vec_type,
                                                               const std::vector<int>& vec_spin)
{
    std::map<ChannelSampleIdSpin, PhysicalValue> rocs;
    const size_t N = vec_mass.size();
    if(vec_value.size() != N || vec_spin.size() != N || vec_channel.size()!=N || vec_type.size()!=N || vec_err.size()!=N)
        throw exception("Incompatible "+name+" roc info in mva tuple.");
    for(size_t n = 0; n < N; ++n){
        const SampleId sample_id(static_cast<SampleType>(vec_type[n]), vec_mass[n]);
        const ChannelSampleIdSpin id{vec_channel[n], sample_id, vec_spin[n]};
        rocs[id] = PhysicalValue(vec_value[n], vec_err[n]);
    }
    return rocs;
}

std::map<ChannelSampleIdSpin, PhysicalValue> GetRocTrainingIntegralMap(const MvaResults& results)
{
    std::map<ChannelSampleIdSpin, PhysicalValue> rocs;
    return rocs = GetRocIntegralMap("training", results.roc_training_value, results.err_roc_training, results.roc_training_channel,
                                    results.roc_training_mass, results.roc_training_type, results.roc_training_spin);
}

std::map<ChannelSampleIdSpin, PhysicalValue> GetRocTestingIntegralMap(const MvaResults& results)
{
    std::map<ChannelSampleIdSpin, PhysicalValue>  rocs;
    return rocs = GetRocIntegralMap("testing", results.roc_testing_value, results.err_roc_testing, results.roc_testing_channel,
                                    results.roc_testing_mass,results.roc_testing_type, results.roc_testing_spin);
}

std::map<ChannelSampleIdSpin, double> GetTestResultsMap(const std::string& name, const std::vector<int>& vec_type,
                                                        const std::vector<std::string>& vec_channel,
                                                        const std::vector<double>& vec_value,
                                                        const std::vector<int>& vec_mass,
                                                        const std::vector<int>& vec_spin)
{
    std::map<ChannelSampleIdSpin, double> test;
    const size_t N = vec_mass.size();
    if(vec_type.size() != N || vec_value.size() != N || vec_spin.size() != N || vec_channel.size() != N)
        throw exception("Incompatible "+name+" info in mva tuple.");
    for(size_t n = 0; n < N; ++n) {
        const SampleId sample_id(static_cast<SampleType>(vec_type[n]), vec_mass[n]);
        const ChannelSampleIdSpin id{vec_channel[n], sample_id, vec_spin[n]};
        test[id] = vec_value[n];
    }
    return test;
}

std::map<ChannelSampleIdSpin, double> GetKSResultsMap(const MvaResults& results)
{
    std::map<ChannelSampleIdSpin, double>  ks;
    return ks = GetTestResultsMap("KS", results.KS_type, results.KS_channel, results.KS_value, results.KS_mass, results.KS_spin);
}

std::map<ChannelSampleIdSpin, double> GetChiResultsMap(const MvaResults& results)
{
    std::map<ChannelSampleIdSpin, double> chi;
    return chi = GetTestResultsMap("chi", results.chi_type, results.chi_channel, results.chi_value, results.chi_mass, results.chi_spin);
}

bool VarRank::operator<(const VarRank& other) const { return position < other.position; }

VarRankMap GetRankingMap(const MvaResults& results)
{
    VarRankMap ranks;
    const size_t N = results.var_name.size();
    if(results.position.size() != N || results.importance.size() != N)
        throw exception("Incompatible ranking info in mva tuple.");
    for(size_t n = 0; n < N; ++n)
        ranks[results.var_name[n]] = VarRank{results.position[n], results.importance[n]};
    return ranks;
}

OptimalSignificanceMap GetOptimalSignificanceMap(const MvaResults& results)
{
    OptimalSignificanceMap significance;
    const size_t N = results.significance_mass.size();
    if(results.significance.size() != N || results.significance_err.size() != N || results.optimal_cut.size() != N
            || results.significance_type.size() != N || results.significance_spin.size() != N || results.significance_channel.size() != N)
        throw exception("Incompatible significance info in mva tuple.");
    for(size_t n = 0; n < N; ++n)
    {
        PhysicalValue sign(results.significance[n], results.significance_err[n]);
        const SampleId sample_id(static_cast<SampleType>(results.significance_type[n]), results.significance_mass[n]);
        ChannelSampleIdSpin id(results.significance_channel[n], sample_id, results.significance_spin[n]);
        significance[id] = OptimalSignificance{results.optimal_cut[n], sign};
    }
    return significance;
}

}
}
