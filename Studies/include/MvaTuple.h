/*! Definition of MvaVariablesStudy class, the main class for Mva studies.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"
#include "hh-bbtautau/Analysis/include/MvaVariables.h"

namespace  analysis {

#define MVA_DATA() \
    VAR(std::vector<std::string>, param_names) \
    VAR(std::vector<double>, param_values) \
    VAR(std::vector<size_t>, param_positions) \
    /**/ \
    VAR(double, ROCIntegral) \
    /**/ \
    VAR(std::vector<double>, roc_value) \
    VAR(std::vector<int>, roc_mass) \
    /**/ \
    VAR(std::vector<double>, KS_value) \
    VAR(std::vector<int>, KS_type) \
    VAR(std::vector<int>, KS_mass) \
    /**/ \
    VAR(std::vector<double>, chi_value) \
    VAR(std::vector<int>, chi_type) \
    VAR(std::vector<int>, chi_mass) \
    /**/ \
    VAR(std::vector<size_t>, position) \
    VAR(std::vector<double>, importance) \
    VAR(std::vector<std::string>, var_name) \
    /**/ \
    VAR(std::vector<double>, optimal_cut) \
    VAR(std::vector<double>, significance) \
    VAR(std::vector<double>, significance_err) \
    VAR(std::vector<int>, significance_mass) \
    VAR(std::string, name) \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(mva_study, MvaResults, MvaTuple, MVA_DATA, "mva_result")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(mva_study, MvaTuple, MVA_DATA)
#undef VAR
#undef MVA_DATA


namespace mva_study{

struct GridParam { size_t position; double value; };
using GridPoint = std::map<std::string, GridParam>;

inline GridPoint GetGridPoint(const MvaResults& results)
{
    GridPoint params;
    const size_t N = results.param_names.size();
    if(results.param_positions.size() != N || results.param_values.size() != N)
        throw exception("Incompatible grid point info in mva tuple.");
    for(size_t n = 0; n < N; ++n)
        params[results.param_names[n]] = GridParam{results.param_positions[n], results.param_values[n]};
    return params;
}

inline std::map<int, double> GetRocIntegralMap(const MvaResults& results)
{
    std::map<int, double> rocs;
    const size_t N = results.roc_mass.size();
    if(results.roc_value.size() != N)
        throw exception("Incompatible roc info in mva tuple.");
    for(size_t n = 0; n < N; ++n)
        rocs[results.roc_mass[n]] = results.roc_value[n];
    return rocs;
}

inline std::map<SampleId, double> GetKSResultsMap(const MvaResults& results)
{
    std::map<SampleId, double> ks;
    const size_t N = results.KS_mass.size();
    if(results.KS_type.size() != N || results.KS_value.size() != N)
        throw exception("Incompatible KS info in mva tuple.");
    for(size_t n = 0; n < N; ++n) {
        const SampleType type = static_cast<SampleType>(results.KS_type[n]);
        const SampleId id(type, results.KS_mass[n]);
        ks[id] = results.KS_value[n];
    }
    return ks;
}

inline std::map<SampleId, double> GetChiResultsMap(const MvaResults& results)
{
    std::map<SampleId, double> chi;
    const size_t N = results.chi_mass.size();
    if(results.chi_type.size() != N || results.chi_value.size() != N)
        throw exception("Incompatible Chi info in mva tuple.");
    for(size_t n = 0; n < N; ++n) {
        const SampleType type = static_cast<SampleType>(results.chi_type[n]);
        const SampleId id(type, results.chi_mass[n]);
        chi[id] = results.chi_value[n];
    }
    return chi;
}


struct VarRank {
    size_t position;
    double importance;
    bool operator<(const VarRank& other) const { return position < other.position; }
};
using VarRankMap = std::map<std::string, VarRank>;

inline VarRankMap GetRankingMap(const MvaResults& results)
{
    VarRankMap ranks;
    const size_t N = results.var_name.size();
    if(results.position.size() != N || results.importance.size() != N)
        throw exception("Incompatible ranking info in mva tuple.");
    for(size_t n = 0; n < N; ++n)
        ranks[results.var_name[n]] = VarRank{results.position[n], results.importance[n]};
    return ranks;
}

struct OptimalSignificance {
    double cut;
    PhysicalValue significance;
};

using OptimalSignificanceMap = std::map<int, OptimalSignificance>;

inline OptimalSignificanceMap GetOptimalSignificanceMap(const MvaResults& results)
{
    OptimalSignificanceMap significance;
    const size_t N = results.significance_mass.size();
    if(results.significance.size() != N || results.significance_err.size() != N || results.optimal_cut.size() != N)
        throw exception("Incompatible significance info in mva tuple.");
    for(size_t n = 0; n < N; ++n)
    {
        PhysicalValue sign(results.significance[n], results.significance_err[n]);
        significance[results.significance_mass[n]] = OptimalSignificance{results.optimal_cut[n], sign};
    }
    return significance;
}
}
}
