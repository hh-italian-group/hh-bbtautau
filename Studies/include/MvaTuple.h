/*! Definition of MvaVariablesStudy class, the main class for Mva studies.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"
#include "hh-bbtautau/Analysis/include/MvaVariables.h"

namespace  analysis {

#define MVA_DATA() \
    VAR(UInt_t, NTrees) \
    VAR(double, shrinkage) \
    VAR(double, BaggedSampleFraction) \
    VAR(double, MaxDepth) \
    VAR(double, MinNodeSize) \
    VAR(double, ROCIntegral) \
    /**/ \
    VAR(std::vector<double>, roc_value) \
    VAR(std::vector<int>, roc_mass) \
    /**/ \
    VAR(std::vector<double>, KS_value) \
    VAR(std::vector<int>, KS_type) \
    VAR(std::vector<int>, KS_mass) \
    /**/ \
    VAR(std::vector<double>, position) \
    VAR(std::vector<double>, importance) \
    VAR(std::vector<std::string>, var_name) \
    /**/ \
    VAR(double, cut) \
    VAR(double, significance) \
    VAR(double, err_significance) \
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

inline std::map<int, double> GetRocIntegralMap(const MvaResults& results)
{
    std::map<int, double> rocs;
    if(results.roc_mass.size() != results.roc_value.size())
        throw exception("Incompatible roc info in mva tuple.");
    for(size_t n = 0; n < results.roc_mass.size(); ++n)
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

struct VarRank { size_t position; double importance; };

inline std::map<std::string, VarRank> GetRankingMap(const MvaResults& results)
{
    std::map<std::string, VarRank> ranks;
    const size_t N = results.var_name.size();
    if(results.position.size() != N || results.importance.size() != N)
        throw exception("Incompatible ranking info in mva tuple.");
    for(size_t n = 0; n < N; ++n)
        ranks[results.var_name[n]] = VarRank{results.position[n], results.importance[n]};
    return ranks;
}

}
}
