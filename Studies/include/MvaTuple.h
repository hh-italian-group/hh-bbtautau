/*! Definition of MvaVariablesStudy class, the main class for Mva studies.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"
#include "hh-bbtautau/Studies/include/MvaMethods.h"

namespace  analysis {

#define MVA_DATA() \
    VAR(std::vector<std::string>, param_names) \
    VAR(std::vector<double>, param_values) \
    VAR(std::vector<size_t>, param_positions) \
    /**/ \
    VAR(std::vector<double>, roc_testing_value) \
    VAR(std::vector<int>, roc_testing_mass) \
    VAR(std::vector<int>, roc_testing_type) \
    VAR(std::vector<int>, roc_testing_spin) \
    VAR(std::vector<std::string>, roc_testing_channel) \
    VAR(std::vector<double>, err_roc_testing) \
    /**/ \
    VAR(std::vector<double>, roc_training_value) \
    VAR(std::vector<int>, roc_training_mass) \
    VAR(std::vector<int>, roc_training_type) \
    VAR(std::vector<int>, roc_training_spin) \
    VAR(std::vector<std::string>, roc_training_channel) \
    VAR(std::vector<double>, err_roc_training) \
    /**/ \
    VAR(std::vector<double>, KS_value) \
    VAR(std::vector<int>, KS_type) \
    VAR(std::vector<int>, KS_mass) \
    VAR(std::vector<int>, KS_spin) \
    VAR(std::vector<std::string>, KS_channel) \
    /**/ \
    VAR(std::vector<double>, chi_value) \
    VAR(std::vector<int>, chi_type) \
    VAR(std::vector<int>, chi_mass) \
    VAR(std::vector<int>, chi_spin) \
    VAR(std::vector<std::string>, chi_channel) \
    /**/ \
    VAR(std::vector<size_t>, position) \
    VAR(std::vector<double>, importance) \
    VAR(std::vector<std::string>, var_name) \
    /**/ \
    VAR(std::vector<double>, optimal_cut) \
    VAR(std::vector<double>, significance) \
    VAR(std::vector<double>, significance_err) \
    VAR(std::vector<int>, significance_mass) \
    VAR(std::vector<int>, significance_type) \
    VAR(std::vector<int>, significance_spin) \
    VAR(std::vector<std::string>, significance_channel) \
    /**/ \
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

GridPoint GetGridPoint(const MvaResults& results);

std::map<ChannelSampleIdSpin, PhysicalValue> GetRocIntegralMap(const std::string& name,
                                                               const std::vector<double>& vec_value,
                                                               const std::vector<double>& vec_err,
                                                               const std::vector<std::string>& vec_channel,
                                                               const std::vector<int>& vec_mass,
                                                               const std::vector<int>& vec_type,
                                                               const std::vector<int>& vec_spin);

std::map<ChannelSampleIdSpin, PhysicalValue> GetRocTrainingIntegralMap(const MvaResults& results);
std::map<ChannelSampleIdSpin, PhysicalValue> GetRocTestingIntegralMap(const MvaResults& results);
std::map<ChannelSampleIdSpin, double> GetTestResultsMap(const std::string& name, const std::vector<int>& vec_type,
                                                        const std::vector<std::string>& vec_channel,
                                                        const std::vector<double>& vec_value,
                                                        const std::vector<int>& vec_mass,
                                                        const std::vector<int>& vec_spin);

std::map<ChannelSampleIdSpin, double>  GetKSResultsMap(const MvaResults& results);
std::map<ChannelSampleIdSpin, double>  GetChiResultsMap(const MvaResults& results);

struct VarRank {
    size_t position;
    double importance;
    bool operator<(const VarRank& other) const;
};

using VarRankMap = std::map<std::string, VarRank>;

VarRankMap GetRankingMap(const MvaResults& results);

using OptimalSignificanceMap = std::map<ChannelSampleIdSpin, OptimalSignificance>;

OptimalSignificanceMap GetOptimalSignificanceMap(const MvaResults& results);

}
}
