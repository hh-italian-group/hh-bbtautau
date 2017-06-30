/*! Definition of MvaVariablesStudy class, the main class for Mva studies.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"

namespace  analysis {
namespace mva_study{

#define MVA_DATA() \
    VAR(UInt_t, NTrees) \
    VAR(double, shrinkage) \
    VAR(double, BaggedSampleFraction) \
    VAR(double, MaxDepth) \
    VAR(double, MinNodeSize) \
    VAR(double, ROCIntegral) \
    VAR(std::vector<double>, KS_value) \
    VAR(std::vector<int>, KS_type) \
    VAR(std::vector<int>, KS_mass) \
    VAR(std::vector<double>, position) \
    VAR(std::vector<double>, importance) \
    VAR(std::vector<std::string>, var_name) \
    VAR(double, cut) \
    VAR(double, significance) \

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(ntuple, MvaResults, MvaTuple, MVA_DATA, "mva_result")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(ntuple, MvaTuple, MVA_DATA)
#undef VAR
#undef MVA_DATA

}
}
