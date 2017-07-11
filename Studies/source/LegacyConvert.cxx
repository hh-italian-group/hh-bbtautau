/*! Study for Mva Training
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <numeric>

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "hh-bbtautau/Studies/include/MvaTuple.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"


struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_file);
    REQ_ARG(std::string, output_file);
};

namespace analysis {
namespace mva_study{

#define MVA_DATA() \
    VAR(UInt_t, NTrees) \
    VAR(double, shrinkage) \
    VAR(double, BaggedSampleFraction) \
    VAR(double, MaxDepth) \
    VAR(double, MinNodeSize) \
    VAR(double, ROCIntegral) \
    VAR(std::vector<double>, roc_value) \
    VAR(std::vector<int>, roc_mass) \
    VAR(std::vector<double>, KS_value) \
    VAR(std::vector<int>, KS_type) \
    VAR(std::vector<int>, KS_mass) \
    VAR(std::vector<double>, position) \
    VAR(std::vector<double>, importance) \
    VAR(std::vector<std::string>, var_name) \
    VAR(double, cut) \
    VAR(double, significance) \
    VAR(double, err_significance) \
    VAR(std::string, name) \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(legacy, MvaResults, MvaTuple, MVA_DATA, "mva_result")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(legacy, MvaTuple, MVA_DATA)
#undef VAR
#undef MVA_DATA

class LegacyConvert{
public:

    LegacyConvert(const Arguments& _args): args(_args), infile(root_ext::OpenRootFile(args.input_file())),
        outfile(root_ext::CreateRootFile(args.output_file()))
    {
    }

    template<typename T>
    static unsigned FindIndex(const std::vector<T>& points, T value)
    {
        for(unsigned n = 0; n < points.size(); ++n) {
            if(std::abs<double>(points.at(n) - value) < 1e-14)
                return n;
        }
        throw exception("Index not found.");
    }

    void Run()
    {
        legacy::MvaTuple legacyTuple("mva_result", infile.get(), true);
        mva_study::MvaTuple tuple("mva_result", outfile.get(), false);

        static std::vector<unsigned> TreeSteps = { 300, 600, 900, 1200 };
        static std::vector<double> shrinkageSteps = { 0.1, 0.4, 0.7, 1. };
        static std::vector<double> MaxDepthSteps = { 2, 3, 4, 5 };
        static std::vector<double> MinNodeSizeSteps = { 0.01, 0.05, 0.09 };
        static std::vector<double> BaggedSampleFractionSteps = { 0.5, 0.75, 1. };

        for(const legacy::MvaResults& legacyResult : legacyTuple) {
            mva_study::MvaResults& results = tuple();
            results.param_names.push_back("NTrees");
            results.param_values.push_back(legacyResult.NTrees);
            results.param_positions.push_back(FindIndex(TreeSteps, legacyResult.NTrees));
            results.param_names.push_back("shrinkage");
            results.param_values.push_back(legacyResult.shrinkage);
            results.param_positions.push_back(FindIndex(shrinkageSteps, legacyResult.shrinkage));
            results.param_names.push_back("MaxDepth");
            results.param_values.push_back(legacyResult.MaxDepth);
            results.param_positions.push_back(FindIndex(MaxDepthSteps, legacyResult.MaxDepth));
            results.param_names.push_back("MinNodeSize");
            results.param_values.push_back(legacyResult.MinNodeSize);
            results.param_positions.push_back(FindIndex(MinNodeSizeSteps, legacyResult.MinNodeSize));
            results.param_names.push_back("BaggedSampleFraction");
            results.param_values.push_back(legacyResult.BaggedSampleFraction);
            results.param_positions.push_back(FindIndex(BaggedSampleFractionSteps, legacyResult.BaggedSampleFraction));
            results.ROCIntegral = legacyResult.ROCIntegral;
            results.KS_mass = legacyResult.KS_mass;
            results.KS_type = legacyResult.KS_type;
            results.KS_value = legacyResult.KS_value;
            results.var_name = legacyResult.var_name;
            results.importance = legacyResult.importance;
            for(double pos : legacyResult.position)
                results.position.push_back(static_cast<unsigned>(pos));
            results.roc_mass = legacyResult.roc_mass;
            results.roc_value = legacyResult.roc_value;
            results.optimal_cut = legacyResult.cut;
            results.significance = legacyResult.significance;
            results.significance_err = legacyResult.err_significance;
            results.name = legacyResult.name;

            tuple.Fill();
        }

        tuple.Write();
    }

private:
    Arguments args;
    std::shared_ptr<TFile> infile, outfile;
};

}
}



PROGRAM_MAIN(analysis::mva_study::LegacyConvert, Arguments) // definition of the main program function
