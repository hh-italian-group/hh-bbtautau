/*! Study for Mva Training
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "hh-bbtautau/Analysis/include/MvaTuple.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"


struct Arguments { // list of all program arguments
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::vector<std::string>, input_file);
};

namespace analysis {
namespace mva_study{

namespace {

class MvaData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    TH1D_ENTRY(NTrees, 4, 150, 1350)
    TH1D_ENTRY(shrinkage, 4, -0.05, 1.15)
    TH1D_ENTRY(BaggedSampleFraction, 3, 0.375, 1.125)
    TH1D_ENTRY(MaxDepth, 4, 1.5, 5.5)
    TH1D_ENTRY(MinNodeSize, 3, -0.01, 0.11)
    TH1D_ENTRY(ROCIntegral, 100, 0., 1.1)
    TH1D_ENTRY(KS_value, 100, 0., 1.1)
    TH1D_ENTRY(KS_mass, 210, 0, 2100)
};

class MVAAnalyzer{
public:

    MVAAnalyzer(const Arguments& _args): args(_args), anaData(args.output_file())
    {
    }

    void Run()
    {

        for (const auto& name : args.input_file()){
            std::shared_ptr<TFile> in_file(root_ext::OpenRootFile(name));
            ntuple::MvaTuple myTree("mva_result", in_file.get(), true);
            for(const ntuple::MvaResults& results : myTree) {
                 if (results.ROCIntegral > 0.9){
                     bool load = false;
                     double count = 0 ;
                     for (size_t j=0; j< results.KS_mass.size(); j++){
                         if (results.KS_value[j]>0.05){
                             count ++;
                             if (!load){
                                 anaData.ROCIntegral(name).Fill(results.ROCIntegral);
                                 anaData.shrinkage(name).Fill(results.shrinkage);
                                 anaData.NTrees(name).Fill(results.NTrees);
                                 anaData.BaggedSampleFraction(name).Fill(results.BaggedSampleFraction);
                                 anaData.MaxDepth(name).Fill(results.MaxDepth);
                                 anaData.MinNodeSize(name).Fill(results.MinNodeSize);
                                 load = true;
                             }
                             anaData.KS_mass(name).Fill(results.KS_mass[j]);
                             anaData.KS_value(name).Fill(results.KS_value[j]);
                         }
                     }
                     std::cout<<count<<std::endl;
                 }
              }
        }

    }

private:
    Arguments args;
    TFile myfile;
    std::shared_ptr<TTree> tree;
    MvaData anaData;
};

}
}
}


PROGRAM_MAIN(analysis::mva_study::MVAAnalyzer, Arguments) // definition of the main program function
