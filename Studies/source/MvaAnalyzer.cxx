/*! Study for Mva Training
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "hh-bbtautau/Analysis/include/Lester_mt2_bisect.h"
#include "hh-bbtautau/Analysis/include/MvaConfiguration.h"
#include "hh-bbtautau/Analysis/include/MvaMethods.h"
#include "hh-bbtautau/Analysis/include/MvaVariables.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TFile.h"
#include <fstream>
#include <random>


struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_file);
    REQ_ARG(std::string, output_file);
    OPT_ARG(size_t, seed, std::numeric_limits<uint_fast32_t>::max());
};

namespace analysis {
namespace mva_study{

namespace {

class MVAAnalyzer{
public:

    MVAAnalyzer(const Arguments& _args): args(_args), outfile(root_ext::CreateRootFile(args.output_file()))

    {
    }

    void Run()
    {

        TFile* myfile = new TFile((args.input_file()).c_str());
        TTree* myTree = (TTree*)myfile->Get("mva_result");
        double myShrinkage, myBaggedSampleFraction, myMaxDepth, myMinNodeSize, myROCIntegral;
        UInt_t myTrees;
        std::vector<double> *myKSvalue, *myKSmass, *myPosition, *myImportance;
        std::vector<std::string> *myVarName;

        gROOT->ProcessLine("#include <vector>");
        myTree->SetBranchAddress("NTrees", &myTrees);
        myTree->SetBranchAddress("shrinkage", &myShrinkage);
        myTree->SetBranchAddress("BaggedSampleFraction", &myBaggedSampleFraction);
        myTree->SetBranchAddress("MaxDepth", &myMaxDepth);
        myTree->SetBranchAddress("MinNodeSize", &myMinNodeSize);
        myTree->SetBranchAddress("ROCIntegral", &myROCIntegral);
        myTree->SetBranchAddress("KS_value", &myKSvalue);
        myTree->SetBranchAddress("KS_mass", &myKSmass);
        myTree->SetBranchAddress("position", &myPosition);
        myTree->SetBranchAddress("importance", &myImportance);
        myTree->SetBranchAddress("var_name", &myVarName);
        std::cout<<"ciao"<<std::endl;
        Long64_t nentries = myTree->GetEntries();


        for (Long64_t i=0;i<nentries;i++) {
             myTree->GetEntry(i);
             std::cout<<myKSmass->size()<<std::endl;
             std::cout<<myKSvalue->size()<<std::endl;
             std::cout<<myPosition->size()<<std::endl;
             std::cout<<myImportance->size()<<std::endl;
             for (size_t j=0; j<myImportance->size(); j++){
                 std::cout<<(*myVarName)[j]<<"  "<<(*myPosition)[j]<<"  "<<(*myImportance)[j]<<std::endl;
             }
             std::cout<<myVarName->size()<<std::endl<<std::endl;

          }






    }

private:
    Arguments args;
    TFile myfile;
    std::shared_ptr<TFile> outfile;
};

}
}
}


PROGRAM_MAIN(analysis::mva_study::MVAAnalyzer, Arguments) // definition of the main program function
//./run.sh MvaTraining --input_path ~/Desktop/tuples --output_file myfile.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name muTau --number_events 10000000

