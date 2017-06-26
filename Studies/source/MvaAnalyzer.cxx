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

    MVAAnalyzer(const Arguments& _args): args(_args), in_file(root_ext::OpenRootFile(args.input_file())),
        out_file(root_ext::CreateRootFile(args.output_file()))

    {
    }

    void Add(std::map<std::string, std::set<double>>& right_methods, const std::string& name, const double& val){
        if (!right_methods[name].count(val)){
            right_methods[name].insert(val);
        }
    }


    void Run()
    {

        auto myTree = (TTree*)(in_file->Get("mva_result"));
        double myShrinkage, myBaggedSampleFraction, myMaxDepth, myMinNodeSize, myROCIntegral;
        UInt_t myNTrees;
        std::vector<double> *myKSvalue, *myKSmass, *myPosition, *myImportance;
        std::vector<std::string> *myVarName;
        std::map<std::string, std::shared_ptr<TH1D>> histo_properties;

        gROOT->ProcessLine("#include <vector>");
        myTree->SetBranchAddress("NTrees", &myNTrees);
        histo_properties["NTrees"] = std::make_shared<TH1D>("NTrees", "NTrees",100,0,1000);
        myTree->SetBranchAddress("shrinkage", &myShrinkage);
        histo_properties["shrinkage"] = std::make_shared<TH1D>("shrinkage", "shrinkage",11,0,1.1);
        myTree->SetBranchAddress("BaggedSampleFraction", &myBaggedSampleFraction);
        histo_properties["BaggedSampleFraction"] = std::make_shared<TH1D>("BaggedSampleFraction", "BaggedSampleFraction",3,0.5,1.25);
        myTree->SetBranchAddress("MaxDepth", &myMaxDepth);
        histo_properties["MaxDepth"] = std::make_shared<TH1D>("MaxDepth", "MaxDepth",3,2,5);
        myTree->SetBranchAddress("MinNodeSize", &myMinNodeSize);
        histo_properties["MinNodeSize"] = std::make_shared<TH1D>("MinNodeSize", "MinNodeSize",11,0,0.1);
        myTree->SetBranchAddress("ROCIntegral", &myROCIntegral);
        histo_properties["ROCIntegral"] = std::make_shared<TH1D>("ROCIntegral", "ROCIntegral",100,0,1.1);
        myTree->SetBranchAddress("KS_value", &myKSvalue);
        histo_properties["KS_value"] = std::make_shared<TH1D>("KS_value", "KS_value",100,0,1.1);
        myTree->SetBranchAddress("KS_mass", &myKSmass);
        histo_properties["KS_mass"] = std::make_shared<TH1D>("KS_mass", "KS_mass",210,0,2100);
        myTree->SetBranchAddress("position", &myPosition);
        myTree->SetBranchAddress("importance", &myImportance);
        myTree->SetBranchAddress("var_name", &myVarName);
        std::cout<<histo_properties.size()<<std::endl;

        Long64_t nentries = myTree->GetEntries();
        std::map<std::string, std::set<double>> right_methods;

        for (Long64_t i = 0 ; i < nentries ; i++) {
             myTree->GetEntry(i);
             if (myROCIntegral > 0.9){
                 bool load = false;
                 double count = 0 ;
                 for (size_t j=0; j< myKSmass->size(); j++){
                     if ((*myKSmass)[j] != 2000 && (*myKSvalue)[j]>0.05){
                         count ++;
                         if (!load){
                             histo_properties.at("ROCIntegral")->Fill(myROCIntegral);
                             Add(right_methods, "shrinkage", myShrinkage);
                             histo_properties.at("shrinkage")->Fill(myShrinkage);
                             Add(right_methods, "NTrees", myNTrees);
                             histo_properties.at("NTrees")->Fill(myNTrees);
                             Add(right_methods, "BaggedSampleFraction", myBaggedSampleFraction);
                             histo_properties.at("BaggedSampleFraction")->Fill(myBaggedSampleFraction);
                             Add(right_methods, "MaxDepth", myMaxDepth);
                             histo_properties.at("MaxDepth")->Fill(myMaxDepth);
                             Add(right_methods, "MinNodeSize", myMinNodeSize);
                             histo_properties.at("MinNodeSize")->Fill(myMinNodeSize);
                             load = true;
                         }
                         histo_properties.at("KS_mass")->Fill((*myKSmass)[j]);
                         histo_properties.at("KS_value")->Fill((*myKSvalue)[j]);
                     }
                 }
                 std::cout<<count<<std::endl;
             }
          }

        for (const auto& met : right_methods)
            std::cout << met.first << " " << met.second.size() <<std::endl;

        for (const auto& histo : histo_properties)
            root_ext::WriteObject(*histo.second, out_file.get());





    }

private:
    Arguments args;
    TFile myfile;
    std::shared_ptr<TFile> in_file;
    std::shared_ptr<TFile> out_file;
    std::shared_ptr<TTree> tree;
};

}
}
}


PROGRAM_MAIN(analysis::mva_study::MVAAnalyzer, Arguments) // definition of the main program function
//./run.sh MvaTraining --input_path ~/Desktop/tuples --output_file myfile.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name muTau --number_events 10000000

