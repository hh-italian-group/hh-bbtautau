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
    TH1D_ENTRY(KS_type, 4, -2, 2)
    TH1D_ENTRY(cut, 200, -1, 1)
    TH1D_ENTRY(significance, 200, 0, 1000)
};

class MVAAnalyzer{
public:

    MVAAnalyzer(const Arguments& _args): args(_args), outfile(root_ext::CreateRootFile(args.output_file())), anaData(outfile)
    {
    }

//    void CreatePositionHisto(std::map<std::string, std::shared_ptr<TH1D>>& histo, const std::vector<double> position, std::vector<std::string> var_name){
//        for (size_t j=0; j< var_name.size(); j++){
//            auto name = var_name[j];
//            if (!histo.count(name)){
//                histo[name] =  std::make_shared<TH1D>((name+"_position").c_str(),(name+"_position").c_str(), var_name.size(), 0.5, var_name.size() + 0.5);
//            }
//            histo[name]->Fill(position[j]);
//        }
//    }

    void Run()
    {

        std::map<std::string, std::shared_ptr<TH1D>> position6, position9, positionKS;
        for (const auto& name : args.input_file()){
            std::cout<<"ciao"<<std::endl;
            std::shared_ptr<TFile> in_file(root_ext::OpenRootFile(name));
            ntuple::MvaTuple myTree("mva_result", in_file.get(), true);
            for(const ntuple::MvaResults& results : myTree) {
                 if (results.ROCIntegral < 0.6){
                     anaData.ROCIntegral(name+"<0.6").Fill(results.ROCIntegral);
                     anaData.shrinkage(name+"<0.6").Fill(results.shrinkage);
                     anaData.NTrees(name+"<0.6").Fill(results.NTrees);
                     anaData.BaggedSampleFraction(name+"<0.6").Fill(results.BaggedSampleFraction);
                     anaData.MaxDepth(name+"<0.6").Fill(results.MaxDepth);
                     anaData.MinNodeSize(name+"<0.6").Fill(results.MinNodeSize);
                     anaData.cut(name+"<0.6").Fill(results.cut);
                     anaData.significance(name+"<0.6").Fill(results.significance);
//                     CreatePositionHisto(position6, results.position, results.var_name);
                     for (size_t j=0; j< results.KS_mass.size(); j++){
                         anaData.KS_mass(name+"<0.6").Fill(results.KS_mass[j]);
                         anaData.KS_value(name+"<0.6").Fill(results.KS_value[j]);
                         anaData.KS_type(name+"<0.6").Fill(results.KS_type[j]);
                     }
                 }
                 if (results.ROCIntegral > 0.9){
                     anaData.ROCIntegral(name+">0.9").Fill(results.ROCIntegral);
                     anaData.shrinkage(name+">0.9").Fill(results.shrinkage);
                     anaData.NTrees(name+">0.9").Fill(results.NTrees);
                     anaData.BaggedSampleFraction(name+">0.9").Fill(results.BaggedSampleFraction);
                     anaData.MaxDepth(name+">0.9").Fill(results.MaxDepth);
                     anaData.MinNodeSize(name+">0.9").Fill(results.MinNodeSize);
                     anaData.cut(name+">0.9").Fill(results.cut);
                     anaData.significance(name+">0.9").Fill(results.significance);
//                     CreatePositionHisto(position9, results.position, results.var_name);
                     for (size_t j=0; j< results.KS_mass.size(); j++){
                         anaData.KS_mass(name+">0.9").Fill(results.KS_mass[j]);
                         anaData.KS_value(name+">0.9").Fill(results.KS_value[j]);
                         anaData.KS_type(name+">0.9").Fill(results.KS_type[j]);
                     }
                 }
                 bool load = false;
                     for (size_t j=0; j< results.KS_mass.size(); j++){
                         if ((results.KS_mass[j] == 2000 && results.KS_type[j]  == 1)|| (results.KS_mass[j] == 0 && results.KS_type[j]  == -1)){
                             if (results.KS_value[j]>0.05){
                                 if (!load){
                                 anaData.ROCIntegral(name+">0.8_KS").Fill(results.ROCIntegral);
                                 anaData.shrinkage(name+">0.8_KS").Fill(results.shrinkage);
                                 anaData.NTrees(name+">0.8_KS").Fill(results.NTrees);
                                 anaData.BaggedSampleFraction(name+">0.8_KS").Fill(results.BaggedSampleFraction);
                                 anaData.MaxDepth(name+">0.8_KS").Fill(results.MaxDepth);
                                 anaData.MinNodeSize(name+">0.8_KS").Fill(results.MinNodeSize);
                                 anaData.cut(name+">0.8_KS").Fill(results.cut);
                                 anaData.significance(name+">0.8_KS").Fill(results.significance);
//                                 CreatePositionHisto(positionKS, results.position, results.var_name);
                                 load = true;
                                 }
                                 anaData.KS_mass(name+">0.8_KS").Fill(results.KS_mass[j]);
                                 anaData.KS_value(name+">0.8_KS").Fill(results.KS_value[j]);
                                 anaData.KS_type(name+">0.8_KS").Fill(results.KS_type[j]);
                             }
                         }
                     }

             }
        }
//        auto directory6 = root_ext::GetDirectory(*outfile.get(), "<0.6");
//        for (const auto& var: position6)
//            root_ext::WriteObject(*var.second, directory6);
//        auto directory9 = root_ext::GetDirectory(*outfile.get(), ">0.9");
//        for (const auto& var: position9)
//            root_ext::WriteObject(*var.second, directory9);
//        auto directoryKS = root_ext::GetDirectory(*outfile.get(), ">0.8KS");
//        for (const auto& var: positionKS)
//            root_ext::WriteObject(*var.second, directoryKS);



    }

private:
    Arguments args;
    TFile myfile;
    std::shared_ptr<TTree> tree;
    std::shared_ptr<TFile> outfile;
    MvaData anaData;
};

}
}
}


PROGRAM_MAIN(analysis::mva_study::MVAAnalyzer, Arguments) // definition of the main program function
