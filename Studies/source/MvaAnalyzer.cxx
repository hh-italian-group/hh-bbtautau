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
    TH1D_ENTRY(ROCIntegral, 100, 0.5, 1.1)
    TH1D_ENTRY(KS_value, 100, 0., 1.1)
    TH1D_ENTRY(KS_mass, 210, 0, 2100)
    TH1D_ENTRY(KS_type, 4, -2, 2)
    TH1D_ENTRY(cut, 200, -1, 1)
    TH1D_ENTRY(significance, 200, 0, 1000)
    TH2D_ENTRY(shrinkage_ROC,  4, -0.05, 1.15, 100, 0.8, 1.)
    TH2D_ENTRY(NTrees_ROC,  4, 150, 1350, 100, 0.5, 1.1)
    TH2D_ENTRY(BaggedSampleFraction_ROC,  3, 0.375, 1.125, 100, 0.8, 1.)
    TH2D_ENTRY(MaxDepth_ROC,  4, 1.5, 5.5, 100, 0.8, 1.)
    TH2D_ENTRY(MinNodeSize_ROC,  3, -0.01, 0.11, 100, 0.8, 1.)
    TH2D_ENTRY(mass_relativeROC, 8, 245, 325, 20,0.9,1.1)

};

class MVAAnalyzer{
public:

    MVAAnalyzer(const Arguments& _args): args(_args), outfile(root_ext::CreateRootFile(args.output_file())), anaData(outfile)
    {
    }

    void CreatePositionHisto(std::map<std::string, std::shared_ptr<TH1D>>& histo, const std::map<std::string, double>& average)
    {
        for (const auto& name :  average){
            if (!histo.count(name.first)){
                histo[name.first] =  std::make_shared<TH1D>((name.first+"_position").c_str(),(name.first+"_position").c_str(), average.size()*2, 0.5, average.size() + 5.);
            }
            histo[name.first]->Fill(name.second);
        }
    }

    std::map<std::string, double> AveragePosition(const std::vector<std::vector<double>>& position, const std::vector<std::vector<std::string>>& var_name)
    {
        std::map<std::string, double> average;
        for(size_t i = 0; i<position.size(); ++i){
            for (size_t j = 0; j< position[i].size(); ++j){
                average[var_name.at(i).at(j)] = average[var_name.at(i).at(j)] + position.at(i).at(j);
            }
        }
        for(auto& av : average){
            av.second = av.second / position.size();
        }
        return average;
    }

    void CreateRanking(const std::map<std::string, std::shared_ptr<TH1D>>& histo)
    {
        std::vector<std::pair<std::string,double>> ranking;
        for(const auto& name : histo)
        {
            ranking.emplace_back(name.first, name.second->GetMean());
        }
        std::sort(ranking.begin(), ranking.end(), [](auto el1, auto el2){
            return el1.second < el2.second;
        });

        for(const auto& rank :  ranking){
            std::cout<<rank.first<<"    "<<rank.second<<std::endl;
        }
    }

    void Run()
    {

        std::map<std::string, std::shared_ptr<TH1D>> positionKS;
        std::map<std::string, std::map<std::string, std::map<double, double>>> param_seed;
        std::map<std::string, double> roc;
        std::map<std::string, std::vector<std::vector<double>>> position;
        std::map<std::string, std::vector<std::vector<std::string>>> var_name;

        for (const auto& name : args.input_file()){
//            std::cout<<"ciao"<<std::endl;
            std::shared_ptr<TFile> in_file(root_ext::OpenRootFile(name));
            ntuple::MvaTuple myTree("mva_result", in_file.get(), true);
            for(const ntuple::MvaResults& results : myTree) {
                std::string method_name = "Grad_BaggedSampleFraction_"+std::to_string(results.BaggedSampleFraction)+"_MaxDepth_"+std::to_string(results.MaxDepth)+
                        "_MinNodeSize_"+std::to_string(results.MinNodeSize)+"_NTrees_"+std::to_string(results.NTrees)+"_shrinkage"+std::to_string(results.shrinkage);
                 bool load = false;
                 if (results.shrinkage > 0.41 || results.MaxDepth > 3) continue;
                 for (size_t j=0; j< results.KS_mass.size(); j++){
                     if ((results.KS_mass[j] == 2000 && results.KS_type[j]  == 1) || (results.KS_mass[j] == 0 && results.KS_type[j]  == -1)){
                         if (results.KS_value[j]>0.05){
                             if (!load){
                                 param_seed[method_name]["shrinkage"][results.shrinkage]++;
                                 param_seed[method_name]["NTrees"][results.NTrees]++;
                                 param_seed[method_name]["BaggedSampleFraction"][results.BaggedSampleFraction]++;
                                 param_seed[method_name]["MaxDepth"][results.MaxDepth]++;
                                 param_seed[method_name]["MinNodeSize"][results.MinNodeSize]++;
                                 roc[method_name] = results.ROCIntegral;
                                 position[method_name].push_back(results.position);
                                 var_name[method_name].push_back(results.var_name);
                                 for (size_t i = 0; i<results.roc_mass.size(); ++i){
                                     anaData.mass_relativeROC().Fill(results.roc_mass[i], results.roc_value[i]/results.ROCIntegral);
                                 }
                             load = true;
                             }
                             anaData.KS_mass().Fill(results.KS_mass[j]);
                             anaData.KS_value().Fill(results.KS_value[j]);
                             anaData.KS_type().Fill(results.KS_type[j]);
                         }
                     }
                 }

             }
        }

        for(const auto& method : param_seed){
            bool passed = true;
            for (const auto& param : method.second){
                for (const auto& value : param.second){
                    if (value.second < (args.input_file().size()-1)){
                        passed = false;
                        continue;
                    }
                    if (param.first == "shrinkage"){
                        anaData.shrinkage().Fill(value.first);
                        anaData.shrinkage_ROC().Fill(value.first, roc[method.first]);
                    }
                    if (param.first == "NTrees"){
                        anaData.NTrees().Fill(value.first);
                        anaData.NTrees_ROC().Fill(value.first, roc[method.first]);
                    }
                    if (param.first == "BaggedSampleFraction"){
                        anaData.BaggedSampleFraction().Fill(value.first);
                        anaData.BaggedSampleFraction_ROC().Fill(value.first, roc[method.first]);
                    }
                    if (param.first == "MaxDepth"){
                        anaData.MaxDepth().Fill(value.first);
                        anaData.MaxDepth_ROC().Fill(value.first, roc[method.first]);
                    }
                    if (param.first == "MinNodeSize"){
                        anaData.MinNodeSize().Fill(value.first);
                        anaData.MinNodeSize_ROC().Fill(value.first, roc[method.first]);
                    }
                }
            }
            if(!passed) continue;
            auto average = AveragePosition(position[method.first], var_name[method.first]);
            CreatePositionHisto(positionKS, average);
        }

        auto directory = root_ext::GetDirectory(*outfile.get(), "RankingPosition");
        for (const auto& var: positionKS)
            root_ext::WriteObject(*var.second, directory);
        CreateRanking(positionKS);
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
