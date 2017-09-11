/*! Study for Mva Training
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <numeric>

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "hh-bbtautau/Studies/include/MvaTuple.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "hh-bbtautau/Studies/include/MvaMethods.h"


struct Arguments { // list of all program arguments
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::vector<std::string>, input_file);
};

namespace analysis {
namespace mva_study{

namespace {

class MvaData : public root_ext::AnalyzerData {
public:
    explicit MvaData(std::shared_ptr<TFile> _outputFile, const std::string& directoryName = "") :
        AnalyzerData(_outputFile, directoryName)
    {
        histo.Emplace("shrinkage", 4, -0.05, 1.15);
        histo.Emplace("NTrees", 4, 150, 1350);
        histo.Emplace("BaggedSampleFraction", 3, 0.375, 1.125);
        histo.Emplace("MaxDepth", 4, 1.5, 5.5);
        histo.Emplace("MinNodeSize", 3, -0.01, 0.11);
        histo.Emplace("ROCIntegral", 50, 0., 1.);
        histo.Emplace("KS_value", 50, 0., 1.);
        histo.Emplace("KS_mass", 210, 0, 2100);
        histo.Emplace("KS_type", 4, -2, 2);
        histo.Emplace("cut", 200, -1, 1);
        histo.Emplace("significance", 100, 0, 30);
        ROC.Emplace("shrinkage", 4, -0.05, 1.15, 50, 0.5, 1.);
        ROC.Emplace("NTrees", 4, 150, 1350, 50, 0.5, 1.);
        ROC.Emplace("BaggedSampleFraction", 3, 0.375, 1.125, 50, 0.5, 1.);
        ROC.Emplace("MaxDepth", 4, 1.5, 5.5, 50, 0.5, 1);
        ROC.Emplace("MinNodeSize", 3, -0.01, 0.11, 50, 0.5, 1.);
        ROC.Emplace("relativeROC", 66, 245, 905, 20,0.7,1.3);
        ROC.Emplace("significance", 100, 0, 30, 50, 0.5, 1.);
        significance.Emplace("err", 100, 0, 1, 100, 0, 30);
        significance.Emplace("shrinkage", 4, -0.05, 1.15, 50, 0.5, 1.);
        significance.Emplace("NTrees", 4, 150, 1350, 50, 0.5, 1.);
        significance.Emplace("BaggedSampleFraction", 3, 0.375, 1.125, 50, 0.5, 1.);
        significance.Emplace("MaxDepth", 4, 1.5, 5.5, 50, 0.5, 1);
        significance.Emplace("MinNodeSize", 3, -0.01, 0.11, 50, 0.5, 1.);
    }
    TH1D_ENTRY(histo, 10, .5, 10.5)
    TH2D_ENTRY(ROC, 10, .5, 10.5, 50, 0.5, 1.)
    TH2D_ENTRY(significance, 10, .5, 10.5, 100, 0, 30)
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

    std::map<std::string, double> AveragePosition(const std::map<std::string, std::vector<VarRank>>& positions)
    {
        std::map<std::string, double> average;
        for(const auto& var : positions) {
            for(const auto& rank : var.second)
                average[var.first] += rank.position;
            average[var.first] /= var.second.size();
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
        static constexpr double KS_cut = 0.05;

        std::map<std::string, std::shared_ptr<TH1D>> histo_position;
        std::map<std::string, GridPoint> method_params;
        std::map<std::string, size_t> method_seed_count;

        std::map<std::string, std::vector<double>> roc;
        std::map<std::string, std::map<int, std::vector<PhysicalValue>>> significance;
        std::map<std::string, std::map<int, std::vector<double>>> optimal_cut;
        std::map<std::string, std::map<std::string, std::vector<VarRank>>> position;
        const size_t n_seeds = args.input_file().size();
        std::cout<<n_seeds<<std::endl;

        for (const auto& name : args.input_file()){
            std::shared_ptr<TFile> in_file(root_ext::OpenRootFile(name));
            MvaTuple myTree("mva_result", in_file.get(), true);
            for(const MvaResults& results : myTree) {
                const std::string method_name = results.name.substr(0, results.name.find_last_of('_'));
                const auto grid_point = GetGridPoint(results);
                std::cout<<grid_point.size()<<std::endl;
                const auto roc_integrals = GetRocIntegralMap(results);
                const auto KS_results = GetKSResultsMap(results);
                const auto ranking = GetRankingMap(results);
                const auto sign = GetOptimalSignificanceMap(results);

                if(!KS_results.count(SampleId::MassTot()))
                    throw exception("Missing total mass info");
                if(!KS_results.count(SampleId::Bkg()))
                    throw exception("Missing bkg info");
                if(KS_results.at(SampleId::MassTot()) <= KS_cut || KS_results.at(SampleId::Bkg()) <= KS_cut) continue;

                method_params[method_name] = grid_point;
                ++method_seed_count[method_name];
                std::cout<<method_seed_count[method_name]<<std::endl;

                roc[method_name].push_back(results.ROCIntegral);
                for(const auto& s : sign)
                {
                    significance[method_name][s.first].emplace_back(s.second.significance);
                    optimal_cut[method_name][s.first].push_back(s.second.cut);
                }
                for(const auto& integral : roc_integrals)
                     anaData.ROC("relativeROC").Fill(integral.first, integral.second/results.ROCIntegral);

                for(const auto& ks : KS_results) {
                    anaData.histo("KS_mass").Fill(ks.first.mass);
                    anaData.histo("KS_type").Fill(static_cast<int>(ks.first.sampleType));
                    anaData.histo("KS_value").Fill(ks.second);
                }

                for(const auto& rank : ranking)
                    position[method_name][rank.first].push_back(rank.second);
             }
        }


        for(const auto& method : method_seed_count) {
            if(method.second < n_seeds - 1) continue;
            const std::string method_name = method.first;
            const double roc_value = std::accumulate(roc[method_name].begin(), roc[method_name].end(), 0.) / roc[method_name].size();
            anaData.histo("ROCIntegral").Fill(roc_value);
            std::map<int, PhysicalValue> significance_value;
            for (const auto& sign : significance[method_name]){
                significance_value[sign.first] = std::accumulate(sign.second.begin(), sign.second.end(), PhysicalValue::Zero) / PhysicalValue(sign.second.size(), 0);
            }

            for(const auto& param : method_params[method_name]) {
                const double value = param.second.value;
                anaData.ROC(param.first).Fill(value, roc_value);
                anaData.significance(param.first).Fill(value, significance_value[mass_tot.mass].GetValue());
                anaData.histo(param.first).Fill(value);
            }
            anaData.ROC("significance").Fill(significance_value[mass_tot.mass].GetValue(), roc_value);
            anaData.significance("err").Fill(significance_value[mass_tot.mass].GetRelativeStatisticalError(), significance_value[mass_tot.mass].GetValue());
            anaData.histo("significance").Fill(significance_value[mass_tot.mass].GetValue());

            std::map<int, double> cut_value;
            for (const auto& cut : optimal_cut[method_name]){
                cut_value[cut.first] = std::accumulate(cut.second.begin(), cut.second.end(), 0.) / cut.second.size();
            }
            anaData.histo("cut").Fill(cut_value[mass_tot.mass]);
            const auto average = AveragePosition(position[method_name]);
            CreatePositionHisto(histo_position, average);
        }

        auto directory = root_ext::GetDirectory(*outfile.get(), "RankingPosition");
        for (const auto& var: histo_position)
            root_ext::WriteObject(*var.second, directory);
        CreateRanking(histo_position);
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
