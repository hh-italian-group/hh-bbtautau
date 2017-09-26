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
    REQ_ARG(std::string, which_test);
    REQ_ARG(double, cut);
    REQ_ARG(bool, cross_validation);
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
        histo.Emplace("shrinkage", 8, 0.04, 0.36);
        histo.Emplace("NTrees", 4, 150, 1350);
        histo.Emplace("BaggedSampleFraction", 3, 0.375, 1.125);
        histo.Emplace("MaxDepth", 4, 1.5, 5.5);
        histo.Emplace("MinNodeSize", 3, -0.01, 0.11);
        histo.Emplace("ROCIntegral_testing", 50, 0., 1.);
        histo.Emplace("ROCIntegral_training", 50, 0., 1.);
        histo.Emplace("KS_value", 50, 0., 1.);
        histo.Emplace("chi_value", 50, 0., 1.);
        histo.Emplace("cut", 200, -1, 1);
        histo.Emplace("significance", 100, 0, 30);
        histo.Emplace("nCuts", 17, 0, 680);

        ROC.Emplace("shrinkage", 5, 0.025, 0.275, 30, 0.7, 1.);
        ROC.Emplace("NTrees", 3, 250, 1150, 30, 0.7, 1.);
        ROC.Emplace("BaggedSampleFraction", 3, 0.375, 1.125, 30, 0.7, 1.);
        ROC.Emplace("MaxDepth", 4, 1.5, 5.5, 30, 0.7, 1.);
        ROC.Emplace("MinNodeSize", 3, 0., 0.06, 30, 0.7, 1.);
        ROC.Emplace("relativeROC_training", 66, 245, 905, 30, 0.7, 1.);
        ROC.Emplace("relativeROC_testing", 66, 245, 905, 30, 0.7, 1.);
        ROC.Emplace("significance", 100, 0, 30, 30, 0.7, 1.);
        ROC.Emplace("nCuts", 3, 0, 120, 30, 0.7, 1.);
        ROC.Emplace("mass", 3, 0, 120, 30, 0.7, 1.);

        sigmaROC.Emplace("shrinkage", 5, 0.025, 0.275, 30, 0., 0.1);
        sigmaROC.Emplace("NTrees", 4, 150, 1350, 30, 0., 0.1);
        sigmaROC.Emplace("BaggedSampleFraction", 3, 0.375, 1.125, 30, 0., 0.1);
        sigmaROC.Emplace("MaxDepth", 4, 1.5, 5.5, 30, 0., 0.1);
        sigmaROC.Emplace("MinNodeSize", 3, -0.01, 0.11, 30, 0., 0.1);
        sigmaROC.Emplace("relativeROC_training", 66, 245, 905, 30, 0., 0.1);
        sigmaROC.Emplace("relativeROC_testing", 66, 245, 905, 30, 0., 0.1);
        sigmaROC.Emplace("significance", 100, 0, 30, 30, 0., 0.1);
        sigmaROC.Emplace("nCuts", 3, 0, 120, 30, 0., 0.1);
        sigmaROC.Emplace("mass", 3, 0, 120, 30, 0., 0.1);

        significance.Emplace("err", 100, 0, 1, 100, 0, 30);
        significance.Emplace("shrinkage", 8, 0.04, 0.36, 100, 0, 30);
        significance.Emplace("NTrees", 4, 150, 1350, 100, 0, 30);
        significance.Emplace("BaggedSampleFraction", 3, 0.375, 1.125, 100, 0, 30);
        significance.Emplace("MaxDepth", 4, 1.5, 5.5, 100, 0, 30);
        significance.Emplace("MinNodeSize", 3, -0.01, 0.11, 100, 0, 30);
        significance.Emplace("nCuts", 17, 0, 680, 100, 0, 30);
    }
    TH1D_ENTRY(histo, 10, .5, 10.5)
    TH2D_ENTRY(ROC, 10, .5, 10.5, 50, 0.5, 1.)
    TH2D_ENTRY(sigmaROC, 10, .5, 10.5, 50, 0.5, 1.)
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
        double cut = args.cut();

        std::map<std::string, std::shared_ptr<TH1D>> histo_position;
        std::map<std::string, GridPoint> method_params;
        std::map<std::string, size_t> method_seed_count;

        std::map<std::string, std::map<ChannelSampleIdSpin, std::vector<PhysicalValue>>> significance;
        std::map<std::string, std::map<ChannelSampleIdSpin, std::vector<double>>> optimal_cut;
        std::map<std::string, std::map<std::string, std::vector<VarRank>>> position;
        std::map<std::string, std::map<ChannelSampleIdSpin, double>> roc_training, roc_testing;
        std::map<std::string, std::map<ChannelSampleIdSpin, std::vector<double>>> vec_roc_training, vec_roc_testing;

        const size_t n_seeds = args.input_file().size();
        std::cout<<n_seeds<<std::endl;

        for (const auto& name : args.input_file()){
            std::shared_ptr<TFile> in_file(root_ext::OpenRootFile(name));
            MvaTuple myTree("mva_result", in_file.get(), true);
            for(const MvaResults& results : myTree) {
                const std::string method_name = results.name.substr(0, results.name.find_last_of('_'));
                const auto grid_point = GetGridPoint(results);               
                const auto KS_results = GetKSResultsMap(results);
                const auto chi_results = GetChiResultsMap(results);
                const auto ranking = GetRankingMap(results);
                const auto sign = GetOptimalSignificanceMap(results);

                if(!KS_results.count(AllSgn))
                    throw exception("Missing total mass info");
                if(!KS_results.count(AllBkg))
                    throw exception("Missing bkg info");

                if (args.which_test() == "KS")
                    if(KS_results.at(AllSgn) <= cut || KS_results.at(AllBkg) <= cut) continue;

                if (args.which_test() == "chi")
                    if(chi_results.at(AllSgn) <= cut || chi_results.at(AllBkg) <= cut) continue;
                method_params[method_name] = grid_point;
                ++method_seed_count[method_name];
                std::cout<<method_seed_count[method_name]<<std::endl;

                roc_training[method_name] = GetRocTrainingIntegralMap(results);
                roc_testing[method_name] = GetRocTestingIntegralMap(results);

                for (const auto& entry: roc_training[method_name])
                    vec_roc_training[method_name][entry.first].push_back(entry.second);
                for (const auto& entry: roc_testing[method_name])
                    vec_roc_testing[method_name][entry.first].push_back(entry.second);

                for(const auto& integral : roc_training[method_name]){
                    if (!integral.first.IsAllChannel() || !integral.first.IsAllSpin()) continue;
                    anaData.ROC("relativeROC_training").Fill(integral.first.sample_id.mass, integral.second/roc_training[method_name][AllSgn]);
                }

                for(const auto& integral : roc_testing[method_name]){
                    if (!integral.first.IsAllChannel() || !integral.first.IsAllSpin()) continue;
                    anaData.ROC("relativeROC_testing").Fill(integral.first.sample_id.mass, integral.second/roc_testing[method_name][AllSgn]);
                }


                for(const auto& s : sign)
                {
                    significance[method_name][s.first].emplace_back(s.second.significance);
                    optimal_cut[method_name][s.first].push_back(s.second.cut);
                }

                for(const auto& ks : KS_results) {
                    anaData.histo("KS_value").Fill(ks.second);
                }
                for(const auto& chi : chi_results) {
                    anaData.histo("chi_value").Fill(chi.second);
                }

                for(const auto& rank : ranking)
                    position[method_name][rank.first].push_back(rank.second);
             }
            std::cout<<std::endl;

        }

        std::cout<<"Quanti metodi?"<<method_seed_count.size()<<std::endl;
        int i = 0;

        std::map<std::string, std::map<ChannelSampleIdSpin, double>> roc_training_value, roc_testing_value;
        std::map<std::string, std::map<ChannelSampleIdSpin, double>> roc_training_err, roc_testing_err;

        MvaTuple mva_tuple(outfile.get(), false);

        for(const auto& method : method_seed_count) {
            if (args.cross_validation())
                if(method.second < (n_seeds-1)) continue;
            const std::string method_name = method.first;
            std::cout<<method_name<<std::endl;

            mva_tuple().name = method_name;

            for (const auto& value: vec_roc_training[method_name]){
                roc_training_value[method_name][value.first] = std::accumulate(value.second.begin(), value.second.end(), 0.) / value.second.size();
                roc_training_err[method_name][value.first] = std::sqrt(stat_estimators::Variance(value.second));

                std::cout<<roc_training_value[method_name][value.first]<<"pm"<<roc_training_err[method_name][value.first]<<std::endl;
                mva_tuple().err_roc_training.push_back(roc_training_err[method_name][value.first]);
                mva_tuple().roc_training_channel.push_back(value.first.channel);
                mva_tuple().roc_training_mass.push_back(value.first.sample_id.mass);
                mva_tuple().roc_training_type.push_back(static_cast<int>(value.first.sample_id.sampleType));
                mva_tuple().roc_training_value.push_back(roc_training_value[method_name][value.first]);
                mva_tuple().roc_training_spin.push_back(value.first.spin);
            }

            for (const auto& value: vec_roc_testing[method_name]){
                roc_testing_value[method_name][value.first] = std::accumulate(value.second.begin(), value.second.end(), 0.) / value.second.size();
                roc_testing_err[method_name][value.first] = stat_estimators::Variance(value.second);
                mva_tuple().err_roc_testing.push_back(roc_testing_err[method_name][value.first]);
                mva_tuple().roc_testing_channel.push_back(value.first.channel);
                mva_tuple().roc_testing_mass.push_back(value.first.sample_id.mass);
                mva_tuple().roc_testing_type.push_back(static_cast<int>(value.first.sample_id.sampleType));
                mva_tuple().roc_testing_value.push_back(roc_testing_value[method_name][value.first]);
                mva_tuple().roc_testing_spin.push_back(value.first.spin);
            }

            double roc = roc_testing_value[method_name][AllSgn];
            double err_roc = roc_testing_err[method_name][AllSgn];
            anaData.histo("ROCIntegral").Fill(roc);
            std::map<ChannelSampleIdSpin, PhysicalValue> significance_value;
            for (const auto& sign : significance[method_name]){
                significance_value[sign.first] = std::accumulate(sign.second.begin(), sign.second.end(), PhysicalValue::Zero) / PhysicalValue(sign.second.size(), 0);
            }

            for(const auto& param : method_params[method_name]) {
                const double value = param.second.value;
                anaData.ROC(param.first).Fill(value, roc);
                anaData.sigmaROC(param.first).Fill(value, err_roc);
                anaData.significance(param.first).Fill(value, significance_value[AllSgn].GetValue());
                anaData.histo(param.first).Fill(value);
            }
            anaData.ROC("significance").Fill(significance_value[AllSgn].GetValue(), roc);
            anaData.significance("err").Fill(significance_value[AllSgn].GetRelativeStatisticalError(), significance_value[AllSgn].GetValue());
            anaData.histo("significance").Fill(significance_value[AllSgn].GetValue());

            std::map<ChannelSampleIdSpin, double> cut_value;
            for (const auto& cut : optimal_cut[method_name]){
                cut_value[cut.first] = std::accumulate(cut.second.begin(), cut.second.end(), 0.) / cut.second.size();
            }
            anaData.histo("cut").Fill(cut_value[AllSgn]);
            const auto average = AveragePosition(position[method_name]);
            CreatePositionHisto(histo_position, average);
            i++;

            mva_tuple.Fill();
        }

        mva_tuple.Write();

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
