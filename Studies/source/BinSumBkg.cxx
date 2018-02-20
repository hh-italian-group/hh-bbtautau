/*! Post-processing of the analysis results.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/EventAnalyzerCore.h"
#include "hh-bbtautau/Analysis/include/AnaTuple.h"
#include "hh-bbtautau/Analysis/include/EventAnalyzerDataCollection.h"
#include "hh-bbtautau/Analysis/include/StackedPlotsProducer.h"
#include "hh-bbtautau/Analysis/include/LimitsInputProducer.h"
#include "AnalysisTools/Core/include/SmartHistogram.h"

namespace analysis {

struct AnalyzerArguments {
    REQ_ARG(std::string, input);

};

class BinSumBkg  {
public:
    using Hist = SmartHistogram<TH1D>;
    using HistPtr = std::shared_ptr<Hist>;

    BinSumBkg(const AnalyzerArguments& _args) : args(_args)
    {

    }

    void Run()
    {
        auto inputFile = root_ext::OpenRootFile(args.input());

        std::vector<std::string> channels = {"eTau", "muTau", "tauTau"};
        std::vector<std::string> categories = {"res1b", "res2b", "boosted"};
        std::vector<std::string> backgrounds = {"TT", "DY_0b", "DY_1b", "DY_2b", "W", "WW","WZ", "ZZ", "ZH", "EWK", "tW" ,"QCD"};

        for(const auto& channel : channels){
            for(const auto& category : categories){

                const std::string hist_dir_name_data = channel+"_"+category+"/data_obs";
                auto hist_data = std::shared_ptr<TH1>
                        (root_ext::TryReadObject<TH1>(*inputFile,hist_dir_name_data));
                if (!hist_data){
                    std::cout << "WARNING! NO bkg histogram found for hist dir name data: "<< hist_dir_name_data <<  std::endl;
                    continue;
                }

                for(const auto& bkg : backgrounds){
                    const std::string hist_dir_name = channel+"_"+category+"/"+bkg;


                    auto hist_bkg = std::shared_ptr<TH1>
                            (root_ext::TryReadObject<TH1>(*inputFile,hist_dir_name));
                    if (!hist_bkg){
                        std::cout << "WARNING! NO bkg histogram found for hist dir name: "<< hist_dir_name <<  std::endl;
                        continue;
                    }

                    if(bkg_sum) {
                        bkg_sum->Add(hist_bkg.get(), 1.);
                    } else {
                        bkg_sum = std::make_shared<Hist>(*hist_bkg);
                    }

                }//loop bkgs

                std::cout << "channel: " << channel << ", category: " << category << ", data entries: " << hist_data->GetEntries() <<
                             ", sum bkg entries: " << bkg_sum->GetEntries() << std::endl;

            } //loop categories
        }//loop channels

    }





private:
    AnalyzerArguments args;
    HistPtr bkg_sum;

};

} // namespace analysis

PROGRAM_MAIN(analysis::BinSumBkg, analysis::AnalyzerArguments)
