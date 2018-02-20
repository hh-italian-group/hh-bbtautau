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
    using Hist = TH1D;
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

                HistPtr bkg_sum;

                for(const auto& bkg : backgrounds){
                    const std::string hist_dir_name = channel+"_"+category+"/"+bkg;


                    auto hist_bkg = std::shared_ptr<TH1D>
                            (root_ext::TryReadObject<TH1D>(*inputFile,hist_dir_name));
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

                if (bkg_sum->GetNbinsX() != hist_data->GetNbinsX())
                    std::cout << "WARNING! bkg histogram has different binning from data histogram" <<  std::endl;
                for(int n = 1; n <= bkg_sum->GetNbinsX(); ++n){
//                    std::cout << "BKG - bin: " << n << ", content: " << bkg_sum->GetBinContent(n) << std::endl;
//                    std::cout << "DATA - bin: " << n << ", content: " << hist_data->GetBinContent(n) << std::endl;
                    if(bkg_sum->GetBinContent(n) == 0 && hist_data->GetBinContent(n) == 0)
                        std::cout << "BAD - both bkg_sum and data are empty - bin: " << n << std::endl;
                    if(bkg_sum->GetBinContent(n) == 0 && hist_data->GetBinContent(n) != 0)
                        std::cout << "CATASTROPHE - bkg_sum is empty and data NO - bin: " << n << ", data content: " <<
                                  hist_data->GetBinContent(n) << std::endl;
                }

            } //loop categories
        }//loop channels

    }





private:
    AnalyzerArguments args;


};

} // namespace analysis

PROGRAM_MAIN(analysis::BinSumBkg, analysis::AnalyzerArguments)
