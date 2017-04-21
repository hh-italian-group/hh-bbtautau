/*! Check SM weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "Instruments/include/SMFileConfigEntryReader.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "Instruments/include/SampleDescriptor.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"


struct Arguments {
    run::Argument<std::string> weight_file{"weight_file", "file to reweight BSM samples"};
    run::Argument<std::string> SM_file{"SM_file", "SM file"};
    run::Argument<std::string> BSM_file{"BSM_file", "BSM file"};
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> output_file{"output_file", "Output root file"};
};

namespace analysis {

namespace sample_merging{

class SMCheckWeightData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    TH1D_ENTRY(lhe_hh_m, 100, 0, 1000)
    TH1D_ENTRY(lhe_hh_cosTheta, 300, -1.5, 1.5)
    //TH2D_ENTRY(lhe_hh_cosTheta_vs_m, 25, 200, 2000, 30, -1.1, 1.1)
    TH2D_ENTRY_CUSTOM(lhe_hh_cosTheta_vs_m, mhh_Bins(), cosTheta_Bins())
    TH2D_ENTRY_CUSTOM(weight, mhh_Bins(), cosTheta_Bins())
};


class SMCheckWeight {
public:
    SMCheckWeight(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.output_file())), anaData(output)
    {
        LoadInputs();
    }

    void Run() {}

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    SMCheckWeightData anaData;

    void LoadInputs()
    {
        auto inputFile_weight = root_ext::OpenRootFile(args.weight_file());
        auto inputFile_SM = root_ext::OpenRootFile(args.SM_file());
        auto inputFile_BSM = root_ext::OpenRootFile(args.BSM_file());

        TH2D *weight = (TH2D*)inputFile_weight->Get("weight_node_BSM");

        Int_t bin_x = weight->GetXaxis()->FindBin(350);
        Int_t bin_y = weight->GetYaxis()->FindBin(0.5);

        for (unsigned n = 0; n < weight->GetNbinsX(); ++n){
            for (unsigned h = 0; h < weight->GetNbinsY(); ++h){
                std::cout << " Bin Content: " << weight->GetBinContent(bin_x,bin_y) << std::endl ;
            }
        }

                std::shared_ptr<ntuple::ExpressTuple> AllEventTuple(new ntuple::ExpressTuple(args.tree_name(), inputFile.get(), true));
                const Long64_t n_entries = AllEventTuple->GetEntries();

                for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on entries

                    AllEventTuple->GetEntry(current_entry);
                    std::shared_ptr<ntuple::ExpressEvent> event(new ntuple::ExpressEvent(AllEventTuple->data()));

                    //std::cout << "Mhh: " << event->lhe_hh_m << ", cosTheta: " << event->lhe_hh_cosTheta << std::endl;
                    anaData.lhe_hh_m(name).Fill(event->lhe_hh_m);
                    anaData.lhe_hh_cosTheta(name).Fill(event->lhe_hh_cosTheta);
                    anaData.lhe_hh_cosTheta_vs_m(name).Fill(event->lhe_hh_m,event->lhe_hh_cosTheta);
                } //end loop on entries



        for (const auto& file_descriptor : file_descriptors) {
            const std::string& name = file_descriptor.second.name;

            for (unsigned n = 0; n <= anaData.lhe_hh_cosTheta_vs_m(name).GetNbinsX()+1; ++n){
                for (unsigned h = 0; h <= anaData.lhe_hh_cosTheta_vs_m(name).GetNbinsY()+1; ++h){
                    if (anaData.lhe_hh_cosTheta_vs_m(name).GetBinContent(n,h) == 0)
                        std::cout << name <<  " - Empty Bin! (" << n << "," << h << ")" << std::endl ;
                }
            }
            anaData.weight(name).CopyContent(anaData.lhe_hh_cosTheta_vs_m(name_sm));
            anaData.weight(name).Divide(&anaData.lhe_hh_cosTheta_vs_m(name));
        }
    }
};

} //namespace sample_merging

} //namespace analysis

PROGRAM_MAIN(analysis::sample_merging::SMCheckWeight, Arguments)
