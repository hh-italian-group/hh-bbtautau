/*! Class to calculate and apply PU reweighting.
Based on PhysicsTools/Utilities/src/LumiReWeighting.cc
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <TH1.h>

#include "h-tautau/Analysis/include/TreeExtractor.h"
#include "AnalysisTools/Core/include/RootExt.h"

class PileUp {
public:
    PileUp(const std::string& MC_File_name, const std::string& Data_File_name, const std::string& reweight_fileName,
                 const std::string& histName, const std::string& _mode, const std::string& _prefix = "none",
                        size_t _maxNumberOfEvents = 0)
        : treeExtractor(_prefix == "none" ? "" : _prefix, MC_File_name, false),
          outputFile(root_ext::CreateRootFile(reweight_fileName)),
          maxNumberOfEvents(_maxNumberOfEvents), mode(_mode)
    {

        auto Data_File = root_ext::OpenRootFile(Data_File_name);
        Data_distr = root_ext::ReadCloneObject(*Data_File, histName);
        root_ext::WriteObject(*Data_distr, outputFile.get());
        Data_distr->Scale( 1.0/ Data_distr->Integral() );
        Data_File->Close();
        nPU_MCdistr = new TH1D("MC_pileup", "MC nPU distribution", Data_distr->GetNbinsX(),
                               Data_distr->GetBinLowEdge(1),
                               Data_distr->GetBinLowEdge(Data_distr->GetNbinsX() + 1));
        nPU_MCdistr->SetDirectory(outputFile.get());
    }


    void Run()
    {
        size_t n = 0;
        analysis::EventDescriptor event;
        for(; !maxNumberOfEvents || n < maxNumberOfEvents; ++n) {
            if(!treeExtractor.ExtractNext(event))
                break;
            const ntuple::Event& eventInfo = event->eventInfo();
            for (unsigned n = 0; n < eventInfo.bunchCrossing.size(); ++n){
                if (eventInfo.bunchCrossing.at(n) == 0){
                    if(mode == "true")
                        nPU_MCdistr->Fill(eventInfo.trueNInt.at(n));
                    else if(mode == "observed")
                        nPU_MCdistr->Fill(eventInfo.nPU.at(n));
                    else
                        throw std::runtime_error("Unknown mode.");
                    break;
                }
            }
        }

        root_ext::WriteObject(*nPU_MCdistr);
        nPU_MCdistr->Scale( 1.0 / nPU_MCdistr->Integral() );

        TH1D* weights = root_ext::CloneObject(*Data_distr, "weights");
        weights->Divide(nPU_MCdistr);
        root_ext::WriteObject(*weights, outputFile.get());
    }

private:
    analysis::TreeExtractor treeExtractor;
    std::shared_ptr<TFile> outputFile;
    size_t maxNumberOfEvents;
    TH1D* Data_distr;
    TH1D* nPU_MCdistr;
    std::string mode;
};
