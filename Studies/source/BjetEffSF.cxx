/*! Study of the b-jet efficiency scale factor properties.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <iostream>
#include <cmath>
#include <set>
#include <list>
#include <locale>

#include <TColor.h>
#include <TLorentzVector.h>

#include "AnalyzerData.h"
#include "h-tautau/Analysis/include/FlatEventInfo.h"
#include "AnalysisMath.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "exception.h"
#include "h-tautau/Analysis/include/Particles.h"
#include "AnalysisTools/Print/include/RootPrintToPdf.h"
#include "ProgressReporter.h"
#include "h-tautau/Analysis/include/Htautau_2015.h"
#include "h-tautau/Analysis/include/SyncTree.h"
#include "hh-bbtautau/Analysis/include/BTagCalibrationStandalone.h"


class BjetStudyData : public root_ext::AnalyzerData {
public:
    BjetStudyData(std::shared_ptr<TFile> outputFile) : root_ext::AnalyzerData(outputFile) {}

    TH1D_ENTRY(HWeight, 40, 0, 2)
};

struct jetEffInfo{
  float pt, eta, CSV;
  int   partonFlavour;
  float eff, SF;

  jetEffInfo(const ntuple::Sync& event, int jetIndex):eff(0.),SF(0.){
      pt  = event.pt_jets.at(jetIndex);
      eta = event.eta_jets.at(jetIndex);
      CSV = event.csv_jets.at(jetIndex);
      partonFlavour = event.partonFlavour_jets.at(jetIndex);
  }

};

typedef std::vector<jetEffInfo>    jetEffInfoVector;
typedef std::pair<bool,jetEffInfo> bTagOutcome;
typedef std::vector<bTagOutcome>   bTagMap;


class BjetEffSF {

public:
  BjetEffSF(const std::string& inputFileName, const std::string& outputFileName,
            const std::string& bTagEffName, const std::string& bjetSFName)
  :inputFile(root_ext::OpenRootFile(inputFileName)),
   outputFile(root_ext::CreateRootFile(outputFileName)),
   bTagEffFile(root_ext::OpenRootFile(bTagEffName)),
   syncTree(new ntuple::SyncTree("sync", inputFile.get(), true)),anaData(outputFile){

     calib        = new btag_calibration::BTagCalibration("CSVv2", bjetSFName);
     reader       = new btag_calibration::BTagCalibrationReader(calib, btag_calibration::BTagEntry::OP_LOOSE, "mujets",
                                                                "central");
     reader_light = new btag_calibration::BTagCalibrationReader(calib, btag_calibration::BTagEntry::OP_LOOSE, "incl",
                                                                "central");

     eff_b = (TH2F*) bTagEffFile->Get("eff_b");
     eff_c = (TH2F*) bTagEffFile->Get("eff_c");
     eff_l = (TH2F*) bTagEffFile->Get("eff_l");

     progressReporter = std::shared_ptr<analysis::tools::ProgressReporter>(
                 new analysis::tools::ProgressReporter(10, std::cout));
   }

  ~BjetEffSF() {}


  void Run(){
    for(Long64_t current_entry = 0; current_entry < syncTree->GetEntries(); ++current_entry) {
        //eventInfoMap.clear();
        using namespace cuts::Htautau_2015;
        bTagMap bTagMedium;
        progressReporter->Report(current_entry);
        syncTree->GetEntry(current_entry);
        const ntuple::Sync& event = syncTree->data();
        //std::cout<<" SyncTree Entry :   "<<current_entry<<"\n";
        for (int jetIndex = 0; jetIndex < event.njets; jetIndex++){
          jetEffInfo jetInfo(event, jetIndex);
          if      (abs (jetInfo.partonFlavour) == 5 ) {
            jetInfo.SF  = reader->eval(btag_calibration::BTagEntry::FLAV_B, jetInfo.eta, jetInfo.pt);
            jetInfo.eff = GetEfficiency(eff_b,jetInfo.pt,TMath::Abs(jetInfo.eta));
          }
          else if (abs (jetInfo.partonFlavour) == 4 ) {
            jetInfo.SF = reader->eval(btag_calibration::BTagEntry::FLAV_C, jetInfo.eta, jetInfo.pt);
            jetInfo.eff = GetEfficiency(eff_c,jetInfo.pt,TMath::Abs(jetInfo.eta));
          }
          else {
            jetInfo.SF = reader_light->eval(btag_calibration::BTagEntry::FLAV_UDSG, jetInfo.eta, jetInfo.pt);
            jetInfo.eff = GetEfficiency(eff_l,jetInfo.pt,TMath::Abs(jetInfo.eta));
          }

          bTagOutcome jetOutcome(false,jetInfo);
          if ( jetOutcome.second.CSV > btag::CSVL ) jetOutcome.first = true;
          bTagMedium.push_back(jetOutcome);

          // std::cout<<"\t  jet number "<<jetIndex<<" Flavour : "<<jetOutcome.second.partonFlavour
          //                                       <<" ---->   pt = "<<jetOutcome.second.pt
          //                                              <<" eta = "<<jetOutcome.second.eta
          //                                              <<" CSV = "<<jetOutcome.second.CSV
          //                                              <<" eff = "<<jetOutcome.second.eff
          //                                              <<"  SF = "<<jetOutcome.second.SF<<std::endl;
        }
        double weight = GetBtagWeight(bTagMedium);
        anaData.HWeight().Fill(weight);
        //std::cout<<" Event btag Weight  --->  "<< weight <<std::endl;

  }

  progressReporter->Report(syncTree->GetEntries(),true);
}

private:
  double GetEfficiency (const TH2F* h, const double pt, const double eta){
    int xBin = h->GetXaxis()->FindFixBin(pt);
    int yBin = h->GetYaxis()->FindFixBin(eta);

    return h->GetBinContent(xBin,yBin);
  }

  double GetBtagWeight (const bTagMap map){
    double MC=1;
    double Data=1;
    //std::cout<<"  -------------   Weight Computation ------------ "<<std::endl;
    for (auto element : map){
      // std::cout<<"\t\t\t Flavour : "<<element.second.partonFlavour
      //                                       <<" ---->   pt = "<<element.second.pt
      //                                              <<" eta = "<<element.second.eta
      //                                              <<" CSV = "<<element.second.CSV
      //                                              <<" eff = "<<element.second.eff
      //                                              <<"  SF = "<<element.second.SF<<std::endl;

      if(element.first){
        MC=MC*element.second.eff;
        Data=Data*(element.second.eff*element.second.SF);
      }
      if(!(element.first)){
        MC=MC*(1-element.second.eff);
        Data=Data*(1-(element.second.eff*element.second.SF));
      }
    }
    //std::cout<<"   Data Weight = "<< Data <<"   MC Weight  = "<<MC<<std::endl;

    if (!MC) return 0;
    return Data/MC;
  }


private:
  std::shared_ptr<TFile> inputFile, outputFile, bTagEffFile;
  std::shared_ptr<ntuple::SyncTree> syncTree;
  BjetStudyData anaData;
  btag_calibration::BTagCalibration *calib;
  btag_calibration::BTagCalibrationReader *reader, *reader_light;

  TH2F* eff_b = 0;
  TH2F* eff_c = 0;
  TH2F* eff_l = 0;

  std::shared_ptr<analysis::tools::ProgressReporter> progressReporter;
};
