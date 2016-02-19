/*!
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <string>
#include <vector>
#include <TString.h>
#include <TStyle.h>
#include <iostream>


using namespace std;



void AnalyzeBKGShape(TString folder_medium, TString folder_loose, TString sample, double max){
	gStyle->SetOptStat(000000);
	
	TH1::SetDefaultSumw2();
	
    //emb loose
    TFile *f_embedded= TFile::Open("/Users/Tita/Desktop/analysis_HH_bbTauTau/src/HHbbTauTau/data/datacards_tautau_ZTTembLoose_chi2/htt_tt.inputs-Hhh-8TeV_m_ttbb_kinfit_KinFitConvergedWithMassWindow.root");
    //MC medium
    //TFile *f_medium= TFile::Open("/Users/Tita/Desktop/analysis_HH_bbTauTau/src/HHbbTauTau/data/datacards_tautau_ZTT_MC_medium/htt_tt.inputs-Hhh-8TeV_m_ttbb_kinfit_KinFitConvergedWithMassWindow.root");
	

    //ZTT emb medium
    TH1D *histoZTT_medium = (TH1D*)f_embedded->Get((folder_medium+"/"+sample));
    //ZTT emb loose
    TH1D *histoZTT_loose = (TH1D*)f_embedded->Get((folder_loose+"/"+sample));
	
		
	TCanvas * c1 = new TCanvas("c1","c1", 800,800);

    histoZTT_medium->Scale(1/(histoZTT_medium->Integral()),"width");
    histoZTT_medium->SetLineColor(kBlue);
    histoZTT_medium->SetMarkerColor(kBlue);
    histoZTT_medium->SetMarkerStyle(20);
	
    histoZTT_loose->Scale(1/(histoZTT_loose->Integral()),"width");
    histoZTT_loose->SetLineColor(kRed);
    histoZTT_loose->SetMarkerColor(kRed);
    histoZTT_loose->SetMarkerStyle(20);

    histoZTT_medium->Draw("");
    histoZTT_loose->Draw("same");
    

    cout<<"--------------------------------------------------------------KS   "<<sample<<"--------------------------------------------------------------"<<endl;

//    TH1D* difference = (TH1D*)histoZTT_loose->Clone("difference");
//    difference->Add(histoZTT_medium,-1);
//    difference->Scale(1/difference->Integral());
    cout<<" Kolmogorov Test "<<histoZTT_medium->KolmogorovTest(histoZTT_loose,"")<<endl;
    histoZTT_medium->SetTitle("ZTT from Embedded");
    histoZTT_medium->GetXaxis()->SetTitle("M_{H} [GeV]");
    histoZTT_medium->GetYaxis()->SetTitleOffset(1.5);
    histoZTT_medium->GetYaxis()->SetTitle("N Events");
    histoZTT_medium->SetMaximum(max);
	
	
		
	
    TLegend* legend = new TLegend(0.65, 0.65, 0.99, 0.9);
	legend->SetFillColor(0);
    legend->SetTextSize(0.02);
	legend->SetEntrySeparation(0.05);
    legend->AddEntry(histoZTT_medium, " ZTT emb Loose 1tag");
    legend->AddEntry(histoZTT_loose, " ZTT emb Loose 2tag");
	
	legend->Draw();
	
    c1->SaveAs("./ZTT_Emb_Loose_1tag_vs_2tag_chi2cut.pdf");
}


void GetShape(){

    AnalyzeBKGShape("tauTau_2jetloose1tag", "tauTau_2jetloose2tag", "ZTT", 0.02);
    //AnalyzeBKGShape("tauTau_2jet2tag", "tauTau_2jetloose2tag", "ZTT", 0.1);
		
}
