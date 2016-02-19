/*!
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <string>
#include <vector>
#include <TString.h>
#include <TStyle.h>
#include <iostream>


using namespace std;



void AnalyzeBKGShape(TString folder, TString datacards_folder1, TString datacards_folder2, TString sample){
	gStyle->SetOptStat(000000);
	
	TH1::SetDefaultSumw2();
	
	TFile *f1= TFile::Open("/Users/Tita/Desktop/analysis_HH_bbTauTau/src/HHbbTauTau/data/"+datacards_folder1+"/htt_tt.inputs-Hhh-8TeV_m_sv.root");
	TFile *f2= TFile::Open("/Users/Tita/Desktop/analysis_HH_bbTauTau/src/HHbbTauTau/data/"+datacards_folder2+"/htt_tt.inputs-Hhh-8TeV_m_sv.root");
	

	
	TH1D *histoQCDplusW = (TH1D*)f1->Get((folder+"/"+sample));
	TH1D *histoQCDminusW = (TH1D*)f2->Get((folder+"/"+sample));
	
	TCanvas * c1 = new TCanvas("c1","c1", 800,800);

	histoQCDplusW->SetLineColor(kBlue);
	histoQCDplusW->SetMarkerColor(kBlue);
	histoQCDplusW->SetMarkerStyle(20);
	

	histoQCDminusW->SetLineColor(kRed);
	histoQCDminusW->SetMarkerColor(kRed);
	histoQCDminusW->SetMarkerStyle(20);

	
	
    histoQCDplusW->Draw("");
    histoQCDminusW->Draw("same");
	histoQCDplusW->SetTitle("M_{#tau#tau}");
	histoQCDplusW->GetXaxis()->SetTitle("M_{#tau#tau} [GeV]");
	histoQCDplusW->GetYaxis()->SetTitleOffset(1.5);
	histoQCDplusW->GetYaxis()->SetTitle("N Events");		
	
	TLegend* legend = new TLegend(0.7, 0.7, 0.99, 0.9);
	legend->SetFillColor(0);
	legend->SetTextSize(0.03);
	legend->SetEntrySeparation(0.05);
	legend->AddEntry(histoQCDplusW, " QCD + Wjets ");
	legend->AddEntry(histoQCDminusW, " QCD - Wjets ");
	
	legend->Draw();
	
	c1->SaveAs("./plots_QCD/QCD_Wjets_"+folder+".eps");
}


void GetShape(){

// compare QCD+Wjets and QCD - Wjets
	AnalyzeBKGShape("tauTau_2jet1tag", "datacards_TauTau_ZTT", "datacards_TauTau_WjetsFull", "QCD");
	AnalyzeBKGShape("tauTau_2jet2tag", "datacards_TauTau_ZTT", "datacards_TauTau_WjetsFull", "QCD");
	
	
		
}
