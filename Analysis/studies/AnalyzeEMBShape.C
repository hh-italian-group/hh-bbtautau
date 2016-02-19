/*!
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <string>
#include <vector>
#include <TString.h>
#include <TStyle.h>
#include <iostream>


using namespace std;



void AnalyzeBKGShape(TString folder, TString sample1, TString sample2){
	gStyle->SetOptStat(000000);
	
	TH1::SetDefaultSumw2();
	
	TFile *f= TFile::Open("/Users/Tita/Desktop/analysis_HH_bbTauTau/src/HHbbTauTau/data/datacards_TauTau_ZTT/htt_tt.inputs-Hhh-8TeV_m_ttbb_kinfit_KinFitConvergedWithMassWindow.root");
	

	
	TH1D *histoDY = (TH1D*)f->Get((folder+"/"+sample1));
	TH1D *histoTT = (TH1D*)f->Get((folder+"/"+sample2));
	
		
	TCanvas * c1 = new TCanvas("c1","c1", 800,800);

	histoDY->SetLineColor(kBlue);
	histoDY->SetMarkerColor(kBlue);
	histoDY->SetMarkerStyle(20);
	

	histoTT->SetLineColor(kRed);
	histoTT->SetMarkerColor(kRed);
	histoTT->SetMarkerStyle(20);

	
	
	
	
    // double error_histo = 0;
//     double integral_histo = histo->IntegralAndError(1, histo->GetNbinsX(),error_histo);
//     cout<<"--------------------------------------------------------------"<<sample<<"--------------------------------------------------------------"<<endl;
//     cout<<"*************"<<folder<<"  "<<sample<<"*************"<<endl;
//     cout<<" MC Shape "<<integral_histo<<"   ----  "<<error_histo<<endl;
//     cout<<" MC Shape Loose B "<<histoBLoose->Integral()<<endl;
//     cout<<" MC Shape Loose B + Relax Tau Iso "<<histoBLooseTauLoose->Integral()<<endl;
    
    
//     histoBLooseTauLoose->Scale(1/(histoBLooseTauLoose->Integral()));
    // histoTT->Scale(1/(histoTT->Integral()));
//     histoDY->Scale(1/(histoDY->Integral()));
	
	cout << "Integral DY = " << histoDY->Integral() << ", Integral TT = " << histoTT->Integral() << endl;
	cout << "percentage = " << histoTT->Integral()/histoDY->Integral() << endl;
    
    histoDY->Draw("");
    histoTT->Draw("same");
    
//     histoBLooseTauLoose->Draw("same");
    
    

    
    // cout<<"--------------------------------------------------------------KS   "<<sample<<"--------------------------------------------------------------"<<endl;
//     cout<<"*************"<<folder<<"  "<<sample<<"*************"<<endl;
//     cout<<" Kolmogorov Test "<<histoDY->KolmogorovTest(histoBLoose,"")<<endl;
	histoDY->SetTitle("M_{H}");
	histoDY->GetXaxis()->SetTitle("M_{H} [GeV]");
	histoDY->GetYaxis()->SetTitleOffset(1.5);
	histoDY->GetYaxis()->SetTitle("N Events");
	
	
	
		
	
	TLegend* legend = new TLegend(0.6, 0.65, 0.99, 0.9);
	legend->SetFillColor(0);
	legend->SetTextSize(0.05);
	legend->SetEntrySeparation(0.05);
	legend->AddEntry(histoDY, " DY embedded ");
	legend->AddEntry(histoTT, " TT embedded ");
	//legend->AddEntry(histoBLooseTauLoose, " No Btagging, TauIso<1 ");
	
	legend->Draw();
	
	c1->SaveAs("./plots_embedded/"+folder+"_embComparison.eps");
}


void GetShape(){

	AnalyzeBKGShape("tauTau_2jet1tag", "DY_emb", "TT_emb");
	AnalyzeBKGShape("tauTau_2jet2tag", "DY_emb", "TT_emb");
		
}
