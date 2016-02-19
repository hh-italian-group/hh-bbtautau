/*!
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <string>
#include <vector>
#include <TString.h>
#include <TStyle.h>
#include <iostream>


using namespace std;



void AnalyzeBKGShape(/*TString datacards_folder,*/ TString folder_0tag, TString folder_1tag, TString folder_2tag,
                     TString sample){
	gStyle->SetOptStat(000000);
	
	TH1::SetDefaultSumw2();
	
//    TFile *f= TFile::Open("/Users/Tita/Desktop/analysis_HH_bbTauTau/Limits/auxiliaries/shapes/Italians/htt_tt.inputs-Hhh-8TeV.root");
//    TFile *f= TFile::Open("/Users/Tita/GoogleDrive/HHbbTauTau/datacards_2015_01_12_tt_all/htt_tt.inputs-Hhh-8TeV_pt_1.root");
    TFile *f= TFile::Open("/Users/Tita/GoogleDrive/HHbbTauTau/datacards_2015_01_12_tt_all/htt_tt.inputs-Hhh-8TeV_pt_2.root");
	

	
    TH1D *histo_0tag = (TH1D*)f->Get((folder_0tag+"/"+sample));
    TH1D *histo_1tag = (TH1D*)f->Get((folder_1tag+"/"+sample));
    TH1D *histo_2tag = (TH1D*)f->Get((folder_2tag+"/"+sample));
	// TH1D *histoBLooseTauLoose = (TH1D*)fLooseBLooseIso->Get((folder+"/"+sample));
	
		
	TCanvas * c1 = new TCanvas("c1","c1", 800,800);

    histo_0tag->SetLineColor(kBlue);
    histo_0tag->SetMarkerColor(kBlue);
    histo_0tag->SetMarkerStyle(20);
	

    histo_1tag->SetLineColor(kRed);
    histo_1tag->SetMarkerColor(kRed);
    histo_1tag->SetMarkerStyle(20);

    histo_2tag->SetLineColor(kBlack);
    histo_2tag->SetMarkerColor(kBlack);
    histo_2tag->SetMarkerStyle(20);

	
	
	// histoBLooseTauLoose->SetLineColor(kGreen);
// 	histoBLooseTauLoose->SetMarkerColor(kGreen);
// 	histoBLooseTauLoose->SetMarkerStyle(20);
	
//    double error_histo = 0;
//    double integral_histo = histo->IntegralAndError(1, histo->GetNbinsX(),error_histo);
//    cout<<"--------------------------------------------------------------"<<sample<<"--------------------------------------------------------------"<<endl;
//    cout<<"*************"<<folder<<"  "<<sample<<"*************"<<endl;
//    cout<<" MC Shape "<<integral_histo<<"   ----  "<<error_histo<<endl;
//    cout<<" MC Shape Loose B "<<histoBLoose->Integral()<<endl;
//     cout<<" MC Shape Loose B + Relax Tau Iso "<<histoBLooseTauLoose->Integral()<<endl;
    
    
//     histoBLooseTauLoose->Scale(1/(histoBLooseTauLoose->Integral()));
//	TH1D* difference = (TH1D*)histoBLoose->Clone("difference");
//	difference->Add(histo,-1);
//	difference->Scale(1/difference->Integral());
    histo_0tag->Scale(1/(histo_0tag->Integral()));
    histo_1tag->Scale(1/(histo_1tag->Integral()));
    histo_2tag->Scale(1/(histo_2tag->Integral()));
	
	
    
    histo_0tag->Draw("");
    histo_1tag->Draw("same");
    histo_2tag->Draw("same");
    
    

    
    cout<<"--------------------------------------------------------------KS   "<<sample<<"  --------------------------------------------------------------"<<endl;
    cout<<"**************************"<<endl;
    cout<<" Kolmogorov Test 0tag-1tag: "<<histo_0tag->KolmogorovTest(histo_1tag,"")<<endl;
    cout<<" Kolmogorov Test 0tag-2tag: "<<histo_0tag->KolmogorovTest(histo_2tag,"")<<endl;
    cout<<" Kolmogorov Test 1tag-2tag: "<<histo_1tag->KolmogorovTest(histo_2tag,"")<<endl;
    cout<<"**************************"<<endl;
    histo_0tag->SetTitle(sample);
    histo_0tag->GetXaxis()->SetTitle("P_{T} [GeV]");
    histo_0tag->GetYaxis()->SetTitleOffset(1.5);
    histo_0tag->GetYaxis()->SetTitle("Normalized Events");
    histo_0tag->SetMaximum(0.4);
	
	
	
		
	
	TLegend* legend = new TLegend(0.6, 0.65, 0.99, 0.9);
	legend->SetFillColor(0);
	legend->SetTextSize(0.03);
	legend->SetEntrySeparation(0.05);
    legend->AddEntry(histo_0tag, "2jet0tag" );
    legend->AddEntry(histo_1tag, "2jet1tag" );
    legend->AddEntry(histo_2tag, "2jet2tag" );
	
	legend->Draw();
	
    c1->SaveAs("./pt_subleadingTau_"+sample+".pdf");
}


void GetShape(){


    AnalyzeBKGShape("tauTau_2jet0tag", "tauTau_2jet1tag", "tauTau_2jet2tag", "ZTT");
    AnalyzeBKGShape("tauTau_2jet0tag", "tauTau_2jet1tag", "tauTau_2jet2tag", "TT");
    AnalyzeBKGShape("tauTau_2jet0tag", "tauTau_2jet1tag", "tauTau_2jet2tag", "QCD");
    AnalyzeBKGShape("tauTau_2jet0tag", "tauTau_2jet1tag", "tauTau_2jet2tag", "data_obs");
    AnalyzeBKGShape("tauTau_2jet0tag", "tauTau_2jet1tag", "tauTau_2jet2tag", "ggHTohhTo2Tau2B300");
//    AnalyzeBKGShape("tauTau_2jet1tag", "ZTT");
//    AnalyzeBKGShape("tauTau_2jet2tag", "ZTT");
	
	
//	AnalyzeBKGShape("datacards_2014_12_18_tautau_Medium","tauTau_2jet1tag", "tauTau_2jetloose1tag", "QCD");
//	AnalyzeBKGShape("datacards_2014_12_18_tautau_Medium","tauTau_2jet2tag", "tauTau_2jetloose2tag", "QCD");
	
//	AnalyzeBKGShape("datacards_2014_12_18_tautau_Medium","tauTau_2jet1tag", "tauTau_2jetloose1tag", "VV");
//	AnalyzeBKGShape("datacards_2014_12_18_tautau_Medium","tauTau_2jet2tag", "tauTau_2jetloose2tag", "VV");
	
	

		
}
