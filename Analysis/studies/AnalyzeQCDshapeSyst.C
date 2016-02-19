/*!
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <string>
#include <vector>
#include <TString.h>
#include <TStyle.h>
#include <iostream>


using namespace std;


//comparing ratio
//void AnalyzeBKGShape(TString folder_0tag, TString folder_1tag, TString folder_2tag, TString sample1, TString sample2){
//	gStyle->SetOptStat(000000);
	
//	TH1::SetDefaultSumw2();
	
//    TFile *f= TFile::Open("/Users/Tita/Desktop/analysis_HH_bbTauTau/src/HHbbTauTau/data/datacards_tautau_try_2/htt_tt.inputs-Hhh-8TeV.root");
	

	
//    TH1D *histo_0tag_QCD = (TH1D*)f->Get((folder_0tag+"/"+sample1));
//    TH1D *histo_0tag_QCDalt = (TH1D*)f->Get((folder_0tag+"/"+sample2));

//    TH1D *histo_1tag_QCD = (TH1D*)f->Get((folder_1tag+"/"+sample1));
//    TH1D *histo_1tag_QCDalt = (TH1D*)f->Get((folder_1tag+"/"+sample2));

//    TH1D *histo_2tag_QCD = (TH1D*)f->Get((folder_2tag+"/"+sample1));
//    TH1D *histo_2tag_QCDalt = (TH1D*)f->Get((folder_2tag+"/"+sample2));




//    TCanvas * c1 = new TCanvas("c1","c1", 800,800);

//    TH1D* histo_0tag_Difference = (TH1D*)histo_0tag_QCD->Clone("histo_0tag_Difference");
//    histo_0tag_Difference->Divide(histo_0tag_QCDalt);
//    //histo_0tag_Difference->Scale(1/(histo_0tag_Difference->Integral()));

//    histo_0tag_Difference->SetLineColor(kBlue);
//    histo_0tag_Difference->SetMarkerColor(kBlue);
//    histo_0tag_Difference->SetMarkerStyle(20);

//    TH1D* histo_1tag_Difference = (TH1D*)histo_1tag_QCD->Clone("histo_1tag_Difference");
//    histo_1tag_Difference->Divide(histo_1tag_QCDalt);
//    //histo_1tag_Difference->Scale(1/(histo_1tag_Difference->Integral()));

//    histo_1tag_Difference->SetLineColor(kRed);
//    histo_1tag_Difference->SetMarkerColor(kRed);
//    histo_1tag_Difference->SetMarkerStyle(20);

//    TH1D* histo_2tag_Difference = (TH1D*)histo_2tag_QCD->Clone("histo_2tag_Difference");
//    histo_2tag_Difference->Divide(histo_2tag_QCDalt);
//    //histo_2tag_Difference->Scale(1/(histo_2tag_Difference->Integral()));

//    histo_2tag_Difference->SetLineColor(kBlack);
//    histo_2tag_Difference->SetMarkerColor(kBlack);
//    histo_2tag_Difference->SetMarkerStyle(20);


//    histo_0tag_Difference->Draw("");
//    histo_1tag_Difference->Draw("same");
//    histo_2tag_Difference->Draw("same");




//    cout<<"--------------------------------------------------------------KS--------------------------------------------------------------"<<endl;
//    cout<<"**************************"<<endl;
//    cout<<" Kolmogorov Test 0tag-1tag: "<<histo_0tag_Difference->KolmogorovTest(histo_1tag_Difference,"")<<endl;
//    cout<<" Kolmogorov Test 0tag-2tag: "<<histo_0tag_Difference->KolmogorovTest(histo_2tag_Difference,"")<<endl;
//    cout<<" Kolmogorov Test 1tag-2tag: "<<histo_1tag_Difference->KolmogorovTest(histo_2tag_Difference,"")<<endl;
//    cout<<"**************************"<<endl;
//    histo_0tag_Difference->SetTitle("QCD shape correlation between categories ");
//    histo_0tag_Difference->GetXaxis()->SetTitle("M_{H} [GeV]");
//    histo_0tag_Difference->GetYaxis()->SetTitleOffset(1.5);
//    histo_0tag_Difference->GetYaxis()->SetTitle("Normalized Events");
//    histo_0tag_Difference->SetMaximum(8);

//    TLegend* legend = new TLegend(0.6, 0.65, 0.99, 0.9);
//    legend->SetFillColor(0);
//    legend->SetTextSize(0.03);
//    legend->SetEntrySeparation(0.05);
//    legend->AddEntry(histo_0tag_Difference, "2jet0tag" );
//    legend->AddEntry(histo_1tag_Difference, "2jet1tag" );
//    legend->AddEntry(histo_2tag_Difference, "2jet2tag" );

//    legend->Draw();

//    c1->SaveAs("./QCD_shape_correlation.pdf");

//}


//void GetShape(){

//    AnalyzeBKGShape("tauTau_2jet0tag", "tauTau_2jet1tag", "tauTau_2jet2tag", "QCD", "QCD_alternative");
		
//}

void AnalyzeBKGShape(TString folder, TString sample1, TString sample2){
    gStyle->SetOptStat(000000);

    TH1::SetDefaultSumw2();

    TFile *f= TFile::Open("/Users/Tita/Desktop/analysis_HH_bbTauTau/src/HHbbTauTau/data/datacards_tautau_try_2/htt_tt.inputs-Hhh-8TeV.root");

    TH1D *histo_QCD = (TH1D*)f->Get((folder+"/"+sample1));
    TH1D *histo_QCDalt = (TH1D*)f->Get((folder+"/"+sample2));

    TCanvas * c1 = new TCanvas("c1","c1", 800,800);


    histo_QCD->Scale(1/(histo_QCD->Integral()),"width");

    histo_QCD->SetLineColor(kBlue);
    histo_QCD->SetMarkerColor(kBlue);
    histo_QCD->SetMarkerStyle(20);

    histo_QCDalt->Scale(1/(histo_QCDalt->Integral()),"width");

    histo_QCDalt->SetLineColor(kRed);
    histo_QCDalt->SetMarkerColor(kRed);
    histo_QCDalt->SetMarkerStyle(20);

    histo_QCD->Draw("");
    histo_QCDalt->Draw("same");

    cout<<"--------------------------------------------------------------KS--------------------------------------------------------------"<<endl;
    cout<<"**************************"<<endl;
    cout<<" Kolmogorov Test: "<<histo_QCD->KolmogorovTest(histo_QCDalt,"")<<endl;

    cout<<"**************************"<<endl;
    histo_QCD->SetTitle("QCD shape ");
    histo_QCD->GetXaxis()->SetTitle("M_{H} [GeV]");
    histo_QCD->GetYaxis()->SetTitleOffset(1.5);
    histo_QCD->GetYaxis()->SetTitle("Normalized Events");
//    histo_0tag_Difference->SetMaximum(8);

    TLegend* legend = new TLegend(0.6, 0.65, 0.99, 0.9);
    legend->SetFillColor(0);
    legend->SetTextSize(0.03);
    legend->SetEntrySeparation(0.05);
    legend->AddEntry(histo_QCD, "QCD" );
    legend->AddEntry(histo_QCDalt, "QCD alt" );

    legend->Draw();

    c1->SaveAs("./QCD_shape_"+folder+".pdf");

}


void GetShape(){

    AnalyzeBKGShape("tauTau_2jet0tag", "QCD", "QCD_alternative");
    AnalyzeBKGShape("tauTau_2jet1tag", "QCD", "QCD_alternative");
    AnalyzeBKGShape("tauTau_2jet2tag", "QCD", "QCD_alternative");

}
