/*!
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <string>
#include <vector>
#include <TString.h>
#include <TStyle.h>
#include <iostream>


using namespace std;



void AnalyzeESshape(/*TString datacards_folder,*/ TString folder_0tag, TString folder_1tag, TString folder_2tag,
                     TString sample){
	gStyle->SetOptStat(000000);
	
	TH1::SetDefaultSumw2();
	
    TFile *f= TFile::Open("/Users/Tita/Desktop/analysis_HH_bbTauTau/Limits/auxiliaries/shapes/Italians/htt_tt.inputs-Hhh-8TeV.root");
	
    TH1D *histo_0tag_Up = (TH1D*)f->Get((folder_0tag+"/"+sample+"_CMS_scale_t_tautau_8TeVUp"));
    TH1D *histo_0tag_Central = (TH1D*)f->Get((folder_0tag+"/"+sample));
    TH1D *histo_0tag_Down = (TH1D*)f->Get((folder_0tag+"/"+sample+"_CMS_scale_t_tautau_8TeVDown"));
    TH1D *histo_1tag_Up = (TH1D*)f->Get((folder_1tag+"/"+sample+"_CMS_scale_t_tautau_8TeVUp"));
    TH1D *histo_1tag_Central = (TH1D*)f->Get((folder_1tag+"/"+sample));
    TH1D *histo_1tag_Down = (TH1D*)f->Get((folder_1tag+"/"+sample+"_CMS_scale_t_tautau_8TeVDown"));
    TH1D *histo_2tag_Up = (TH1D*)f->Get((folder_2tag+"/"+sample+"_CMS_scale_t_tautau_8TeVUp"));
    TH1D *histo_2tag_Central = (TH1D*)f->Get((folder_2tag+"/"+sample));
    TH1D *histo_2tag_Down = (TH1D*)f->Get((folder_2tag+"/"+sample+"_CMS_scale_t_tautau_8TeVDown"));

	
		
	TCanvas * c1 = new TCanvas("c1","c1", 800,800);

    TH1D* histo_0tag_Up_Central = (TH1D*)histo_0tag_Up->Clone("histo_0tag_Up_Central");
    histo_0tag_Up_Central->Divide(histo_0tag_Central);

    histo_0tag_Up_Central->SetLineColor(kBlue);
    histo_0tag_Up_Central->SetMarkerColor(kBlue);
    histo_0tag_Up_Central->SetMarkerStyle(20);

    TH1D* histo_1tag_Up_Central = (TH1D*)histo_1tag_Up->Clone("histo_1tag_Up_Central");
    histo_1tag_Up_Central->Divide(histo_1tag_Central);

    histo_1tag_Up_Central->SetLineColor(kRed);
    histo_1tag_Up_Central->SetMarkerColor(kRed);
    histo_1tag_Up_Central->SetMarkerStyle(20);

    TH1D* histo_2tag_Up_Central = (TH1D*)histo_2tag_Up->Clone("histo_2tag_Up_Central");
    histo_2tag_Up_Central->Divide(histo_2tag_Central);

    histo_2tag_Up_Central->SetLineColor(kBlack);
    histo_2tag_Up_Central->SetMarkerColor(kBlack);
    histo_2tag_Up_Central->SetMarkerStyle(20);
	
    
    
    histo_0tag_Up_Central->Draw("");
    histo_1tag_Up_Central->Draw("same");
    histo_2tag_Up_Central->Draw("same");
    
    

    
    cout<<"--------------------------------------------------------------KS   "<<sample<<"  --------------------------------------------------------------"<<endl;
    cout<<"**************************"<<endl;
    cout<<" Kolmogorov Test 0tag-1tag: "<<histo_0tag_Up_Central->KolmogorovTest(histo_1tag_Up_Central,"")<<endl;
    cout<<" Kolmogorov Test 0tag-2tag: "<<histo_0tag_Up_Central->KolmogorovTest(histo_2tag_Up_Central,"")<<endl;
    cout<<" Kolmogorov Test 1tag-2tag: "<<histo_1tag_Up_Central->KolmogorovTest(histo_2tag_Up_Central,"")<<endl;
    cout<<"**************************"<<endl;
    histo_0tag_Up_Central->SetTitle("Scale Up/ Central - "+sample);
    histo_0tag_Up_Central->GetXaxis()->SetTitle("M_{H} [GeV]");
    histo_0tag_Up_Central->GetYaxis()->SetTitleOffset(1.5);
    histo_0tag_Up_Central->GetYaxis()->SetTitle("Normalized Events");
    histo_0tag_Up_Central->SetMaximum(3.5);
		
	TLegend* legend = new TLegend(0.6, 0.65, 0.99, 0.9);
	legend->SetFillColor(0);
	legend->SetTextSize(0.03);
	legend->SetEntrySeparation(0.05);
    legend->AddEntry(histo_0tag_Up_Central, "2jet0tag" );
    legend->AddEntry(histo_1tag_Up_Central, "2jet1tag" );
    legend->AddEntry(histo_2tag_Up_Central, "2jet2tag" );
	
	legend->Draw();
	
    c1->SaveAs("./UP_CENTRAL_"+sample+".pdf");

    TCanvas * c2 = new TCanvas("c2","c2", 800,800);

    TH1D* histo_0tag_Down_Central = (TH1D*)histo_0tag_Down->Clone("histo_0tag_Down_Central");
    histo_0tag_Down_Central->Divide(histo_0tag_Central);

    histo_0tag_Down_Central->SetLineColor(kBlue);
    histo_0tag_Down_Central->SetMarkerColor(kBlue);
    histo_0tag_Down_Central->SetMarkerStyle(20);

    TH1D* histo_1tag_Down_Central = (TH1D*)histo_1tag_Down->Clone("histo_1tag_Down_Central");
    histo_1tag_Down_Central->Divide(histo_1tag_Central);

    histo_1tag_Down_Central->SetLineColor(kRed);
    histo_1tag_Down_Central->SetMarkerColor(kRed);
    histo_1tag_Down_Central->SetMarkerStyle(20);

    TH1D* histo_2tag_Down_Central = (TH1D*)histo_2tag_Down->Clone("histo_2tag_Down_Central");
    histo_2tag_Down_Central->Divide(histo_2tag_Central);

    histo_2tag_Down_Central->SetLineColor(kBlack);
    histo_2tag_Down_Central->SetMarkerColor(kBlack);
    histo_2tag_Down_Central->SetMarkerStyle(20);



    histo_0tag_Down_Central->Draw("");
    histo_1tag_Down_Central->Draw("same");
    histo_2tag_Down_Central->Draw("same");




    cout<<"--------------------------------------------------------------KS   "<<sample<<"  --------------------------------------------------------------"<<endl;
    cout<<"**************************"<<endl;
    cout<<" Kolmogorov Test 0tag-1tag: "<<histo_0tag_Down_Central->KolmogorovTest(histo_1tag_Down_Central,"")<<endl;
    cout<<" Kolmogorov Test 0tag-2tag: "<<histo_0tag_Down_Central->KolmogorovTest(histo_2tag_Down_Central,"")<<endl;
    cout<<" Kolmogorov Test 1tag-2tag: "<<histo_1tag_Down_Central->KolmogorovTest(histo_2tag_Down_Central,"")<<endl;
    cout<<"**************************"<<endl;
    histo_0tag_Down_Central->SetTitle("Scale Down/ Central - "+sample);
    histo_0tag_Down_Central->GetXaxis()->SetTitle("M_{H} [GeV]");
    histo_0tag_Down_Central->GetYaxis()->SetTitleOffset(1.5);
    histo_0tag_Down_Central->GetYaxis()->SetTitle("Normalized Events");
    histo_0tag_Down_Central->SetMaximum(2.5);

    TLegend* legend = new TLegend(0.6, 0.65, 0.99, 0.9);
    legend->SetFillColor(0);
    legend->SetTextSize(0.03);
    legend->SetEntrySeparation(0.05);
    legend->AddEntry(histo_0tag_Down_Central, "2jet0tag" );
    legend->AddEntry(histo_1tag_Down_Central, "2jet1tag" );
    legend->AddEntry(histo_2tag_Down_Central, "2jet2tag" );

    legend->Draw();

    c2->SaveAs("./DOWN_CENTRAL_"+sample+".pdf");
}


void GetShape(){


    AnalyzeESshape("tauTau_2jet0tag", "tauTau_2jet1tag", "tauTau_2jet2tag", "ZTT");
    AnalyzeESshape("tauTau_2jet0tag", "tauTau_2jet1tag", "tauTau_2jet2tag", "TT");
    AnalyzeESshape("tauTau_2jet0tag", "tauTau_2jet1tag", "tauTau_2jet2tag", "QCD");
    AnalyzeESshape("tauTau_2jet0tag", "tauTau_2jet1tag", "tauTau_2jet2tag", "ggHTohhTo2Tau2B300");
	
	

		
}
