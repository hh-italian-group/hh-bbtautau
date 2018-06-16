#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TTreeReader.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TROOT.h"
#include "TChain.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"

#include <vector>
#include <map>
#include <fstream>
void plot_LO_NLO(){
	/*TChain chain("muMu");
	chain.Add("/home/rbhattac/HH_BB_tau_tau_full_trees/DYJetsToLL_M-50_ext1.root");
	chain.Add("/home/rbhattac/HH_BB_tau_tau_full_trees/DYJetsToLL_M-50_ext2.root");
	TTreeReader reader_LO(&chain);
	*/
	TFile* f_LO = new TFile("Skimmed_muMu/DYJetsToLL_M-50.root");
	TTreeReader reader_LO("muMu",f_LO);

	TFile* f_NLO = new TFile("Skimmed_muMu/DYJetsToLL_M-50_NLO.root");
	TTreeReader reader_NLO("muMu",f_NLO);


	TTreeReaderArray<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>> genParticles_p4_LO(reader_LO,"genParticles_p4");
	TTreeReaderValue<float> genEventWeight_LO(reader_LO,"genEventWeight");
	TTreeReaderValue<unsigned int> lhe_n_partons_LO(reader_LO,"lhe_n_partons");
	TTreeReaderValue<unsigned int> lhe_n_b_partons_LO(reader_LO,"lhe_n_b_partons");
	TTreeReaderValue<unsigned int> lhe_n_c_partons_LO(reader_LO,"lhe_n_c_partons");
	TTreeReaderValue<float> lhe_HT_LO(reader_LO,"lhe_HT");
	TTreeReaderValue<double> weight_dy_LO(reader_LO,"weight_dy");


	TTreeReaderArray<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>> genParticles_p4_NLO(reader_NLO,"genParticles_p4");
	TTreeReaderValue<float> genEventWeight_NLO(reader_NLO,"genEventWeight");
	TTreeReaderValue<unsigned int> lhe_n_partons_NLO(reader_NLO,"lhe_n_partons");
	TTreeReaderValue<unsigned int> lhe_n_b_partons_NLO(reader_NLO,"lhe_n_b_partons");
	TTreeReaderValue<unsigned int> lhe_n_c_partons_NLO(reader_NLO,"lhe_n_c_partons");
	TTreeReaderValue<float> lhe_HT_NLO(reader_NLO,"lhe_HT");


	std::string categories[] = {"inclusive","0Jet","1Jet","1Jet_0bJet","1Jet_1bJet","2Jet","2Jet_0bJet","2Jet_1bJet","2Jet_2bJet"};

    

    std::map<std::string,TH1D*> pt_weight_histo_map;
    TFile* NLO_weight_file = new TFile("comp_LO_NLO_7.root");
    for(int i=0;i<sizeof(categories)/sizeof(std::string);i++){
    	std::string category = categories[i];
    	TH1D* h1;
    	NLO_weight_file->GetObject(("h_ratio_pt"+category).c_str(),h1);
    	pt_weight_histo_map[category] = h1;
    }

    std::map<std::string,double> fractional_weight_map;
    fractional_weight_map["0Jet"] = 0.93;
    fractional_weight_map["1Jet_0bJet"] = 1.02;
    fractional_weight_map["1Jet_1bJet"] = 1.38;
    fractional_weight_map["2Jet_0bJet"] = 0.99;
    fractional_weight_map["2Jet_1bJet"] = 1.15;
    fractional_weight_map["2Jet_2bJet"] = 1.39;

	std::map<std::string,TH1D*> LO_pt_histograms;
    std::map<std::string,TH1D*> NLO_pt_histograms;

    std::map<std::string,TH2D*> LO_pt_eta_histograms;
    std::map<std::string,TH2D*> NLO_pt_eta_histograms;

    std::map<std::string,TH1D*> LO_eta_histograms;
    std::map<std::string,TH1D*> NLO_eta_histograms;
    
    std::map<std::string,TH1D*> LO_n_partons_histograms;
    std::map<std::string,TH1D*> NLO_n_partons_histograms;

    std::map<std::string,TH1D*> LO_n_b_partons_histograms;
    std::map<std::string,TH1D*> NLO_n_b_partons_histograms;

    std::map<std::string,TH1D*> LO_n_c_partons_histograms;
    std::map<std::string,TH1D*> NLO_n_c_partons_histograms;

    std::map<std::string,TH1D*> LO_lhe_HT_histograms;
    std::map<std::string,TH1D*> NLO_lhe_HT_histograms;

    int bin_no = 27;
    double edges[28] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,350,500};

    int bin_no_0Jet = 9;
    double edges_0Jet[10] = {0,10,20,30,40,50,60,70,80,500};

    int bin_no_lhe_HT = 13;
    double edeges_lhe_HT[14] = {0,50,100,150,200,250,300,350,400,450,500,600,700,1000};


    for(int i=0;i<sizeof(categories)/sizeof(std::string);i++){
    	std::string category = categories[i];
    	if(category=="0Jet"){
			LO_pt_histograms[category] = new TH1D(("hist_LO_pt"+category).c_str(),("LO hist for pt "+category).c_str(),bin_no_0Jet,edges_0Jet);
			NLO_pt_histograms[category] = new TH1D(("hist_NLO_pt"+category).c_str(),("NLO hist for pt "+category).c_str(),bin_no_0Jet,edges_0Jet);

			LO_pt_eta_histograms[category] = new TH2D(("hist_2D_LO_pt_eta"+category).c_str(),("hist_2D_LO_pt_eta"+category).c_str(),bin_no_0Jet,edges_0Jet,60,-3,3);
			NLO_pt_eta_histograms[category] = new TH2D(("hist_2D_NLO_pt_eta"+category).c_str(),("hist_2D_NLO_pt_eta"+category).c_str(),bin_no_0Jet,edges_0Jet,60,-3,3);

    	}
    	else{
    		LO_pt_histograms[category] = new TH1D(("hist_LO_pt"+category).c_str(),("LO hist for pt "+category).c_str(),bin_no,edges);
    		NLO_pt_histograms[category] = new TH1D(("hist_NLO_pt"+category).c_str(),("NLO hist for pt "+category).c_str(),bin_no,edges);

    		LO_pt_eta_histograms[category] = new TH2D(("hist_2D_LO_pt_eta"+category).c_str(),("hist_2D_LO_pt_eta"+category).c_str(),bin_no,edges,60,-3,3);
			NLO_pt_eta_histograms[category] = new TH2D(("hist_2D_NLO_pt_eta"+category).c_str(),("hist_2D_NLO_pt_eta"+category).c_str(),bin_no,edges,60,-3,3);
    	}

    	LO_eta_histograms[category] = new TH1D(("hist_LO_eta"+category).c_str(),("LO hist for eta "+category).c_str(),60,-3,3);
    	NLO_eta_histograms[category] = new TH1D(("hist_NLO_eta"+category).c_str(),("NLO_hist_for_eta"+category).c_str(),60,-3,3);



    	LO_n_partons_histograms[category] = new TH1D(("hist_LO_n_partons"+category).c_str(),("LO hist for n_partons "+category).c_str(),6,-0.5,5.5);
    	NLO_n_partons_histograms[category] = new TH1D(("hist_NLO_n_partons"+category).c_str(),("NLO hist for n_partons "+category).c_str(),6,-0.5,5.5);
    	
    	LO_n_b_partons_histograms[category] = new TH1D(("hist_LO_n_b_partons"+category).c_str(),("LO hist for n_b_partons "+category).c_str(),6,-0.5,5.5);
    	NLO_n_b_partons_histograms[category] = new TH1D(("hist_NLO_n_b_partons"+category).c_str(),("NLO hist for n_b_partons "+category).c_str(),6,-0.5,5.5);
    	
    	LO_n_c_partons_histograms[category] = new TH1D(("hist_LO_n_c_partons"+category).c_str(),("LO hist for n_c_partons "+category).c_str(),6,-0.5,5.5);
    	NLO_n_c_partons_histograms[category] = new TH1D(("hist_NLO_n_c_partons"+category).c_str(),("NLO hist for n_c_partons "+category).c_str(),6,-0.5,5.5);

    	LO_lhe_HT_histograms[category] = new TH1D(("hist_LO_lhe_HT"+category).c_str(),("LO hist for lhe_HT "+category).c_str(),bin_no_lhe_HT,edeges_lhe_HT);
    	NLO_lhe_HT_histograms[category] = new TH1D(("hist_NLO_lhe_HT"+category).c_str(),("NLO hist for lhe_HT "+category).c_str(),bin_no_lhe_HT,edeges_lhe_HT);
    }
    TH1D* gen_particle_size_LO = new TH1D("gen_particle_size_LO","gen_particle_size_LO",6,-0.5,5.5);
    TH1D* gen_particle_sie_NLO = new TH1D("gen_particle_sie_NLO","gen_particle_sie_NLO",6,-0.5,5.5);

	while(reader_LO.Next()){
		float genEventWeight = *genEventWeight_LO;
		float lhe_n_partons = *lhe_n_partons_LO;
		float lhe_n_b_partons = *lhe_n_b_partons_LO;
		float lhe_n_c_partons = *lhe_n_c_partons_LO;
		double lhe_HT = *lhe_HT_LO;
		float weight_dy = *weight_dy_LO;
		std::string category = "";
		if(lhe_n_partons > 2) continue;
		
		if(lhe_n_partons==0) category = "0Jet";
		else if(lhe_n_partons==1){
			if(lhe_n_b_partons==0) category = "1Jet_0bJet";
			else if (lhe_n_b_partons==1) category = "1Jet_1bJet";
		}
		else if(lhe_n_partons == 2){
			if(lhe_n_b_partons==0) category = "2Jet_0bJet";
			else if (lhe_n_b_partons==1) category = "2Jet_1bJet";
			else if (lhe_n_b_partons==2) category = "2Jet_2bJet";
		}
		double pt_weight = 1;
		double fractional_weight = fractional_weight_map[category];
		gen_particle_size_LO->Fill(genParticles_p4_LO.GetSize());
		for(int i=0;i<genParticles_p4_LO.GetSize(); ++i){
			double pt;
			if(genParticles_p4_LO[i].Pt()>400) pt=450;
			else pt = genParticles_p4_LO[i].Pt();
			double eta = genParticles_p4_LO[i].Eta();
			pt_weight = pt_weight_histo_map[category]->GetBinContent(pt_weight_histo_map[category]->FindBin(pt));
			LO_pt_histograms["inclusive"]->Fill(pt,weight_dy*fractional_weight*pt_weight);
			LO_eta_histograms["inclusive"]->Fill(eta,weight_dy*fractional_weight*pt_weight);
			LO_pt_eta_histograms["inclusive"]->Fill(pt,eta,weight_dy*fractional_weight*pt_weight);
		
			LO_pt_histograms[category]->Fill(pt,weight_dy*fractional_weight*pt_weight);
			LO_eta_histograms[category]->Fill(eta,weight_dy*fractional_weight*pt_weight);
			LO_pt_eta_histograms[category]->Fill(pt,eta,weight_dy*fractional_weight*pt_weight);

			if(lhe_n_partons==1){
				LO_pt_histograms["1Jet"]->Fill(pt,weight_dy*fractional_weight*pt_weight);
				LO_eta_histograms["1Jet"]->Fill(eta,weight_dy*fractional_weight*pt_weight);
				LO_pt_eta_histograms["1Jet"]->Fill(pt,eta,weight_dy*fractional_weight*pt_weight);

			}
			if(lhe_n_partons==2){
			 	LO_pt_histograms["2Jet"]->Fill(pt,weight_dy*fractional_weight*pt_weight);
			 	LO_eta_histograms["2Jet"]->Fill(eta,weight_dy*fractional_weight*pt_weight);
				LO_pt_eta_histograms["2Jet"]->Fill(pt,eta,weight_dy*fractional_weight*pt_weight);
			}
		}
		LO_n_partons_histograms["inclusive"]->Fill(lhe_n_partons,weight_dy*fractional_weight*pt_weight);
		LO_n_b_partons_histograms["inclusive"]->Fill(lhe_n_b_partons,weight_dy*fractional_weight*pt_weight);
		LO_n_c_partons_histograms["inclusive"]->Fill(lhe_n_c_partons,weight_dy*fractional_weight*pt_weight);
		LO_lhe_HT_histograms["inclusive"]->Fill(lhe_HT,weight_dy*fractional_weight*pt_weight);

		LO_n_partons_histograms[category]->Fill(lhe_n_partons,weight_dy*fractional_weight*pt_weight);
		LO_n_b_partons_histograms[category]->Fill(lhe_n_b_partons,weight_dy*fractional_weight*pt_weight);
		LO_n_c_partons_histograms[category]->Fill(lhe_n_c_partons,weight_dy*fractional_weight*pt_weight);
		LO_lhe_HT_histograms[category]->Fill(lhe_HT,weight_dy*fractional_weight*pt_weight);

		if(lhe_n_partons==1){
			LO_n_partons_histograms["1Jet"]->Fill(lhe_n_partons,weight_dy*pt_weight*fractional_weight);
			LO_n_b_partons_histograms["1Jet"]->Fill(lhe_n_b_partons,weight_dy*pt_weight*fractional_weight);
			LO_n_c_partons_histograms["1Jet"]->Fill(lhe_n_c_partons,weight_dy*pt_weight*fractional_weight);
			LO_lhe_HT_histograms["1Jet"]->Fill(lhe_HT,weight_dy*pt_weight*fractional_weight);
		}
		
		if(lhe_n_partons==2){
			LO_n_partons_histograms["2Jet"]->Fill(lhe_n_partons,weight_dy*pt_weight*fractional_weight);
		 	LO_n_b_partons_histograms["2Jet"]->Fill(lhe_n_b_partons,weight_dy*pt_weight*fractional_weight);
	 		LO_n_c_partons_histograms["2Jet"]->Fill(lhe_n_c_partons,weight_dy*pt_weight*fractional_weight);
		 	LO_lhe_HT_histograms["2Jet"]->Fill(lhe_HT,weight_dy*pt_weight*fractional_weight);
		}
	}


	while(reader_NLO.Next()){
		float genEventWeight = *genEventWeight_NLO;
		float lhe_n_partons = *lhe_n_partons_NLO;
		float lhe_n_b_partons = *lhe_n_b_partons_NLO;
		float lhe_n_c_partons = *lhe_n_c_partons_NLO;
		float lhe_HT = *lhe_HT_NLO;
		if(lhe_n_partons > 2) continue;
		std::string category = "";
		if(lhe_n_partons==0) category = "0Jet";
		else if(lhe_n_partons==1){	
			if(lhe_n_b_partons==0) category = "1Jet_0bJet";
			else if (lhe_n_b_partons==1) category = "1Jet_1bJet";
		}
		else if(lhe_n_partons == 2){	
			if(lhe_n_b_partons==0) category = "2Jet_0bJet";
			else if (lhe_n_b_partons==1) category = "2Jet_1bJet";	
			else if (lhe_n_b_partons==2) category = "2Jet_2bJet";
		}
		gen_particle_sie_NLO->Fill(genParticles_p4_NLO.GetSize());
		for(int i=0;i<genParticles_p4_NLO.GetSize(); ++i){
			std::cout<<"B"<<std::endl;
			double pt;
			if(genParticles_p4_NLO[i].Pt()>400) pt=450;
			else pt = genParticles_p4_NLO[i].Pt();
			double eta = genParticles_p4_NLO[i].Eta();
			
		    NLO_pt_histograms["inclusive"]->Fill(pt,genEventWeight);
		    NLO_eta_histograms["inclusive"]->Fill(eta,genEventWeight);
		    NLO_pt_eta_histograms["inclusive"]->Fill(pt,eta,genEventWeight);
		
			NLO_pt_histograms[category]->Fill(pt,genEventWeight);
			NLO_eta_histograms[category]->Fill(eta,genEventWeight);
		    NLO_pt_eta_histograms[category]->Fill(pt,eta,genEventWeight);

			if(lhe_n_partons==1){
				NLO_pt_histograms["1Jet"]->Fill(pt,genEventWeight);
				NLO_eta_histograms["1Jet"]->Fill(eta,genEventWeight);
		    	NLO_pt_eta_histograms["1Jet"]->Fill(pt,eta,genEventWeight);
			}
			if(lhe_n_partons==2){
			 	NLO_pt_histograms["2Jet"]->Fill(pt,genEventWeight);
			 	NLO_eta_histograms["2Jet"]->Fill(eta,genEventWeight);
			 	NLO_pt_eta_histograms["2Jet"]->Fill(pt,eta,genEventWeight);
			}
		}
		NLO_n_partons_histograms["inclusive"]->Fill(lhe_n_partons,genEventWeight);
		NLO_n_b_partons_histograms["inclusive"]->Fill(lhe_n_b_partons,genEventWeight);
		NLO_n_c_partons_histograms["inclusive"]->Fill(lhe_n_c_partons,genEventWeight);
		NLO_lhe_HT_histograms["inclusive"]->Fill(lhe_HT,genEventWeight);

		NLO_n_partons_histograms[category]->Fill(lhe_n_partons,genEventWeight);
		NLO_n_b_partons_histograms[category]->Fill(lhe_n_b_partons,genEventWeight);
		NLO_n_c_partons_histograms[category]->Fill(lhe_n_c_partons,genEventWeight);
		NLO_lhe_HT_histograms[category]->Fill(lhe_HT,genEventWeight);

		if(lhe_n_partons==1){
			NLO_n_partons_histograms["1Jet"]->Fill(lhe_n_partons,genEventWeight);
			NLO_n_b_partons_histograms["1Jet"]->Fill(lhe_n_b_partons,genEventWeight);
			NLO_n_c_partons_histograms["1Jet"]->Fill(lhe_n_c_partons,genEventWeight);
			NLO_lhe_HT_histograms["1Jet"]->Fill(lhe_HT,genEventWeight);
		}
		if(lhe_n_partons==2){
			NLO_n_partons_histograms["2Jet"]->Fill(lhe_n_partons,genEventWeight);
		 	NLO_n_b_partons_histograms["2Jet"]->Fill(lhe_n_b_partons,genEventWeight);
	 		NLO_n_c_partons_histograms["2Jet"]->Fill(lhe_n_c_partons,genEventWeight);
		 	NLO_lhe_HT_histograms["2Jet"]->Fill(lhe_HT,genEventWeight);
		}

	}
	TFile* f_out = new TFile("comp_LO_NLO_8.root","RECREATE");
	f_out->cd();
	gen_particle_size_LO->Write();
	gen_particle_sie_NLO->Write();
	double inclusive_LO_integral = LO_pt_histograms["inclusive"]->Integral();
	double inclusive_NLO_integral = NLO_pt_histograms["inclusive"]->Integral();
	std::ofstream LO_output;
	LO_output.open("LO_output_2.txt");
	LO_output<<std::setw(10)<<" "<<std::setw(20)<<"Integral"<<std::setw(20)<<"Fraction"<<std::endl;


	std::ofstream NLO_output;
	NLO_output.open("NLO_output_2.txt");
	NLO_output<<std::setw(10)<<"Category"<<std::setw(20)<<"Integral"<<std::setw(20)<<"Fraction"<<std::endl;

  	for (int i = 0; i<sizeof(categories)/sizeof(std::string); ++i)
  	{
  		std::cout<<"C"<<std::endl;
		std::string category = categories[i];

		TH1D* hist_LO = LO_pt_histograms[category];
		TH1D* hist_NLO = NLO_pt_histograms[category];

		double fraction_LO = hist_LO->Integral()/inclusive_LO_integral;
		double fraction_NLO = hist_NLO->Integral()/inclusive_NLO_integral;

		LO_output<<std::setw(10)<<category<<std::setw(20)<<hist_LO->Integral()<<std::setw(20)<<fraction_LO<<std::endl;
		NLO_output<<std::setw(10)<<category<<std::setw(20)<<hist_NLO->Integral()<<std::setw(20)<<fraction_NLO<<std::endl;

		hist_LO->Scale(1./hist_LO->Integral());
		hist_NLO->Scale(1./hist_NLO->Integral());
		hist_LO->SetLineWidth(4);
    	hist_LO->SetLineColor(kRed);
    	hist_NLO->SetLineWidth(3);
    	//hist_NLO->SetMaximum(1.1*std::max(hist_NLO->GetBinContent(hist_NLO->GetMaximumBin()), hist_LO->GetBinContent(hist_LO->GetMaximumBin())));

		TH1D *h_ratio = (TH1D*)hist_NLO->Clone(("h_ratio_pt"+category).c_str());
		h_ratio->SetStats(0);
		h_ratio->Divide(hist_LO);


    	TCanvas *c = new TCanvas(("c_pt"+category).c_str(), ("canvas for"+ category).c_str(), 800, 800);
  	
  		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   		pad1->Draw();             // Draw the upper pad: pad1
   		pad1->cd();               // pad1 becomes the current pad
		if(hist_LO->GetBinContent(hist_LO->GetMaximumBin()) > hist_NLO->GetBinContent(hist_NLO->GetMaximumBin())){
			hist_LO->Draw("HIST");
			hist_NLO->Draw("HIST SAMES");
		}
		else{
			hist_NLO->Draw("HIST");
			hist_LO->Draw("HIST SAMES");
		}
   		TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   		axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   		axis->SetLabelSize(15);
   		axis->Draw();

 		c->cd();          // Go back to the main canvas before defining pad2
   		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   		pad2->Draw();
   		pad2->cd();       // pad2 becomes the current pad
   		h_ratio->Draw("ep");
   		hist_LO->Write();
   		hist_NLO->Write();
   		h_ratio->Write();
   		c->Write();
    }

    for (int i = 0; i<sizeof(categories)/sizeof(std::string); ++i)
  	{
  		std::cout<<"D"<<std::endl;
		std::string category = categories[i];

		TH1D* hist_LO = LO_n_partons_histograms[category];
		TH1D* hist_NLO = NLO_n_partons_histograms[category];

		hist_LO->Scale(1./hist_LO->Integral());
		hist_NLO->Scale(1./hist_NLO->Integral());
		hist_LO->SetLineWidth(4);
    	hist_LO->SetLineColor(kRed);
    	hist_NLO->SetLineWidth(3);
    	//hist_NLO->SetMaximum(1.1*std::max(hist_NLO->GetBinContent(hist_NLO->GetMaximumBin()), hist_LO->GetBinContent(hist_LO->GetMaximumBin())));

		TH1D *h_ratio = (TH1D*)hist_NLO->Clone(("h_ratio_n_partons"+category).c_str());
		h_ratio->SetStats(0);
		h_ratio->Divide(hist_LO);


    	TCanvas *c = new TCanvas(("c_n_partons"+category).c_str(), ("canvas for"+ category).c_str(), 800, 800);
  	
  		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   		pad1->Draw();             // Draw the upper pad: pad1
   		pad1->cd();               // pad1 becomes the current pad
		if(hist_LO->GetBinContent(hist_LO->GetMaximumBin()) > hist_NLO->GetBinContent(hist_NLO->GetMaximumBin())){
			hist_LO->Draw("HIST");
			hist_NLO->Draw("HIST SAMES");
		}
		else{
			hist_NLO->Draw("HIST");
			hist_LO->Draw("HIST SAMES");
		}
   		TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   		axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   		axis->SetLabelSize(15);
   		axis->Draw();

 		c->cd();          // Go back to the main canvas before defining pad2
   		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   		pad2->Draw();
   		pad2->cd();       // pad2 becomes the current pad
   		h_ratio->Draw("ep");
   		hist_LO->Write();
   		hist_NLO->Write();
   		h_ratio->Write();
   		c->Write();
    }

    for (int i = 0; i<sizeof(categories)/sizeof(std::string); ++i)
  	{
  		std::cout<<"E"<<std::endl;
		std::string category = categories[i];

		TH1D* hist_LO = LO_n_b_partons_histograms[category];
		TH1D* hist_NLO = NLO_n_b_partons_histograms[category];

		hist_LO->Scale(1./hist_LO->Integral());
		hist_NLO->Scale(1./hist_NLO->Integral());
		hist_LO->SetLineWidth(4);
		hist_LO->SetLineColor(kRed);
    	hist_NLO->SetLineWidth(3);
    	//hist_NLO->SetMaximum(1.1*std::max(hist_NLO->GetBinContent(hist_NLO->GetMaximumBin()), hist_LO->GetBinContent(hist_LO->GetMaximumBin())));

		TH1D *h_ratio = (TH1D*)hist_NLO->Clone(("h_ratio_n_b_partons"+category).c_str());
		h_ratio->SetStats(0);
		h_ratio->Divide(hist_LO);


    	TCanvas *c = new TCanvas(("c_n_b_partons"+category).c_str(), ("canvas for"+ category).c_str(), 800, 800);
  	
  		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   		pad1->Draw();             // Draw the upper pad: pad1
   		pad1->cd();               // pad1 becomes the current pad
		if(hist_LO->GetBinContent(hist_LO->GetMaximumBin()) > hist_NLO->GetBinContent(hist_NLO->GetMaximumBin())){
			hist_LO->Draw("HIST");
			hist_NLO->Draw("HIST SAMES");
		}
		else{
			hist_NLO->Draw("HIST");
			hist_LO->Draw("HIST SAMES");
		}
   		TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   		axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   		axis->SetLabelSize(15);
   		axis->Draw();

 		c->cd();          // Go back to the main canvas before defining pad2
   		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   		pad2->Draw();
   		pad2->cd();       // pad2 becomes the current pad
   		h_ratio->Draw("ep");
   		hist_LO->Write();
   		hist_NLO->Write();
   		h_ratio->Write();
   		c->Write();
    }

    for (int i = 0; i<sizeof(categories)/sizeof(std::string); ++i)
  	{
  		std::cout<<"F"<<std::endl;
		std::string category = categories[i];

		TH1D* hist_LO = LO_n_c_partons_histograms[category];
		TH1D* hist_NLO = NLO_n_c_partons_histograms[category];

		hist_LO->Scale(1./hist_LO->Integral());
		hist_NLO->Scale(1./hist_NLO->Integral());
		hist_LO->SetLineWidth(4);
    	hist_LO->SetLineColor(kRed);
    	hist_NLO->SetLineWidth(3);
    	//hist_NLO->SetMaximum(1.1*std::max(hist_NLO->GetBinContent(hist_NLO->GetMaximumBin()), hist_LO->GetBinContent(hist_LO->GetMaximumBin())));

		TH1D *h_ratio = (TH1D*)hist_NLO->Clone(("h_ratio_n_c_partons"+category).c_str());
		h_ratio->SetStats(0);
		h_ratio->Divide(hist_LO);


    	TCanvas *c = new TCanvas(("c_n_c_partons"+category).c_str(), ("canvas for"+ category).c_str(), 800, 800);
  	
  		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   		pad1->Draw();             // Draw the upper pad: pad1
   		pad1->cd();               // pad1 becomes the current pad
		if(hist_LO->GetBinContent(hist_LO->GetMaximumBin()) > hist_NLO->GetBinContent(hist_NLO->GetMaximumBin())){
			hist_LO->Draw("HIST");
			hist_NLO->Draw("HIST SAMES");
		}
		else{
			hist_NLO->Draw("HIST");
			hist_LO->Draw("HIST SAMES");
		}
   		TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   		axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   		axis->SetLabelSize(15);
   		axis->Draw();

 		c->cd();          // Go back to the main canvas before defining pad2
   		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   		pad2->Draw();
   		pad2->cd();       // pad2 becomes the current pad
   		h_ratio->Draw("ep");
   		hist_LO->Write();
   		hist_NLO->Write();
   		h_ratio->Write();
   		c->Write();
    }


    for (int i = 0; i<sizeof(categories)/sizeof(std::string); ++i)
  	{
  		std::cout<<"G"<<std::endl;
		std::string category = categories[i];

		TH1D* hist_LO = LO_lhe_HT_histograms[category];
		TH1D* hist_NLO = NLO_lhe_HT_histograms[category];

		hist_LO->Scale(1./hist_LO->Integral());
		hist_NLO->Scale(1./hist_NLO->Integral());
		hist_LO->SetLineWidth(4);
    	hist_LO->SetLineColor(kRed);
    	hist_NLO->SetLineWidth(3);
    	//hist_NLO->SetMaximum(1.1*std::max(hist_NLO->GetBinContent(hist_NLO->GetMaximumBin()), hist_LO->GetBinContent(hist_LO->GetMaximumBin())));

		TH1D *h_ratio = (TH1D*)hist_NLO->Clone(("h_ratio_lhe_HT"+category).c_str());
		h_ratio->SetStats(0);
		h_ratio->Divide(hist_LO);


    	TCanvas *c = new TCanvas(("c_lhe_HT"+category).c_str(), ("canvas for"+ category).c_str(), 800, 800);
  	
  		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   		pad1->Draw();             // Draw the upper pad: pad1
   		pad1->cd();               // pad1 becomes the current pad
		if(hist_LO->GetBinContent(hist_LO->GetMaximumBin()) > hist_NLO->GetBinContent(hist_NLO->GetMaximumBin())){
			hist_LO->Draw("HIST");
			hist_NLO->Draw("HIST SAMES");
		}
		else{
			hist_NLO->Draw("HIST");
			hist_LO->Draw("HIST SAMES");
		}
   		TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   		axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   		axis->SetLabelSize(15);
   		axis->Draw();

 		c->cd();          // Go back to the main canvas before defining pad2
   		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   		pad2->Draw();
   		pad2->cd();       // pad2 becomes the current pad
   		h_ratio->Draw("ep");
   		hist_LO->Write();
   		hist_NLO->Write();
   		h_ratio->Write();
   		c->Write();
    }


  for (int i = 0; i<sizeof(categories)/sizeof(std::string); ++i)
  	{
  		std::cout<<"G"<<std::endl;
		std::string category = categories[i];

		TH1D* hist_LO = LO_eta_histograms[category];
		TH1D* hist_NLO = NLO_eta_histograms[category];

		hist_LO->Scale(1./hist_LO->Integral());
		hist_NLO->Scale(1./hist_NLO->Integral());
		hist_LO->SetLineWidth(4);
    	hist_LO->SetLineColor(kRed);
    	hist_NLO->SetLineWidth(3);
    	//hist_NLO->SetMaximum(1.1*std::max(hist_NLO->GetBinContent(hist_NLO->GetMaximumBin()), hist_LO->GetBinContent(hist_LO->GetMaximumBin())));

		TH1D *h_ratio = (TH1D*)hist_NLO->Clone(("h_ratio_eta"+category).c_str());
		h_ratio->SetStats(0);
		h_ratio->Divide(hist_LO);


    	TCanvas *c = new TCanvas(("c_eta"+category).c_str(), ("canvas for"+ category).c_str(), 800, 800);
  	
  		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   		pad1->Draw();             // Draw the upper pad: pad1
   		pad1->cd();               // pad1 becomes the current pad
		if(hist_LO->GetBinContent(hist_LO->GetMaximumBin()) > hist_NLO->GetBinContent(hist_NLO->GetMaximumBin())){
			hist_LO->Draw("HIST");
			hist_NLO->Draw("HIST SAMES");
		}
		else{
			hist_NLO->Draw("HIST");
			hist_LO->Draw("HIST SAMES");
		}
   		TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   		axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   		axis->SetLabelSize(15);
   		axis->Draw();

 		c->cd();          // Go back to the main canvas before defining pad2
   		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   		pad2->Draw();
   		pad2->cd();       // pad2 becomes the current pad
   		h_ratio->Draw("ep");
   		hist_LO->Write();
   		hist_NLO->Write();
   		h_ratio->Write();
   		c->Write();
    }

      for (int i = 0; i<sizeof(categories)/sizeof(std::string); ++i)
  	{
  		std::cout<<"G"<<std::endl;
		std::string category = categories[i];

		TH2D* hist_LO = LO_pt_eta_histograms[category];
		TH2D* hist_NLO = NLO_pt_eta_histograms[category];

		hist_LO->Scale(1./hist_LO->Integral());
		hist_NLO->Scale(1./hist_NLO->Integral());
		/*hist_LO->SetLineWidth(4);
    	hist_LO->SetLineColor(kRed);
    	hist_NLO->SetLineWidth(3);*/
    	//hist_NLO->SetMaximum(1.1*std::max(hist_NLO->GetBinContent(hist_NLO->GetMaximumBin()), hist_LO->GetBinContent(hist_LO->GetMaximumBin())));

		TH2D *h_ratio = (TH2D*)hist_NLO->Clone(("h_ratio_pt_eta"+category).c_str());
		h_ratio->SetStats(0);
		h_ratio->Divide(hist_LO);


    	/*TCanvas *c = new TCanvas(("c_lhe_HT"+category).c_str(), ("canvas for"+ category).c_str(), 800, 800);
  	
  		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   		pad1->Draw();             // Draw the upper pad: pad1
   		pad1->cd();               // pad1 becomes the current pad
		if(hist_LO->GetBinContent(hist_LO->GetMaximumBin()) > hist_NLO->GetBinContent(hist_NLO->GetMaximumBin())){
			hist_LO->Draw("HIST");
			hist_NLO->Draw("HIST SAMES");
		}
		else{
			hist_NLO->Draw("HIST");
			hist_LO->Draw("HIST SAMES");
		}
   		TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   		axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   		axis->SetLabelSize(15);
   		axis->Draw();

 		c->cd();          // Go back to the main canvas before defining pad2
   		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   		pad2->Draw();
   		pad2->cd(); */      // pad2 becomes the current pad
   		h_ratio->Draw("ep");
   		hist_LO->Write();
   		hist_NLO->Write();
   		h_ratio->Write();
   		//c->Write();
    }



    f_out->Close();

    f_NLO->Close();

    LO_output.close();
    NLO_output.close();

    
}