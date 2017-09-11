#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TChain.h>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include <iterator>
#include <cmath>
#include <TLorentzVector.h>
#include "Math/GenVector/LorentzVector.h"
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TROOT.h"
#include "TStopwatch.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

#define PI 3.14159265358979323846

// hh-bbtautau framework
#include "h-tautau/Analysis/include/EventInfo.h"


/*struct bjet
{
	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>> jet;
	double csv;
	bool operator > (const bjet& a) const
	{
		return csv > a.csv;
	}
};*/

namespace analysis {

// class to evalluate the BDT score
class BDT_reader
{
	using EventInfo = analysis::EventInfo<ElectronCandidate>;

	public:
		BDT_reader(TString file);
		~BDT_reader();
		double BDT_score(EventInfo& event);
	
	private:
		// Declare MVAreader, weight file, MVAmethodname
		TMVA::Reader *reader = new TMVA::Reader( "!Color:Silent" );
		TString weight_file;
		//TString weight_file = "hh-bbtautau/Analysis/data/TMVAClassification_BDT_full_mass.weights.xml";
		TString methodname = "BDT_full_mass";
	
		// Declare variables for MVAreader
		Float_t var1, var2, var3, var4, var5, var6, var7, var8, var9;
		Float_t var10, var11, var12, var13, var14, var15, var16, var17, var18, var19;
		Float_t var20, var21, var22, var23, var24, var25, var26;
	
		// Fill variables methods
		Float_t get_njets20(EventInfo& event);
		Float_t get_njets30(EventInfo& event);
		Float_t get_njets50(EventInfo& event);
		Float_t get_HT(EventInfo& event);
		Float_t get_HT50(EventInfo& event);
		Float_t get_HTfull(EventInfo& event);
		Float_t get_centrality(EventInfo& event);
		Float_t get_dphi_mumet(EventInfo& event);
		Float_t get_dphi_taumu(EventInfo& event);
		Float_t get_dphi_metsv(EventInfo& event);
		Float_t get_dR_bb(EventInfo& event);
		Float_t get_dR_taumu(EventInfo& event);
		Float_t get_dR_bbsv(EventInfo& event);
		Float_t get_jet_max_eta(EventInfo& event);
		Float_t get_lead_jet_pt(EventInfo& event);
		Float_t get_lead_jet_eta(EventInfo& event);
		Float_t get_sublead_jet_pt(EventInfo& event);
		Float_t get_dphi_bbmet(EventInfo& event);
		Float_t get_dphi_bbsv(EventInfo& event);
	
		// Other methods
		Float_t deltaPhi(TLorentzVector& v1, TLorentzVector& v2);
		Float_t deltaEta(TLorentzVector& v1, TLorentzVector& v2);
		Float_t deltaR(TLorentzVector& v1, TLorentzVector& v2);
		std::vector<bjet> find_bjets(EventInfo& event);
};


// ***************** Member functions definitions *****************
// Constructor
BDT_reader::BDT_reader(TString file)
{
	std::cout<<" --- Class BDT_reader is being created ---"<<std::endl;
	std::cout<<file<<std::endl;
	
	weight_file = file;
	
	//reader->AddVariable("HT",&var1);
	//reader->AddVariable("HT50",&var2);
	reader->AddVariable("HTfull",&var3);
	//reader->AddVariable("centrality",&var4);
	reader->AddVariable("dphi_mumet",&var5);
	//reader->AddVariable("dphi_taumu",&var6);
	reader->AddVariable("dphi_metsv",&var7);
	reader->AddVariable("dR_bb",&var8);
	reader->AddVariable("dR_bbsv",&var9);
	reader->AddVariable("dR_taumu",&var10);
	//reader->AddVariable("jet_max_eta",&var11);
	//reader->AddVariable("lead_jet_eta",&var12);
	//reader->AddVariable("lead_jet_pt",&var13);
	//reader->AddVariable("p4_1.Eta()",&var14);
	//reader->AddVariable("njets20",&var15);
	//reader->AddVariable("njets30",&var16);
	//reader->AddVariable("njets50",&var17);
	//reader->AddVariable("pfMET_p4.Pt()",&var18);
	reader->AddVariable("pfmt_1",&var19);
	reader->AddVariable("pfmt_2",&var20);
	//reader->AddVariable("sublead_jet_pt",&var21);
	//reader->AddVariable("p4_2.Pt()",&var22);
	//reader->AddVariable("p4_2.Eta()",&var23);
	reader->AddVariable("dphi_bbmet",&var24);
	reader->AddVariable("dphi_bbsv",&var25);
	reader->AddVariable("channel",&var26);
	
	// --- Book the MVA methods
	reader->BookMVA(methodname, weight_file);
}

// Destructor
BDT_reader::~BDT_reader(void)
{
	std::cout<<" --- Class BDT_reader is being destroyed ---"<<std::endl;
	delete reader;
}

// Calculate the BDT value
double BDT_reader::BDT_score(EventInfo& event)
{
	//std::cout<<" --- EVENT ---\n";

	double BDT_value=0.0;
	
	var1 = get_HT(event);
	var2 = get_HT50(event);
	var3 = get_HTfull(event);
	var4 = get_centrality(event);
	var5 = get_dphi_mumet(event);
	var6 = get_dphi_taumu(event);
	var7 = get_dphi_metsv(event);
	var8 = get_dR_bb(event);
	var9 = get_dR_bbsv(event);
	var10 = get_dR_taumu(event);
	var11 = get_jet_max_eta(event);
	var12 = get_lead_jet_eta(event);
	var13 = get_lead_jet_pt(event);
	var14 = event->p4_1.Eta();
	var15 = get_njets20(event);
	var16 = get_njets30(event);
	var17 = get_njets50(event);
	var18 = event->pfMET_p4.Pt();
	var19 = event->pfmt_1;
	var20 = event->pfmt_2;
	var21 = get_sublead_jet_pt(event);
	var22 = event->p4_2.Pt();
	var23 = event->p4_2.Eta();
	var24 = get_dphi_bbmet(event);
	var25 = get_dphi_bbsv(event);
	var26 = 1;
	
	BDT_value = reader->EvaluateMVA(methodname);
	
	return BDT_value;
}

// delta phi
Float_t BDT_reader::deltaPhi(TLorentzVector& v1, TLorentzVector& v2)
{
	Float_t res = v1.Phi() - v2.Phi();
	if (res > PI)  res = res - 2*PI;
	if (res < -PI) res = res + 2*PI;
	return std::abs(res);
}

// deltaETA
Float_t BDT_reader::deltaEta(TLorentzVector& v1, TLorentzVector& v2)
{
	Float_t res = v1.Eta() - v2.Eta();
	return std::abs(res);
}

// deltaR
Float_t BDT_reader::deltaR(TLorentzVector& v1, TLorentzVector& v2)
{
	Float_t dphi=deltaPhi(v1,v2);
	Float_t deta=deltaEta(v1,v2);
	return std::sqrt((dphi*dphi) + (deta*deta));
}

// find bjets
std::vector<bjet> BDT_reader::find_bjets(EventInfo& event)
{
	std::vector<bjet> bjets;
	for (unsigned int i = 0; i< (event->jets_p4).size(); i++)
	{
		if (event->jets_p4[i].Pt()>20. && std::abs((event->jets_p4[i].Eta())) < 2.4 && event->jets_csv[i]>0.65)
		{
			bjet temp_bjet;
			temp_bjet.jet = event->jets_p4[i];
			temp_bjet.csv = event->jets_csv[i];
			bjets.push_back(temp_bjet);
			//std::cout<<" - bjet(jet,csv,i):"<< temp_bjet.jet << " " << temp_bjet.csv << " " << i <<std::endl;
		}
	}
	return bjets;
}

// HT
Float_t BDT_reader::get_HT(EventInfo& event)
{
	Float_t HT = 0;
	for (unsigned int i=0; i< (event->jets_p4).size(); i++)
	{
		if (event->jets_p4[i].Pt() > 20.)
		{
			HT = HT + event->jets_p4[i].Pt();
		}
	}
	return HT;
}

// HT50
Float_t BDT_reader::get_HT50(EventInfo& event)
{
	Float_t HT50 = 0;
	for (unsigned int i = 0; i< (event->jets_p4).size(); i++)
	{
		if (event->jets_p4[i].Pt() > 50.)
		{
			HT50 = HT50 + event->jets_p4[i].Pt();
		}
	}
	return HT50;
}


// HTfull
Float_t BDT_reader::get_HTfull(EventInfo& event)
{
	Float_t HTfull = 0;
	for (unsigned int i = 0; i< (event->jets_p4).size(); i++)
	{
		if (event->jets_p4[i].Pt() > 20.)
		{
			HTfull = HTfull + event->jets_p4[i].Pt();
		}
	}
	HTfull = HTfull + (event->p4_1.Pt()) + (event->p4_2.Pt());
	return HTfull;
}

// centrality
Float_t BDT_reader::get_centrality(EventInfo& event)
{
	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>> tot_jet;
	for (unsigned int i=0; i< (event->jets_p4).size(); i++)
	{
		tot_jet = tot_jet + event->jets_p4[i];
	}
	return tot_jet.Pt() / get_HT(event);
}


// dphi_mumet
Float_t BDT_reader::get_dphi_mumet(EventInfo& event)
{
	TLorentzVector mu;
	TLorentzVector met;
	mu.SetPtEtaPhiM(event->p4_1.Pt(), event->p4_1.Eta(), event->p4_1.Phi(), event->p4_1.M());
	met.SetPtEtaPhiM(event->pfMET_p4.Pt(), event->pfMET_p4.Eta(), event->pfMET_p4.Phi(), event->pfMET_p4.M());

	return deltaPhi(mu, met);
}

// dphi_taumu
Float_t BDT_reader::get_dphi_taumu(EventInfo& event)
{
	TLorentzVector mu;
	TLorentzVector tau;
	mu.SetPtEtaPhiM(event->p4_1.Pt(), event->p4_1.Eta(), event->p4_1.Phi(), event->p4_1.M());
	tau.SetPtEtaPhiM(event->p4_2.Pt(), event->p4_2.Eta(), event->p4_2.Phi(), event->p4_2.M());
	
	return deltaPhi(mu, tau);
}

// dphi_metsv
Float_t BDT_reader::get_dphi_metsv(EventInfo& event)
{
	TLorentzVector sv;
	TLorentzVector met;
	sv.SetPtEtaPhiM(event->SVfit_p4.Pt(), event->SVfit_p4.Eta(), event->SVfit_p4.Phi(), event->SVfit_p4.M());
	met.SetPtEtaPhiM(event->pfMET_p4.Pt(), event->pfMET_p4.Eta(), event->pfMET_p4.Phi(), event->pfMET_p4.M());
	
	return deltaPhi(sv, met);
}

// dphi_bbmet
Float_t BDT_reader::get_dphi_bbmet(EventInfo& event)
{
	Float_t res = 0.;
	TLorentzVector met;
	met.SetPtEtaPhiM(event->pfMET_p4.Pt(), event->pfMET_p4.Eta(), event->pfMET_p4.Phi(), event->pfMET_p4.M());
	std::vector<bjet> bjets = find_bjets(event);
	
	if (bjets.size() >= 2)
	{
		std::sort(bjets.begin(), bjets.end(), std::greater<bjet>());
		TLorentzVector b1, b2, bsum;
		b1.SetPtEtaPhiM(bjets[0].jet.Pt(), bjets[0].jet.Eta(), bjets[0].jet.Phi(), bjets[0].jet.M());
		b2.SetPtEtaPhiM(bjets[1].jet.Pt(), bjets[1].jet.Eta(), bjets[1].jet.Phi(), bjets[1].jet.M());
		bsum = b1+b2;
		res = deltaPhi(bsum, met);
	}
	else
	{
		res =  -99;
	}
	return res;
}

// dphi_bbsv
Float_t BDT_reader::get_dphi_bbsv(EventInfo& event)
{
	Float_t res = 0.;
	TLorentzVector sv;
	sv.SetPtEtaPhiM(event->SVfit_p4.Pt(), event->SVfit_p4.Eta(), event->SVfit_p4.Phi(), event->SVfit_p4.M());
	
	std::vector<bjet> bjets = find_bjets(event);
	
	if (bjets.size() >= 2)
	{
		std::sort(bjets.begin(), bjets.end(), std::greater<bjet>());
		TLorentzVector b1, b2, bsum;
		b1.SetPtEtaPhiM(bjets[0].jet.Pt(), bjets[0].jet.Eta(), bjets[0].jet.Phi(), bjets[0].jet.M());
		b2.SetPtEtaPhiM(bjets[1].jet.Pt(), bjets[1].jet.Eta(), bjets[1].jet.Phi(), bjets[1].jet.M());
		bsum = b1+b2;
		res = deltaPhi(bsum, sv);
	}
	else
	{
		res =  -99;
	}
	return res;
}

// dR_bb
Float_t BDT_reader::get_dR_bb(EventInfo& event)
{
	Float_t res = 0.;
	std::vector<bjet> bjets = find_bjets(event);
	/*for (unsigned int i=0; i < bjets.size(); i++)
	{
		std::cout<<"   bjet: "<<bjets[i].jet<<" "<<bjets[i].csv<<std::endl;
	}*/
	if (bjets.size() >= 2)
	{
		std::sort(bjets.begin(), bjets.end(), std::greater<bjet>());
		/*for (unsigned int i=0; i < bjets.size(); i++)
		 {
		 std::cout<<"   order bjet: "<<bjets[i].jet<<" "<<bjets[i].csv<<std::endl;
		 }*/
		TLorentzVector b1, b2;
		b1.SetPtEtaPhiM(bjets[0].jet.Pt(), bjets[0].jet.Eta(), bjets[0].jet.Phi(), bjets[0].jet.M());
		b2.SetPtEtaPhiM(bjets[1].jet.Pt(), bjets[1].jet.Eta(), bjets[1].jet.Phi(), bjets[1].jet.M());
		res = deltaR(b1,b2);
	}
	else
	{
		res = -99.;
	}
	return res;
}

// dR_taumu
Float_t BDT_reader::get_dR_taumu(EventInfo& event)
{
	TLorentzVector mu;
	TLorentzVector tau;
	mu.SetPtEtaPhiM(event->p4_1.Pt(), event->p4_1.Eta(), event->p4_1.Phi(), event->p4_1.M());
	tau.SetPtEtaPhiM(event->p4_2.Pt(), event->p4_2.Eta(), event->p4_2.Phi(), event->p4_2.M());
	
	return deltaR(mu, tau);
}

// dR_bbsv
Float_t BDT_reader::get_dR_bbsv(EventInfo& event)
{
	Float_t res = 0.;
	TLorentzVector sv;
	sv.SetPtEtaPhiM(event->SVfit_p4.Pt(), event->SVfit_p4.Eta(), event->SVfit_p4.Phi(), event->SVfit_p4.M());

	std::vector<bjet> bjets = find_bjets(event);
	
	if (bjets.size() >= 2)
	{
		std::sort(bjets.begin(), bjets.end(), std::greater<bjet>());
		TLorentzVector b1, b2, bsum;
		b1.SetPtEtaPhiM(bjets[0].jet.Pt(), bjets[0].jet.Eta(), bjets[0].jet.Phi(), bjets[0].jet.M());
		b2.SetPtEtaPhiM(bjets[1].jet.Pt(), bjets[1].jet.Eta(), bjets[1].jet.Phi(), bjets[1].jet.M());
		bsum = b1+b2;
		res = deltaR(bsum, sv);
	}
	else
	{
		res =  -99;
	}
	return res;
}

// jet max eta
Float_t BDT_reader::get_jet_max_eta(EventInfo& event)
{
	Float_t jet_max_eta = 0;
	for (unsigned int i = 0; i< (event->jets_p4).size(); i++)
	{
		if (std::abs(event->jets_p4[i].Eta()) > std::abs(jet_max_eta))
		{
			jet_max_eta = std::abs(event->jets_p4[i].Eta());
		}
	}
	return jet_max_eta;

}

// njets20
Float_t BDT_reader::get_njets20(EventInfo& event)
{
	Float_t njets20 = 0;
	for (unsigned int i = 0; i< (event->jets_p4).size(); i++)
	{
		if (event->jets_p4[i].Pt() > 20.)
		{
			njets20 = njets20 +1;
		}
	}
	return njets20;
}

// njets30
Float_t BDT_reader::get_njets30(EventInfo& event)
{
	Float_t njets30 = 0;
	for (unsigned int i = 0; i< (event->jets_p4).size(); i++)
	{
		if (event->jets_p4[i].Pt() > 30.)
		{
			njets30 = njets30 +1;
		}
	}
	return njets30;
}

// njets50
Float_t BDT_reader::get_njets50(EventInfo& event)
{
	Float_t njets50 = 0;
	for (unsigned int i = 0; i< (event->jets_p4).size(); i++)
	{
		if (event->jets_p4[i].Pt() > 50.)
		{
			njets50 = njets50 +1;
		}
	}
	return njets50;
}

// lead jet pt
Float_t BDT_reader::get_lead_jet_pt(EventInfo& event)
{
	Float_t lead_jet_pt = 0;
	for (unsigned int i = 0; i< (event->jets_p4).size(); i++)
	{
		if (event->jets_p4[i].Pt() > lead_jet_pt)
		{
			lead_jet_pt = event->jets_p4[i].Pt();
		}
	}
	return lead_jet_pt;
}

// lead jet eta
Float_t BDT_reader::get_lead_jet_eta(EventInfo& event)
{
	Float_t lead_jet_eta = 0;
	Float_t lead_jet_pt = 0;
	for (unsigned int i = 0; i< (event->jets_p4).size(); i++)
	{
		if (event->jets_p4[i].Pt() > lead_jet_pt)
		{
			lead_jet_eta = event->jets_p4[i].Eta();
		}
	}
	return lead_jet_eta;
}

// sublead_jet_pt
Float_t BDT_reader::get_sublead_jet_pt(EventInfo& event)
{
	Float_t sublead_jet_pt = 0;
	sublead_jet_pt = event->jets_p4[1].Pt();
	return sublead_jet_pt;
}

} //namespace analysis
