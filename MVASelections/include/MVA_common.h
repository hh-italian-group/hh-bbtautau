/*!
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <TChain.h>
#include <TH1.h>
#include <TMVA/Factory.h>
#include <TMVA/Tools.h>
#include <TLorentzVector.h>

#include "MVA_variables.h"

typedef std::vector<size_t> IndexVector;

struct FlatTree {

    TChain* chain;

    std::vector<Float_t> *pt_BJets;
    std::vector<Float_t> *eta_BJets;
    std::vector<Float_t> *phi_BJets;
    std::vector<Float_t> *energy_BJets;
    std::vector<Float_t> *csv_BJets;

    Float_t mvamet_phi, mvamet, drLT, ptsv, etasv, phisv;
    Float_t pt1, pt2, eta1, eta2, phi1, phi2, e1, e2;
    Float_t byCombinedIsolationDeltaBetaCorrRaw3Hits1,byCombinedIsolationDeltaBetaCorrRaw3Hits2,pfRelIso1,mt1,mt2,weight;
    Int_t   Q1,Q2;
    Bool_t  againstMuonTight2,againstElectronLooseMVA2;


    FlatTree() : chain(new TChain("flatTree")), pt_BJets(0), eta_BJets(0), phi_BJets(0), energy_BJets(0), csv_BJets(0)
    {
        chain->SetBranchAddress( "pt_1", &pt1);
        chain->SetBranchAddress( "pt_2", &pt2);
        chain->SetBranchAddress( "eta_1", &eta1);
        chain->SetBranchAddress( "eta_2", &eta2);
        chain->SetBranchAddress( "phi_1", &phi1);
        chain->SetBranchAddress( "phi_2", &phi2);
        chain->SetBranchAddress( "energy_1", &e1);
        chain->SetBranchAddress( "energy_2", &e2);

        chain->SetBranchAddress( "pt_Bjets", &pt_BJets);
        chain->SetBranchAddress( "eta_Bjets", &eta_BJets);
        chain->SetBranchAddress( "phi_Bjets", &phi_BJets);
        chain->SetBranchAddress( "energy_Bjets", &energy_BJets);
        chain->SetBranchAddress( "csv_Bjets", &csv_BJets);
        chain->SetBranchAddress("mvametphi", &mvamet_phi);
        chain->SetBranchAddress("mvamet", &mvamet);
        chain->SetBranchAddress("DeltaR_leptons", &drLT);
        chain->SetBranchAddress("pt_sv", &ptsv);
        chain->SetBranchAddress("eta_sv", &etasv);
        chain->SetBranchAddress("phi_sv", &phisv);

        chain->SetBranchAddress("againstMuonTight_2",&againstMuonTight2);
        chain->SetBranchAddress("againstElectronLooseMVA_2",&againstElectronLooseMVA2);
        chain->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_1",&byCombinedIsolationDeltaBetaCorrRaw3Hits1);
        chain->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_2",&byCombinedIsolationDeltaBetaCorrRaw3Hits2);

        chain->SetBranchAddress("pfRelIso_1",&pfRelIso1);
        chain->SetBranchAddress("q_1",&Q1);
        chain->SetBranchAddress("q_2",&Q2);
        chain->SetBranchAddress("mt_1",&mt1);
        chain->SetBranchAddress("mt_2",&mt2);
        chain->SetBranchAddress("weight",&weight);
    }

    ~FlatTree() { delete chain; }

    Int_t Add(const std::string& name) { return chain->Add(name.c_str()); }
    Long64_t GetEntriesFast() const { return chain->GetEntriesFast(); }
    Int_t GetEntry(Long64_t entry) { return chain->GetEntry(entry); }
    std::string GetShortFileName()
    {
        const std::string full_file_name = chain->GetFile()->GetName();
        const size_t pos = full_file_name.find_last_of('/');
        if(pos == std::string::npos)
            return full_file_name;
        return full_file_name.substr(pos + 1);
    }

    IndexVector CollectJets() const
    {
        static const Double_t eta_cut = 2.4;
        IndexVector jet_indexes;
        for(size_t k = 0; k < eta_BJets->size(); ++k) {
            if(TMath::Abs(eta_BJets->at(k)) < eta_cut)
                jet_indexes.push_back(k);
        }
        return jet_indexes;
    }

    IndexVector CollectMediumBJets() const
    {
        static const Double_t medium_csv_wp = 0.679;
        return CollectBJets(medium_csv_wp);
    }

    IndexVector CollectTightBJets() const
    {
        static const Double_t tight_csv_wp = 0.898;
        return CollectBJets(tight_csv_wp);
    }

    IndexVector CollectBJets(Double_t csv_cut) const
    {
        static const Double_t eta_cut = 2.4;
        IndexVector jet_indexes;
        for(size_t k = 0; k < eta_BJets->size(); ++k) {
            if(TMath::Abs(eta_BJets->at(k)) < eta_cut && csv_BJets->at(k) > csv_cut)
                jet_indexes.push_back(k);
        }
        return jet_indexes;
    }

    TLorentzVector GetBJetMomentum(size_t bjet_index) const
    {
        TLorentzVector momentum;
        momentum.SetPtEtaPhiE(pt_BJets->at(bjet_index), eta_BJets->at(bjet_index),
                              phi_BJets->at(bjet_index), energy_BJets->at(bjet_index));
        return momentum;
    }

    TLorentzVector GetLeptonMomentum(size_t lepton_id) const
    {
        TLorentzVector momentum;
        if(lepton_id == 1)
            momentum.SetPtEtaPhiE(pt1, eta1, phi1, e1);
        else if(lepton_id == 2)
            momentum.SetPtEtaPhiE(pt2, eta2, phi2, e2);
        else
            throw std::runtime_error("Unknown lepton id");
        return momentum;
    }

    TLorentzVector GetMetMomentum() const
    {
        TLorentzVector momentum;
        momentum.SetPtEtaPhiE(mvamet, 0, mvamet_phi, mvamet);
        return momentum;
    }

    MVA_Selections::EventCategory DetermineEventCategory() const
    {
        const IndexVector jet_indexes = CollectJets();
        if(jet_indexes.size() == 1) {
            const IndexVector bjet_indexes = CollectTightBJets();
            if(bjet_indexes.size() == 0) return MVA_Selections::OneJet_ZeroBtag;
            if(bjet_indexes.size() == 1) return MVA_Selections::OneJet_OneBtag;
        }

        if(jet_indexes.size() >= 2) {
            const IndexVector bjet_indexes = CollectMediumBJets();
            if(bjet_indexes.size() == 0) return MVA_Selections::TwoJets_ZeroBtag;
            if(bjet_indexes.size() == 1) return MVA_Selections::TwoJets_OneBtag;
            if(bjet_indexes.size() >= 2) return MVA_Selections::TwoJets_TwoBtag;
        }

        return MVA_Selections::UnknownCategory;
    }

};

inline TH1F* MakeTH1F(const std::string& name, const std::string& prefix, Int_t nbinsx,Double_t xlow,Double_t xup)
{
    const std::string full_name = name + "_" + prefix;
    return new TH1F(full_name.c_str(), full_name.c_str(), nbinsx, xlow, xup);
}

class TH1F_Ptr {
public:
    TH1F_Ptr() : hist(0) {}
    TH1F_Ptr(TH1F* _hist) : hist(_hist) {}
    TH1F* operator-> () { return hist; }

    TH1F_Ptr& operator= (TH1F* _hist)
    {
        if(hist) {

            hist->Write();
            delete hist;
        }
        hist = _hist;
        return *this;
    }

    ~TH1F_Ptr()
    {
        if(hist) {
            hist->Write();
            delete hist;
        }
    }

private:
    TH1F_Ptr(const TH1F_Ptr& /*other*/) { }
    TH1F_Ptr& operator= (const TH1F_Ptr& /*other*/) { return *this; }

private:
    TH1F* hist;
};

struct MVA_Histograms {
    TH1F_Ptr pt_l1, pt_l2, pt_b1, pt_b2, dR_bb, dPhi_bb_MET, dR_ll,
             pt_Htt, dR_Hbb_Htt, pt_Hbb, dPhi_MET_tt, pt_H, mT2, mT1, pt_Htt_MET;

    MVA_Histograms(bool is_signal)
    {
        const std::string prefix = is_signal ? "signal" : "bkg";
        pt_l1 = MakeTH1F("hPt1",prefix,200,0,200);
        pt_l2 = MakeTH1F("hPt2",prefix,200,0,200);
        pt_b1 = MakeTH1F("hPtb1",prefix,300, 0, 300);
        pt_b2 = MakeTH1F("hPtb2",prefix,300, 0, 300);
        dR_bb = MakeTH1F("hDRbb",prefix,100, 0, 10);
        dPhi_bb_MET = MakeTH1F("hDPhiBBMET",prefix,100, -3, 3);
        dR_ll = MakeTH1F("hDRll",prefix,100, 0, 10);
        pt_Htt = MakeTH1F("hPtHtt",prefix,600, 0, 600);
        dR_Hbb_Htt= MakeTH1F("hDRHBBTT",prefix,100, 0, 10);
        pt_Hbb= MakeTH1F("hPtHBB",prefix,600, 0, 600);
        dPhi_MET_tt= MakeTH1F("hDeltaPhi_METTT",prefix,100, -3, 3);
        pt_H= MakeTH1F("hPtH",prefix,1000, 0, 1000);
        mT2= MakeTH1F("hmT2",prefix,600, 0, 600);
        mT1= MakeTH1F("hmT1",prefix,600, 0, 600);
        pt_Htt_MET= MakeTH1F("hPtHttMET",prefix,600, 0, 600);
    }
};

inline void AddEvent(TMVA::Factory *factory, const std::vector<Double_t>& vars, bool is_signal, Double_t weight)
{
    if (gRandom->Uniform(0, 1) < 0.5) {
        if(is_signal)
            factory->AddSignalTrainingEvent(vars, weight);
        else
            factory->AddBackgroundTrainingEvent(vars, weight);
    } else {
        if(is_signal)
            factory->AddSignalTestEvent(vars, weight);
        else
            factory->AddBackgroundTestEvent(vars, weight);
    }
}

inline Double_t GetFileScaleFactor(const std::string& file_name)
{
    static std::map<std::string, Double_t> file_name_map;
    if(!file_name_map.size()) {
        file_name_map["ggH_hh_bbtautau_260.root"] = 0.0358936;
        file_name_map["ggH_hh_bbtautau_270.root"] = 0.03057;
        file_name_map["ggH_hh_bbtautau_280.root"] = 0.0284073;
        file_name_map["ggH_hh_bbtautau_290.root"] = 0.0201117;
        file_name_map["ggH_hh_bbtautau_300.root"] = 0.0253141;
        file_name_map["ggH_hh_bbtautau_310.root"] = 0.0221424;
        file_name_map["ggH_hh_bbtautau_320.root"] = 0.0215317;
        file_name_map["ggH_hh_bbtautau_330.root"] = 0.0199818;
        file_name_map["ggH_hh_bbtautau_340.root"] = 0.0207792;
        file_name_map["ggH_hh_bbtautau_350.root"] = 0.0179311;

        file_name_map["tt_hadr.root"] = 0.0719394;
        file_name_map["tt_lept.root"] = 0.0429666;
        file_name_map["tt_semi.root"] = 0.0862387;
    }
    if(!file_name_map.count(file_name)) {
        std::cerr << "ERROR: unknown file name: " << file_name << std::endl;
        throw std::runtime_error("Unknown file name.");
    }
    return file_name_map[file_name];
}

bool ApplyFullEventSelection(const FlatTree* tree, std::vector<Double_t>& vars, MVA_Selections::EventCategory selectedCategory,
                             MVA_Selections::MvaMethod method, MVA_Histograms& h);

void ApplySelection(TMVA::Factory *factory, FlatTree* tree, MVA_Selections::EventCategory selectedCategory,
                    MVA_Selections::MvaMethod method, bool is_signal)
{
    MVA_Histograms h(is_signal);

    std::string file_name;
    double file_scale_factor;

    for (Long64_t i = 0; i < tree->GetEntriesFast(); i++) {
        tree->GetEntry(i);
        std::string current_file_name = tree->GetShortFileName();
        if(file_name != current_file_name) {
            file_scale_factor = GetFileScaleFactor(current_file_name);
            std::cout << "New file: " << current_file_name << ", SF = " << file_scale_factor << std::endl;
            file_name = current_file_name;
        }
        if(tree->DetermineEventCategory() != selectedCategory)
            continue;
        std::vector<Double_t> vars;
        const bool event_passed = ApplyFullEventSelection(tree, vars, selectedCategory, method, h);
        if(event_passed)
            AddEvent(factory, vars, is_signal, tree->weight * file_scale_factor);
    }
}

void TMVAtestAndTraining(TMVA::Factory *factory, MVA_Selections::MvaMethod method)
{
    factory->PrepareTrainingAndTestTree("", "","SplitMode=Random");

    if(method == MVA_Selections::BDT)
        factory->BookMethod( TMVA::Types::kBDT, "BDT",
                        "!H:!V:NTrees=850:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

    if(method == MVA_Selections::BDTMitFisher)
        factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
                            "!H:!V:NTrees=50:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

    if(method == MVA_Selections::BDTD)
        factory->BookMethod( TMVA::Types::kBDT, "BDTD","!H:!V:NTrees=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();
}

inline void AddVariablesToMVA(TMVA::Factory* factory, const MVA_Selections::str_vector& vars)
{
    for(size_t n = 0; n < vars.size(); ++n)
        factory->AddVariable(vars[n]);
}

inline void SetMvaInput(std::vector<Double_t>& vars, const MVA_Selections::var_map& var_names,
                        const std::string name, Double_t value)
{
    if(var_names.count(name))
       vars.at(var_names.find(name)->second) = value;
}
