/*!
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

#include "../include/MVA_common.h"

bool ApplyFullEventSelection(const FlatTree* tree, std::vector<Double_t>& vars, MVA_Selections::EventCategory selectedCategory,
                             MVA_Selections::MvaMethod method, MVA_Histograms& h)
{
    if(!tree->againstMuonTight2 || !tree->againstElectronLooseMVA2 ||
            tree->byCombinedIsolationDeltaBetaCorrRaw3Hits2 >= 1.5
            || tree->pfRelIso1 >= 0.1
            || tree->Q1 * tree->Q2 != -1 || tree->mt1 >= 30)
        return false;

    if(selectedCategory != MVA_Selections::TwoJets_TwoBtag && selectedCategory != MVA_Selections::TwoJets_OneBtag &&
            selectedCategory != MVA_Selections::TwoJets_ZeroBtag)
        throw std::runtime_error("Unsupported event category.");

    const IndexVector bjet_indexes = tree->CollectJets(); //to run categories 2jets_*tag
    if(bjet_indexes.size() < 2)
        throw std::runtime_error("Unconsistent category information.");

    const TLorentzVector b1 = tree->GetBJetMomentum(bjet_indexes.at(0));
    const TLorentzVector b2 = tree->GetBJetMomentum(bjet_indexes.at(1));
    const TLorentzVector BB = b1 + b2;
    const TLorentzVector l1 = tree->GetLeptonMomentum(1);
    const TLorentzVector l2 = tree->GetLeptonMomentum(2);
    const TLorentzVector TT = l1 + l2;
    const TLorentzVector H = BB + TT;
    const TLorentzVector MET = tree->GetMetMomentum();
    const TLorentzVector TT_MET = TT + MET;

    const MVA_Selections::ParamId key(MVA_Selections::MuTau, selectedCategory, method);
    const MVA_Selections::var_map& var_names = MVA_Selections::Input_Variables_Map(key);
    vars.assign(var_names.size(), 0);
    SetMvaInput(vars, var_names, "pt_mu", l1.Pt());
    SetMvaInput(vars, var_names, "pt_tau", l2.Pt());
    SetMvaInput(vars, var_names, "pt_b1", b1.Pt());
    SetMvaInput(vars, var_names, "pt_b2", b2.Pt());
    SetMvaInput(vars, var_names, "DR_bb", b1.DeltaR(b2));
    SetMvaInput(vars, var_names, "DPhi_BBMET", MET.DeltaPhi(BB));
    SetMvaInput(vars, var_names, "DR_ll", l1.DeltaR(l2));
    SetMvaInput(vars, var_names, "Pt_Htt", TT.Pt());
    SetMvaInput(vars, var_names, "DR_HBBHTT", TT.DeltaR(BB));
    SetMvaInput(vars, var_names, "Pt_Hbb", BB.Pt());
    SetMvaInput(vars, var_names, "DeltaPhi_METTT", MET.DeltaPhi(TT));
    SetMvaInput(vars, var_names, "PtH", H.Pt());
    SetMvaInput(vars, var_names, "mT2", tree->mt2);
    SetMvaInput(vars, var_names, "mT1", tree->mt1);
    SetMvaInput(vars, var_names, "Pt_Htt_MET", TT_MET.Pt());

    h.pt_l1->Fill(tree->pt1);
    h.pt_l2->Fill(tree->pt2);
    h.pt_b1->Fill(b1.Pt());
    h.pt_b2->Fill(b2.Pt());
    h.dR_bb->Fill(b1.DeltaR(b2));
    h.dPhi_bb_MET->Fill(MET.DeltaPhi(BB));
    h.dR_ll->Fill(tree->drLT);
    h.pt_Htt->Fill(TT.Pt());
    h.dR_Hbb_Htt->Fill(TT.DeltaR(BB));
    h.pt_Hbb->Fill(BB.Pt());
    h.dPhi_MET_tt->Fill(MET.DeltaPhi(TT));
    h.pt_H->Fill(H.Pt());
    h.mT2->Fill(tree->mt2);
    h.mT1->Fill(tree->mt1);
    h.pt_Htt_MET->Fill(TT_MET.Pt());

    return true;
}

void MVA_mutau(const std::string& filePath, const std::string& category_name, const std::string& method_name,
               const std::string& suffix)
{
    std::cout << "==> Start TMVAClassification" << std::endl;

    MVA_Selections::EventCategory eventCategory = MVA_Selections::EventCategoryFromString(category_name);
    MVA_Selections::MvaMethod method = MVA_Selections::MvaMethodFromString(method_name);
    std::cout << "Training for event category = '" << category_name << "', mva method '" << method_name << "'" << std::endl;

    const std::string output_name = "TMVA_mutau_" + category_name + "_" + method_name + "_" + suffix;
    TFile* outputFile = TFile::Open( (output_name + ".root").c_str(), "RECREATE" );

    TMVA::Factory *factory = new TMVA::Factory( output_name, outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

    const MVA_Selections::ParamId key(MVA_Selections::MuTau, eventCategory, method);
    AddVariablesToMVA(factory, MVA_Selections::Input_Variables(key));

    FlatTree* sigTree = new FlatTree();
    FlatTree* bkgTree = new FlatTree();

    sigTree->Add(filePath+"ggH_hh_bbtautau_*.root");
    bkgTree->Add(filePath+"tt_*.root");

    outputFile->cd();

    ApplySelection(factory, sigTree, eventCategory, method, true);
    ApplySelection(factory, bkgTree, eventCategory, method, false);

    TMVAtestAndTraining(factory, method);

    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;

    delete factory;
    delete sigTree;
    delete bkgTree;
}


