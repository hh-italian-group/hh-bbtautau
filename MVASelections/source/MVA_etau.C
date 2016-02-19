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


bool ApplyFullEventSelection(const FlatTree* tree, std::vector<Double_t>& vars, FlatTree::EventCategory selectedCategory,
                             MVA_Histograms& h)
{
    if(tree->byCombinedIsolationDeltaBetaCorrRaw3Hits2 >= 1.5
            || tree->pfRelIso1 >= 0.1
            || tree->Q1 * tree->Q2 != -1 || tree->mt1 >= 30)
        return false;

    if(selectedCategory != FlatTree::TwoJets_TwoBtag && selectedCategory != FlatTree::TwoJets_OneBtag &&
            selectedCategory != FlatTree::TwoJets_ZeroBtag)
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

    vars.assign(15, 0);
    vars[0] = tree->pt1;
    vars[1] = tree->pt2;
    vars[2] = b1.Pt();
    vars[3] = b2.Pt();
    vars[4] = b1.DeltaR(b2);
    vars[5] = MET.DeltaPhi(BB);
    vars[6] = tree->drLT;
    vars[7]=  TT.Pt();
    vars[8] = TT.DeltaR(BB);
    vars[9] = BB.Pt();
    vars[10] = MET.DeltaPhi(TT);
    vars[11] = H.Pt();
    vars[12] = tree->mt2;
    vars[13] = tree->mt1;
    vars[14] = TT_MET.Pt();

    h.pt_l1->Fill(vars[0]);
    h.pt_l2->Fill(vars[1]);
    h.pt_b1->Fill(vars[2]);
    h.pt_b2->Fill(vars[3]);
    h.dR_bb->Fill(vars[4] );
    h.dPhi_bb_MET->Fill(vars[5] );
    h.dR_ll->Fill(vars[6] );
    h.pt_Htt->Fill(vars[7] );
    h.dR_Hbb_Htt->Fill(vars[8] );
    h.pt_Hbb->Fill(vars[9] );
    h.dPhi_MET_tt->Fill(vars[10] );
    h.pt_H->Fill(vars[11] );
    h.mT2->Fill(vars[12] );
    h.mT1->Fill(vars[13]);
    h.pt_Htt_MET->Fill(vars[14]);

    return true;
}

void MVA_etau(const std::string& filePath, const std::string& category_name, const std::string& suffix)
{
    std::cout << "==> Start TMVAClassification" << std::endl;

    FlatTree::EventCategory eventCategory = FlatTree::EventCategoryFromString(category_name);
    std::cout << "Training for event category: " << category_name << std::endl;

    const std::string output_name = "TMVA_etau_" + category_name + "_" + suffix;
    TFile* outputFile = TFile::Open( (output_name + ".root").c_str(), "RECREATE" );

    TMVA::Factory *factory = new TMVA::Factory( output_name, outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");


    factory->AddVariable("pt_mu", 'F');
    factory->AddVariable("pt_tau", 'F');
    factory->AddVariable("pt_b1", 'F');
    factory->AddVariable("pt_b2", 'F');
    factory->AddVariable("DR_bb", 'F');
    factory->AddVariable("DPhi_BBMET", 'F');
    factory->AddVariable("DR_ll", 'F');
    factory->AddVariable("Pt_Htt", 'F');
    factory->AddVariable("DR_HBBHTT", 'F');
    factory->AddVariable("Pt_Hbb", 'F');
    factory->AddVariable("DeltaPhi_METTT", 'F');
    factory->AddVariable("PtH", 'F');
    factory->AddVariable("mT2", 'F');
    factory->AddVariable("mT1", 'F');
    factory->AddVariable("Pt_Htt_MET", 'F');

    FlatTree* sigTree = new FlatTree();
    FlatTree* bkgTree = new FlatTree();

    sigTree->Add(filePath+"ggH_hh_bbtautau_*.root");
    bkgTree->Add(filePath+"tt_*.root");

    outputFile->cd();

    ApplySelection(factory, sigTree, eventCategory, true);
    ApplySelection(factory, bkgTree, eventCategory, false);

    TMVAtestAndTraining(factory);

    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;

    delete factory;
    delete sigTree;
    delete bkgTree;
}





