/*! Definition of a containers that store analyzers output.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/variadic.hpp>
#include <boost/bimap.hpp>
#include "AnalysisTools/Core/include/SmartTree.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "EventAnalyzerDataId.h"
#include "h-tautau/Core/include/TauIdResults.h"

namespace analysis {

#define CREATE_VAR(r, type, name) VAR(type, name)
#define VAR_LIST(type, ...) BOOST_PP_SEQ_FOR_EACH(CREATE_VAR, type, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

#define P4_DATA(name) \
    VAR_LIST(float, name##_pt, name##_eta, name##_phi, name##_m) \
    /**/

#define P4_DATA_EX(name) \
    VAR_LIST(float, name##_pt, name##_eta, name##_phi, name##_m) \
    VAR_LIST(float, name##_pt_error, name##_eta_error, name##_phi_error, name##_m_error) \
    /**/

#define TAU_DATA(name) \
    P4_DATA(name) \
    VAR_LIST(float, name##_iso, name##_DeepTauVSe, name##_DeepTauVSmu, name##_DeepTauVSjet) \
    VAR_LIST(int, name##_q, name##_gen_match, name##_decay_mode) \
    /**/

#define JET_DATA(name) \
    P4_DATA(name) \
    VAR_LIST(float, name##_DeepFlavour, name##_DeepFlavour_CvsL, name##_DeepFlavour_CvsB, \
                    name##_HHbtag, name##_resolution) \
    VAR_LIST(int, name##_valid, name##_hadronFlavour) \
    /**/

#define FAT_JET_DATA(name) \
    P4_DATA(name) \
    VAR(float, name##_m_softDrop) \
    VAR(int, name##_valid) \
    /**/

#define UNC_VAR(name) \
    VAR(float, unc_##name##_Up) \
    VAR(float, unc_##name##_Down) \
    /**/

#define WVAR(name) \
    VAR(float, weight_##name) \
    UNC_VAR(name) \
    /**/

#define CREATE_UNC_VAR(r, _, name) UNC_VAR(name)
#define UNC_VAR_LIST(_, ...) BOOST_PP_SEQ_FOR_EACH(CREATE_UNC_VAR, _, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

#define CREATE_WVAR(r, _, name) WVAR(name)
#define WVAR_LIST(_, ...) BOOST_PP_SEQ_FOR_EACH(CREATE_WVAR, _, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))


#define ANA_EVENT_DATA() \
    VAR(UInt_t, dataset) /* hashed dataset name */ \
    VAR(UInt_t, event_region) /* hashed event region name */ \
    VAR(Int_t, unc_source) /* uncertainty source */ \
    VAR(Int_t, unc_scale) /* uncertainty scale */ \
    VAR(float, weight) /* full event weight of the event (excluding category-dependent weights) */ \
    VAR_LIST(UInt_t, run, lumi) /* run ID and lumi ID */ \
    VAR(ULong64_t, evt)  /* event ID */ \
    VAR(int, channelId) /* Channel: 0 - eTau, 1 - muTau, 2 - tauTau, 3 - muMu */ \
    VAR(bool, is_data) /* whatever event is a real data or MC */ \
    VAR_LIST(bool, has_b_pair, has_VBF_pair) /* event has two selected b (VBF) jets */ \
    VAR_LIST(bool, pass_VBF_trigger, pass_only_VBF_trigger) /* even passes (only) VBF trigger */ \
    VAR(int, num_central_jets) /* number of selected central jets in the event */ \
    VAR_LIST(int, num_btag_Loose, num_btag_Medium, num_btag_Tight) /* number of b tagged jets among the two \
                                                                        selected b jets */ \
    VAR_LIST(int, nbtag_Loose, nbtag_Medium, nbtag_Tight) /* number of b tagged jets in the event */ \
    VAR(bool, is_vbf) /* event passes VBF selection */ \
    VAR(bool, is_boosted) /* event passes boosted selection */ \
    TAU_DATA(tau1) /* first tau candidate */ \
    TAU_DATA(tau2) /* second tau candidate */ \
    JET_DATA(b1) /* first b candidate */ \
    JET_DATA(b2) /* second b candidate */ \
    JET_DATA(VBF1) /* first VBF jet candidate */ \
    JET_DATA(VBF2) /* second VBF candidate */ \
    JET_DATA(central_jet1) /* 1st additional central jet in the event */ \
    JET_DATA(central_jet2) /* 2nd additional central jet in the event */ \
    JET_DATA(central_jet3) /* 3rd additional central jet in the event */ \
    JET_DATA(central_jet4) /* 4th additional central jet in the event */ \
    JET_DATA(central_jet5) /* 5th additional central jet in the event */ \
    JET_DATA(forward_jet1) /* 1st additional forward jet in the event */ \
    JET_DATA(forward_jet2) /* 2nd additional forward jet in the event */ \
    JET_DATA(forward_jet3) /* 3rd additional forward jet in the event */ \
    JET_DATA(forward_jet4) /* 4th additional forward jet in the event */ \
    JET_DATA(forward_jet5) /* 5th additional forward jet in the event */ \
    FAT_JET_DATA(fat_jet) /* AK8 jet that matches with the two selected b candidates */ \
    VAR_LIST(float, MET_pt, MET_phi, MET_cov_00, MET_cov_01, MET_cov_11) /* MET momentum and covariance */ \
    VAR(int, SVfit_valid) /* SVfit has valid result */ \
    P4_DATA_EX(SVfit) /* SVfit momentum and its error */ \
    VAR_LIST(float, SVfit_mt, SVfit_mt_error) /* Transverse mass estimated by SVfit and its error */ \
    VAR(int, kinFit_convergence) /* HHKinFit convergence flag: -1 - not converged, 1 - converged, \
                                    2 - converged with b momentum at the limit, \
                                    3 - converged with tau momentum at the limit, \
                                    4 - converged with both b and tau momenta at the limits */ \
    VAR_LIST(float, kinFit_m, kinFit_chi2) /* fourbody mass computed by HHKinfit and chi2 of the fit */ \
    VAR(float, MT2) /* s-transverse mass */ \
    VAR_LIST(float, npv, HT_total, HT_otherjets, lhe_HT, n_jets, n_jets_eta24, n_jets_eta24_eta5, \
                    n_selected_gen_jets, n_selected_gen_bjets, genJets_nTotal, \
                    jets_nTotal_hadronFlavour_b, jets_nTotal_hadronFlavour_c) /* other global event quantities */ \
    WVAR_LIST(_, btag_Loose, btag_Medium, btag_Tight) \
    VAR_LIST(float, weight_btag_IterativeFit, weight_btag_IterativeFit_withJES) \
             /* btag weights for various working points estimated by different algorithms and their unc variations */ \
    UNC_VAR_LIST(_, \
             EleTriggerUnc, MuonTriggerUnc, \
             TauTriggerUnc_DM0, TauTriggerUnc_DM1, TauTriggerUnc_DM10, TauTriggerUnc_DM11, \
             TauVSjetSF_DM0, TauVSjetSF_DM1, TauVSjetSF_3prong, \
             TauVSjetSF_pt20to25, TauVSjetSF_pt25to30, TauVSjetSF_pt30to35, TauVSjetSF_pt35to40, TauVSjetSF_ptgt40, \
             TauVSeSF_barrel, TauVSeSF_endcap, \
             TauVSmuSF_etaLt0p4, TauVSmuSF_eta0p4to0p8, TauVSmuSF_eta0p8to1p2, TauVSmuSF_eta1p2to1p7, \
             TauVSmuSF_etaGt1p7, \
             EleIdIsoUnc, MuonIdIsoUnc, \
             TopPt, L1_prefiring, PileUp, PileUpJetId_eff, PileUpJetId_mistag, \
             TauCustomSF_DM0, TauCustomSF_DM1, TauCustomSF_DM10, TauCustomSF_DM11, VBFTriggerUnc_jets, \
             VBFTauTriggerUnc_DM0, VBFTauTriggerUnc_DM1, VBFTauTriggerUnc_3prong, \
             btag_lf, btag_hf, btag_hfstats1, btag_hfstats2, btag_lfstats1, btag_lfstats2, btag_cferr1, btag_cferr2, \
             btag_JES) \
             /* effects of various uncertainty sources on the event weight */ \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(bbtautau, AnaEvent, AnaTuple, ANA_EVENT_DATA, "events")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(bbtautau, AnaTuple, ANA_EVENT_DATA)
#undef VAR
#undef ANA_EVENT_DATA

#undef VAR_LIST
#undef CREATE_VAR
#undef TAU_DATA
#undef JET_DATA
#undef FAT_JET_DATA
#undef P4_DATA
#undef P4_DATA_EX
#undef UNC_VAR
#undef CREATE_UNC_VAR
#undef UNC_VAR_LIST
#undef WVAR
#undef CREATE_WVAR
#undef WVAR_LIST

#define ANA_AUX_DATA() \
    VAR(std::vector<unsigned>, dataset_hashes) /* hash code for each dataset stored in AnaTuple */ \
    VAR(std::vector<std::string>, dataset_names) /* name of each datasetd stored in AnaTuple */ \
    VAR(std::vector<unsigned>, region_hashes) /* hash code for each event region stored in AnaTuple */ \
    VAR(std::vector<std::string>, region_names) /* name of each datasetd stored in AnaTuple */ \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(bbtautau, AnaAux, AnaAuxTuple, ANA_AUX_DATA, "aux")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(bbtautau, AnaAuxTuple, ANA_AUX_DATA)
#undef VAR
#undef ANA_AUX_DATA

namespace bbtautau {
class AnaTupleWriter {
public:
    using DataId = EventAnalyzerDataId;
    using NameMap = std::map<std::string, unsigned>;
    using DataIdMap = std::map<DataId, double>;

    struct CategoriesFlags {
        size_t num_jets;
        std::map<DiscriminatorWP, size_t> num_btag;
        std::map<DiscriminatorWP, size_t> nbtag;
        bool is_vbf, is_boosted;
        bool pass_vbf_trigger, pass_only_vbf_trigger;
        const FatJetCandidate* fat_jet_cand;
    };

    struct BTagWeights {
        std::map<DiscriminatorWP, std::map<UncertaintyScale, float>> weights;
        std::map<std::pair<UncertaintySource, UncertaintyScale>, float> weights_iter;
        std::map<UncertaintyScale, float> weights_iter_jes;
    };


    AnaTupleWriter(const std::string& file_name, Channel channel, bool _runKinFit, bool _runSVfit,
                   bool _allow_calc_svFit);
    ~AnaTupleWriter();
    void AddEvent(EventInfo& event, const DataIdMap& dataIds, const CategoriesFlags& categories_flags,
                  const BTagWeights& btag_weights, const std::map<UncertaintySource,
                  std::map<UncertaintyScale, float>>& uncs_weight_map);

private:
    static void FillUncWeightVec(const std::map<UncertaintyScale, float>& weights_in,
                                 float* weight, float* unc_up, float* unc_down);
    static void FillUncIter(const BTagWeights& btag_weights, UncertaintySource source,
                            float* unc_up, float* unc_down);

private:
    std::shared_ptr<TFile> file;
    AnaTuple tuple;
    AnaAuxTuple aux_tuple;
    bool runKinFit, runSVfit, allow_calc_svFit;
    NameMap known_datasets, known_regions;
};

} // namespace bbtautau
} // namespace analysis
