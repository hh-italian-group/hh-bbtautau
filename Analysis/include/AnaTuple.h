/*! Definition of a containers that store analyzers output.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/variadic.hpp>
#include <boost/bimap.hpp>
#include <ROOT/RDataFrame.hxx>
#include "AnalysisTools/Core/include/SmartTree.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "EventAnalyzerDataId.h"
#include "h-tautau/Core/include/TauIdResults.h"
#include "h-tautau/Cuts/include/hh_bbtautau_Run2.h"

namespace analysis {

#define CREATE_VAR(r, type, name) VAR(type, name)
#define VAR_LIST(type, ...) BOOST_PP_SEQ_FOR_EACH(CREATE_VAR, type, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

#define P4_DATA(name) \
    VAR_LIST(float, name##_pt, name##_eta, name##_phi, name##_m) \
    /**/

#define TAU_DATA(name) \
    P4_DATA(name) \
    VAR_LIST(float, name##_iso, name##_DeepTauVSe, name##_DeepTauVSmu, name##_DeepTauVSjet) \
    VAR_LIST(int, name##_q, name##_gen_match, name##_decay_mode) \
    /**/

#define JET_DATA(name) \
    P4_DATA(name) \
    VAR_LIST(float, name##_DeepFlavour, name##_DeepFlavour_CvsL, name##_DeepFlavour_CvsB, \
                    name##_HHbtag) \
    VAR_LIST(int, name##_valid, name##_hadronFlavour) \

#define FAT_JET_DATA(name) \
    P4_DATA(name) \
    VAR(float, name##_m_softDrop) \
    VAR(int, name##_valid) \

#define ANA_EVENT_DATA() \
    VAR(std::vector<size_t>, dataIds) /* EventAnalyzerDataId */ \
    VAR(std::vector<double>, all_weights) /* all weight */ \
    VAR(std::vector<float>, all_mva_scores) /* all mva scores */ \
    VAR(std::vector<float>, btag_weight_Loose) /* btag weight Loose wp: central,up,down */ \
    VAR(std::vector<float>, btag_weight_Medium) /* btag weight Medium wp: central,up,down */ \
    VAR(std::vector<float>, btag_weight_Tight) /* btag weight Tight wp: central,up,down */ \
    VAR_LIST(bool, has_b_pair, has_VBF_pair) /* has 2 jets */ \
    VAR(bool, pass_VBF_trigger) \
    VAR(bool, is_central_es) \
    VAR(ULong64_t, sample_id) \
    VAR(double, weight) /* weight */ \
    VAR(float, mva_score) /* mva score */ \
    VAR(ULong64_t, evt)  /* event */ \
    VAR(Int_t, channelId) /* Channel: eTau, muTau or tauTau */ \
    VAR(Int_t, num_central_jets) /* num jets */ \
    VAR(Int_t, num_btag_loose)  \
    VAR(Int_t, num_btag_medium) \
    VAR(Int_t, num_btag_tight)   \
    VAR(bool, is_vbf) \
    VAR(bool, is_boosted) \
    VAR_LIST(UInt_t, run, lumi) \
    TAU_DATA(tau1) \
    TAU_DATA(tau2) \
    JET_DATA(b1) \
    JET_DATA(b2) \
    JET_DATA(VBF1) \
    JET_DATA(VBF2) \
    JET_DATA(central_jet1) \
    JET_DATA(central_jet2) \
    JET_DATA(central_jet3) \
    JET_DATA(central_jet4) \
    JET_DATA(central_jet5) \
    JET_DATA(forward_jet1) \
    JET_DATA(forward_jet2) \
    JET_DATA(forward_jet3) \
    JET_DATA(forward_jet4) \
    JET_DATA(forward_jet5) \
    FAT_JET_DATA(fat_jet) \
    VAR_LIST(float, MET_pt, MET_phi) /* MET */ \
    VAR(int, SVfit_valid) \
    P4_DATA(SVfit) \
    VAR_LIST(float, SVfit_pt_error, SVfit_eta_error, SVfit_phi_error, SVfit_m_error, SVfit_mt, SVfit_mt_error) \
    VAR(int, kinFit_convergence) \
    VAR_LIST(float, kinFit_m, kinFit_chi2) \
    VAR(float, MT2) \
    VAR_LIST(float, npv, HT_total, HT_otherjets, lhe_HT, n_jets, n_jets_eta24, n_jets_eta24_eta5, \
                    n_selected_gen_jets, n_selected_gen_bjets, genJets_nTotal, \
                    jets_nTotal_hadronFlavour_b, jets_nTotal_hadronFlavour_c) \
    VAR_LIST(std::vector<float>, unc_EleTriggerUnc, unc_MuonTriggerUnc, unc_TauTriggerUnc_DM0, unc_TauTriggerUnc_DM1, \
                                 unc_TauTriggerUnc_DM10, unc_TauTriggerUnc_DM11, unc_TauVSjetSF_DM0, unc_TauVSjetSF_DM1, \
                                 unc_TauVSjetSF_3prong, unc_TauVSjetSF_pt20to25, unc_TauVSjetSF_pt25to30, \
                                 unc_TauVSjetSF_pt30to35, unc_TauVSjetSF_pt35to40, unc_TauVSjetSF_ptgt40, \
                                 unc_TauVSeSF_barrel, unc_TauVSeSF_endcap, unc_TauVSmuSF_etaLt0p4, \
                                 unc_TauVSmuSF_eta0p4to0p8, unc_TauVSmuSF_eta0p8to1p2, unc_TauVSmuSF_eta1p2to1p7, \
                                 unc_TauVSmuSF_etaGt1p7, unc_EleIdIsoUnc, unc_MuonIdIsoUnc, unc_TopPt, unc_L1_prefiring, \
                                 unc_PileUp, unc_PileUpJetId_eff, unc_PileUpJetId_mistag, unc_TauCustomSF_DM0, \
                                 unc_TauCustomSF_DM1, unc_TauCustomSF_DM10, unc_TauCustomSF_DM11) 
    /**/

namespace bbtautau {

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
struct AnaEvent : public root_ext::detail::BaseDataClass {
    ANA_EVENT_DATA()
    // void UpdateSecondaryVariables();
};
#undef VAR

#define VAR(type, name) AddBranch(#name, _data->name);

class AnaTuple : public root_ext::detail::BaseSmartTree<AnaEvent> {
public:
    static const std::string& Name() { static const std::string name = "events"; return name; }
    AnaTuple(TDirectory* directory, bool readMode, const std::set<std::string>& disabled_branches = {},
             const std::set<std::string>& enabled_branches = {})
        : BaseSmartTree(Name(), directory, readMode, disabled_branches, enabled_branches) { Initialize(); }
    AnaTuple(const std::string& name, TDirectory* directory, bool readMode,
             const std::set<std::string>& disabled_branches = {}, const std::set<std::string>& enabled_branches = {})
        : BaseSmartTree(name, directory, readMode, disabled_branches,enabled_branches) { Initialize(); }

    const float* GetVarAddress(const std::string& var_name) const
    {
        auto iter = var_pointers.find(var_name);
        if(iter != var_pointers.end())
            return iter->second;
        return &get<float>(var_name);
    }

    const std::set<std::string>& GetSecondaryVariables() const { return var_names; }

private:
    void Initialize()
    {
        ANA_EVENT_DATA()
        if (GetEntries() > 0)
            GetEntry(0);
    }

    void AddVarPointer(const std::string& name, const float* ptr)
    {
        if(var_pointers.count(name))
            throw exception("AnaTuple: duplicated variable '%1%'.") % name;
        var_pointers[name] = ptr;
        var_names.insert(name);
    }

    std::map<std::string, const float* > var_pointers;
    std::set<std::string> var_names;
};

#undef VAR
#undef ANA_EVENT_DATA
#undef VAR_LIST
#undef CREATE_VAR
#undef TAU_DATA
#undef JET_DATA
#undef FAT_JET_DATA
#undef P4_DATA
}

#define ANA_AUX_DATA() \
    VAR(std::vector<size_t>, dataIds) /* EventAnalyzerDataId code */ \
    VAR(std::vector<std::string>, dataId_names) /* EventAnalyzerDataId name */ \
    VAR(std::vector<size_t>, sampleIds) /* sample code */ \
    VAR(std::vector<std::string>, sampleId_names) /* sample name */ \
    VAR(std::vector<unsigned>, mva_selections) /* SelectionCut for mva selection */ \
    VAR(std::vector<double>, mva_min) /* Minimal value for the given mva SelectionCut */ \
    VAR(std::vector<double>, mva_max) /* Maximal value for the given mva SelectionCut */ \
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
    using DataIdBiMap = boost::bimap<DataId, size_t>;
    using DataIdMap = std::map<DataId, std::tuple<double, double>>;
    using SampleIdBiMap = boost::bimap<std::string, size_t>;
    using Range = ::analysis::Range<double>;
    using RangeMap = std::map<SelectionCut, Range>;

    struct CategoriesFlags {
        size_t num_jets;
        size_t num_btag_loose;
        size_t num_btag_medium;
        size_t num_btag_tight;
        bool is_vbf;
        bool is_boosted;
        const FatJetCandidate* fat_jet_cand;
    };

    AnaTupleWriter(const std::string& file_name, Channel channel, bool _runKinFit, bool _runSVfit, bool _allow_calc_svFit);
    ~AnaTupleWriter();
    void AddEvent(EventInfo& event, const DataIdMap& dataIds, const bool pass_VBF_trigger,
                  const CategoriesFlags& categories_flags,
                  const std::map<DiscriminatorWP, std::map<UncertaintyScale, float>>& btag_weights,
                  const std::map<UncertaintySource, std::map<UncertaintyScale, float>>& uncs_weight_map);

private:
    std::shared_ptr<TFile> file;
    AnaTuple tuple;
    AnaAuxTuple aux_tuple;
    bool runKinFit;
    bool runSVfit;
    bool allow_calc_svFit;
    DataIdBiMap known_data_ids;
    SampleIdBiMap known_sample_ids;
    RangeMap mva_ranges;
};

class AnaTupleReader {
public:
    using DataId = EventAnalyzerDataId;
    using Hash = size_t;
    using DataIdBiMap = boost::bimap<DataId, Hash>;
    using DataIdMap = std::map<DataId, std::tuple<double, double>>;
    using NameSet = std::set<std::string>;
    using Range = ::analysis::Range<double>;
    using RangeMap = std::map<SelectionCut, Range>;
    using RDF = ROOT::RDF::RNode;

    static const NameSet BoolBranches, IntBranches;

    struct category_storage{
        float num_jets;
        int num_btag_loose;
        int num_btag_medium;
        int num_btag_tight;
        bool is_vbf;
        bool is_boosted;
        //std::vector<size_t> dataId;
        //std::vector<double> weight;
        //const FatJetCandidate* fat_jet_cand;

        category_storage(){
            num_jets = 999;
            num_btag_loose = -999;
            num_btag_medium = -999;
            num_btag_tight = -999;
            is_vbf = false;
            is_boosted = false;
            //dataId = ;
            //weight = 999;
        }

        category_storage(float _num_jets, int _num_btag_loose, int _num_btag_medium, int _num_btag_tight,
                         bool _is_vbf, bool _is_boosted):
             num_jets(_num_jets), num_btag_loose(_num_btag_loose), num_btag_medium(_num_btag_medium),
             num_btag_tight(_num_btag_tight), is_vbf(_is_vbf), is_boosted(_is_boosted) {}

    };

    AnaTupleReader(const std::string& file_name, Channel channel, NameSet& active_var_names,
                   const std::vector<std::string>& input_friends);
    size_t GetNumberOfEntries() const;
    const DataId& GetDataIdByHash(Hash hash) const;
    const RDF& GetDataFrame() const;
    const std::list<RDF>& GetSkimmedDataFrames() const;
    float GetNormalizedMvaScore(const DataId& dataId, float raw_score) const;
    const std::map<std::string, std::set<std::string>>& GetParametricVariables() const;

private:
    void DefineBranches(const NameSet& active_var_names, bool all);
    void ExtractDataIds(const AnaAux& aux);
    void ExtractMvaRanges(const AnaAux& aux);

    static std::vector<std::shared_ptr<TFile>> OpenFiles(const std::string& file_name,
                                                         const std::vector<std::string>& input_friends);
    static std::vector<std::shared_ptr<TTree>> ReadTrees(Channel channel,
                                                         const std::vector<std::shared_ptr<TFile>>& files);

private:
    std::vector<std::shared_ptr<TFile>> files;
    std::vector<std::shared_ptr<TTree>> trees;
    ROOT::RDataFrame dataFrame;
    RDF df;
    std::list<RDF> skimmed_df;
    DataIdBiMap known_data_ids;
    NameSet var_branch_names;
    RangeMap mva_ranges;
    Range mva_target_range{0., 0.99999};
    std::map<std::string, std::set<std::string>> parametric_vars;
};

struct HyperPoint {
    boost::optional<int> spin;
    boost::optional<double> mass;
    boost::optional<double> kl;

    std::string ToString();
};

} // namespace bbtautau
} // namespace analysis
