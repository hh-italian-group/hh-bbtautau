/*! Definition of a containers that store analyzers output.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/variadic.hpp>
#include <boost/bimap.hpp>
#include "AnalysisTools/Core/include/SmartTree.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "EventAnalyzerDataId.h"
#include "h-tautau/Core/include/TauIdResults.h"

namespace analysis {

#define CREATE_VAR(r, type, name) VAR(type, name)
#define VAR_LIST(type, ...) BOOST_PP_SEQ_FOR_EACH(CREATE_VAR, type, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

#define CREATE_SVAR(r, type, name) SVAR(type, name)
#define SVAR_LIST(type, ...) BOOST_PP_SEQ_FOR_EACH(CREATE_SVAR, type, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

#define P4_DATA(name) \
    VAR_LIST(float, name##_pt, name##_eta, name##_phi, name##_m) \
    /**/

#define TAU_DATA(name) \
    P4_DATA(name) \
    VAR_LIST(float, name##_iso, name##_DeepTauVSe, name##_DeepTauVSmu, name##_DeepTauVSjet) \
    VAR_LIST(int, name##_q, name##_gen_match) \
    SVAR_LIST(float, name##_q_float, name##_gen_match_float) \
    /**/

#define JET_DATA(name) \
    P4_DATA(name) \
    VAR_LIST(float, name##_CSV, name##_DeepCSV, name##_DeepFlavour, name##_DeepFlavour_CvsL, name##_DeepFlavour_CvsB, \
                    name##_HHbtag) \
    VAR_LIST(int, name##_valid, name##_hadronFlavour) \
    SVAR_LIST(float, name##_hadronFlavour_float) \

#define ANA_EVENT_DATA() \
    VAR(std::vector<size_t>, dataIds) /* EventAnalyzerDataId */ \
    VAR(std::vector<double>, all_weights) /* all weight */ \
    VAR(std::vector<float>, all_mva_scores) /* all mva scores */ \
    VAR_LIST(bool, has_b_pair, has_VBF_pair) /* has 2 jets */ \
    VAR(bool, pass_VBF_trigger) \
    VAR(bool, is_central_es) \
    VAR(ULong64_t, sample_id) \
    VAR(double, weight) /* weight */ \
    VAR(float, mva_score) /* mva score */ \
    VAR(ULong64_t, evt)  /* event */ \
    VAR_LIST(UInt_t, run, lumi) \
    TAU_DATA(tau1) \
    TAU_DATA(tau2) \
    JET_DATA(b1) \
    JET_DATA(b2) \
    JET_DATA(VBF1) \
    JET_DATA(VBF2) \
    VAR_LIST(float, MET_pt, MET_phi) /* MET */ \
    VAR(int, SVfit_valid) \
    SVAR(float, SVfit_valid_float) \
    P4_DATA(SVfit) \
    VAR_LIST(float, SVfit_pt_error, SVfit_eta_error, SVfit_phi_error, SVfit_m_error, SVfit_mt, SVfit_mt_error) \
    VAR(int, kinFit_convergence) \
    SVAR(float, kinFit_convergence_float) \
    VAR_LIST(float, kinFit_m, kinFit_chi2) \
    VAR(float, MT2) \
    VAR_LIST(float, npv, HT_total, HT_otherjets, lhe_HT, n_jets, n_jets_eta24, n_jets_eta24_eta5, \
                    n_selected_gen_jets, n_selected_gen_bjets, genJets_nTotal, \
                    jets_nTotal_hadronFlavour_b, jets_nTotal_hadronFlavour_c) \
    SVAR_LIST(float, m_tt_vis, pt_H_tt, eta_H_tt, phi_H_tt, pt_H_tt_MET, mt_1, mt_2, dR_l1l2, abs_dphi_l1MET, \
                     dphi_htautauMET, dR_l1l2MET, dR_l1l2Pt_htautau, mass_l1l2MET, pt_l1l2MET, MT_htautau, p_zeta, \
                     p_zetavisible, mt_tot) \
    SVAR_LIST(float, m_bb, pt_H_bb, dphi_hbbhtautau, deta_hbbhtautau, costheta_METhbb, dR_b1b2, dR_b1b2_boosted, \
                     dR_lj, mass_top1, mass_top2) \
    SVAR_LIST(float, hh_btag_b1b2, hh_btag_b1_minus_b2, hh_btag_VBF1VBF2) \
    /**/

namespace bbtautau {

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
#define SVAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
struct AnaEvent : public root_ext::detail::BaseDataClass {
    ANA_EVENT_DATA()
    void UpdateSecondaryVariables();
};
#undef VAR
#undef SVAR

#define VAR(type, name) AddBranch(#name, _data->name);
#define SVAR(type, name) AddVarPointer(#name, &_data->name);

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
#undef SVAR
#undef ANA_EVENT_DATA
#undef VAR_LIST
#undef CREATE_VAR
#undef SVAR_LIST
#undef CREATE_SVAR
#undef TAU_DATA
#undef JET_DATA
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

    AnaTupleWriter(const std::string& file_name, Channel channel, bool _runKinFit, bool _runSVfit, bool _allow_calc_svFit);
    ~AnaTupleWriter();
    void AddEvent(EventInfo& event, const DataIdMap& dataIds, const bool& pass_vbf_trigger);

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

    AnaTupleReader(const std::string& file_name, Channel channel, NameSet& active_var_names);
    const DataId& GetDataIdByHash(Hash hash) const;
    const DataId& GetDataIdByIndex(size_t dataId_index) const;
    AnaTuple& GetAnaTuple();
    void UpdateSecondaryBranches(const DataId& dataId, size_t dataId_index);

private:
    void ExtractDataIds(const AnaAux& aux);
    void ExtractMvaRanges(const AnaAux& aux);
    float GetNormalizedMvaScore(const DataId& dataId, float raw_score) const;

private:
    std::shared_ptr<TFile> file;
    std::shared_ptr<AnaTuple> tuple;
    DataIdBiMap known_data_ids;
    NameSet var_branch_names;
    RangeMap mva_ranges;
    Range mva_target_range{0., 0.99999};
};
} // namespace bbtautau
} // namespace analysis
