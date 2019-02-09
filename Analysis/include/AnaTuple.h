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

namespace analysis {

#define CREATE_VAR(r, type, name) VAR(type, name)
#define VAR_LIST(type, ...) BOOST_PP_SEQ_FOR_EACH(CREATE_VAR, type, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

#define ANA_EVENT_DATA() \
    VAR(std::vector<size_t>, dataIds) /* EventAnalyzerDataId */ \
    VAR(std::vector<double>, all_weights) /* all weight */ \
    VAR(std::vector<float>, all_mva_scores) /* all mva scores */ \
    VAR(bool, has_2jets) /* has 2 jets */ \
    VAR(double, weight) /* weight */ \
    VAR(float, mva_score) /* mva score */ \
    VAR_LIST(float, m_ttbb, m_ttbb_kinfit, m_sv, MT2, mt_tot, deta_hbbhtautau, dphi_hbbhtautau, m_tt_vis, pt_H_tt, \
             pt_H_tt_MET, pt_1, eta_1, iso_1, mt_1, pt_2, eta_2, iso_2, mt_2, dR_l1l2, abs_dphi_l1MET, \
             dphi_htautauMET, dR_l1l2MET, dR_l1l2Pt_htautau, mass_l1l2MET, pt_l1l2MET, MT_htautau, npv, MET, phiMET, \
             pt_MET, m_bb, pt_H_bb, pt_b1, eta_b1, csv_b1, deepcsv_b1, pt_b2, eta_b2, csv_b2, deepcsv_b2, costheta_METhbb, dR_b1b2, \
             dR_b1b2_boosted, HT_otherjets, mass_top1, mass_top2, p_zeta, p_zetavisible, HT_total, HT_otherjets_gen,\
             HT_total_gen, n_selected_gen_jets, n_selected_gen_bjets, n_selected_gen_notbjets, genJets_nTotal, \
             jets_nTotal_hadronFlavour_b, jets_nTotal_hadronFlavour_c, n_jets, gen_match_1, gen_match_2) \
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

#define ANA_AUX_DATA() \
    VAR(std::vector<size_t>, dataIds) /* EventAnalyzerDataId code */ \
    VAR(std::vector<std::string>, dataId_names) /* EventAnalyzerDataId name */ \
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
    using Range = ::analysis::Range<double>;
    using RangeMap = std::map<SelectionCut, Range>;

    AnaTupleWriter(const std::string& file_name, Channel channel, bool _runKinFit);
    ~AnaTupleWriter();
    void AddEvent(EventInfoBase& event, const DataIdMap& dataIds);

private:
    std::shared_ptr<TFile> file;
    AnaTuple tuple;
    AnaAuxTuple aux_tuple;
    bool runKinFit;
    DataIdBiMap known_data_ids;
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
    NameSet ExtractAllBranchNames(Channel channel) const;
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
