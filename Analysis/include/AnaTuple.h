/*! Definition of a containers that store analyzers output.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/variadic.hpp>
#include "AnalysisTools/Core/include/SmartTree.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "EventAnalyzerDataId.h"

namespace analysis {

enum class Sample_Index { Signal_NonRes = -125, Signal_Radion_M250 = -250,
                          Signal_Radion_M260 = -260, Signal_Radion_M270 = -270, Signal_Radion_M280 = -280,
                          Signal_Radion_M300 = -300, Signal_Radion_M320 = -320, Signal_Radion_M340 = -340,
                          Signal_Radion_M350 = -350, Signal_Radion_M400 = -400, Signal_Radion_M450 = -450,
                          Signal_Radion_M500 = -500, Signal_Radion_M550 = -550, Signal_Radion_M600 = -600,
                          Signal_Radion_M650 = -650, Signal_Radion_M750 = -750, Signal_Radion_M800 = -800,
                          Signal_Radion_M900 = -900, Data = 0, TT = 1, DY = 2, Wjets= 3, SM_Higgs = 4, other_bkg = 5 };
ENUM_NAMES(Sample_Index) = {
        { Sample_Index::Signal_NonRes, "Signal_NonRes" },
        { Sample_Index::Signal_Radion_M250, "Signal_Radion_M250" },
        { Sample_Index::Signal_Radion_M260, "Signal_Radion_M260" },
        { Sample_Index::Signal_Radion_M270, "Signal_Radion_M270" },
        { Sample_Index::Signal_Radion_M280, "Signal_Radion_M280" },
        { Sample_Index::Signal_Radion_M300, "Signal_Radion_M300" },
        { Sample_Index::Signal_Radion_M320, "Signal_Radion_M320" },
        { Sample_Index::Signal_Radion_M340, "Signal_Radion_M340" },
        { Sample_Index::Signal_Radion_M350, "Signal_Radion_M350" },
        { Sample_Index::Signal_Radion_M400, "Signal_Radion_M400" },
        { Sample_Index::Signal_Radion_M450, "Signal_Radion_M450" },
        { Sample_Index::Signal_Radion_M500, "Signal_Radion_M500" },
        { Sample_Index::Signal_Radion_M550, "Signal_Radion_M550" },
        { Sample_Index::Signal_Radion_M600, "Signal_Radion_M600" },
        { Sample_Index::Signal_Radion_M650, "Signal_Radion_M650" },
        { Sample_Index::Signal_Radion_M750, "Signal_Radion_M750" },
        { Sample_Index::Signal_Radion_M800, "Signal_Radion_M800" },
        { Sample_Index::Signal_Radion_M900, "Signal_Radion_M900" },
        { Sample_Index::Data, "Data" },
        { Sample_Index::TT, "TT" }, { Sample_Index::DY, "DY" }, { Sample_Index::Wjets, "Wjets" },
        { Sample_Index::SM_Higgs, "SM_Higgs" }, { Sample_Index::other_bkg, "other_bkg" }
};

enum class Region_Index { Other = 0, OS_Isolated = 1, OS_AntiIsolated = 2, SS_Isolated = 3, SS_AntiIsolated = 4 };

#define CREATE_VAR(r, type, name) VAR(type, name)
#define VAR_LIST(type, ...) BOOST_PP_SEQ_FOR_EACH(CREATE_VAR, type, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

#define ANA_EVENT_DATA() \
    VAR(std::vector<size_t>, dataIds) /* EventAnalyzerDataId */ \
    VAR(std::vector<double>, all_weights) /* all weight */ \
    VAR(std::vector<float>, all_mva_scores) /* all mva scores */ \
    VAR(bool, has_2jets) /* has 2 jets */ \
    VAR(double, weight) /* weight */ \
    VAR(int, sample_id) /* sample_id */ \
    VAR(float, mva_score) /* mva score */ \
    VAR(uint16_t, iso_wp_1) \
    VAR(uint16_t, iso_wp_2) \
    VAR(int, region_id) \
    VAR_LIST(float, m_ttbb, m_ttbb_kinfit, m_sv, pt_sv, eta_sv, phi_sv, MT2, mt_tot, deta_hbbhtautau, dphi_hbbhtautau, m_tt_vis, pt_H_tt, \
             pt_H_tt_MET, pt_1, eta_1, phi_1, E_1, iso_1, q_1, mt_1, pt_2, eta_2, phi_2, E_2, iso_2, q_2, mt_2, dR_l1l2, abs_dphi_l1MET, \
             dphi_htautauMET, dR_l1l2MET, dR_l1l2Pt_htautau, mass_l1l2MET, pt_l1l2MET, MT_htautau, npv, MET, phiMET, \
             pt_MET, m_bb, pt_H_bb, pt_b1, eta_b1, phi_b1, E_b1, csv_b1, pt_b2, eta_b2, phi_b2, E_b2, csv_b2, costheta_METhbb, dR_b1b2, \
             dR_b1b2_boosted, HT_otherjets, HT_total, n_jets) \
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

    AnaTupleWriter(const std::string& file_name, Channel channel) :
        file(root_ext::CreateRootFile(file_name, ROOT::kLZ4, 5)), tuple(ToString(channel), file.get(), false),
        aux_tuple(file.get(), false)
    {
    }

    ~AnaTupleWriter()
    {
        for(const auto& id : known_data_ids.left) {
            aux_tuple().dataIds.push_back(id.second);
            aux_tuple().dataId_names.push_back(id.first.GetName());
        }
        for(const auto& range : mva_ranges) {
            aux_tuple().mva_selections.push_back(static_cast<unsigned>(range.first));
            aux_tuple().mva_min.push_back(range.second.min());
            aux_tuple().mva_max.push_back(range.second.max());
        }
        aux_tuple.Fill();
        aux_tuple.Write();
        tuple.Write();
    }

    static Sample_Index GetSampleId(const std::string& sample_name)
    {
        const auto vector_name = analysis::SplitValueList(sample_name,true,"_");
        const auto iter_nonres = std::find(vector_name.begin(), vector_name.end(), "NonRes");
        const auto iter_signal = std::find(vector_name.begin(), vector_name.end(), "Radion");
        const auto iter_data = std::find(vector_name.begin(), vector_name.end(), "Data");
        const auto iter_DY = std::find(vector_name.begin(), vector_name.end(), "DY");
        if (iter_nonres != vector_name.end())
            return Sample_Index::Signal_NonRes;
        else if (iter_signal != vector_name.end())
            return analysis::Parse<Sample_Index>("Signal_Radion_"+ vector_name.at(2));
        else if (iter_data != vector_name.end())
            return Sample_Index::Data;
        else if (sample_name == "TT")
            return Sample_Index::TT;
        else if (iter_DY != vector_name.end())
            return Sample_Index::DY;
        else if (sample_name == "Wjets")
            return Sample_Index::Wjets;
        else if (sample_name == "ZH")
            return Sample_Index::SM_Higgs;
        else
            return Sample_Index::other_bkg;
    }

    static Region_Index GetRegionId(const EventRegion& region)
    {
        static const std::map<EventRegion, Region_Index> regions = {
            { EventRegion::OS_Isolated(), Region_Index::OS_Isolated },
            { EventRegion::OS_AntiIsolated(), Region_Index::OS_AntiIsolated },
            { EventRegion::SS_Isolated(), Region_Index::SS_Isolated },
            { EventRegion::SS_AntiIsolated(), Region_Index::SS_AntiIsolated },
        };
        const auto iter = regions.find(region);
        if(iter == regions.end()) return Region_Index::Other;
        return iter->second;
    }

    template<typename EventInfo>
    void AddEvent(EventInfo& event, const DataIdMap& dataIds)
    {
        using namespace ROOT::Math::VectorUtil;
        static constexpr float def_val = std::numeric_limits<float>::lowest();

        if(!dataIds.size()) return;
        Region_Index region_id = Region_Index::Other;
        for(const auto& entry : dataIds) {
            if(!known_data_ids.left.count(entry.first)) {
                const size_t hash = std::hash<std::string>{}(entry.first.GetName());
                if(known_data_ids.right.count(hash))
                    throw exception("Duplicated hash for event id '%1%' and '%2%'.") % entry.first
                        %  known_data_ids.right.at(hash);
                known_data_ids.insert({entry.first, hash});
            }

            SelectionCut mva_cut;
            if(entry.first.Get<EventSubCategory>().TryGetLastMvaCut(mva_cut)) {
                mva_ranges[mva_cut] = mva_ranges[mva_cut].Extend(std::get<1>(entry.second));
            }

            tuple().dataIds.push_back(known_data_ids.left.at(entry.first));
            tuple().all_weights.push_back(std::get<0>(entry.second));
            tuple().all_mva_scores.push_back(static_cast<float>(std::get<1>(entry.second)));

            const Region_Index entry_region_id = GetRegionId(entry.first.Get<EventRegion>());
            if(entry_region_id != Region_Index::Other) {
                if(region_id != Region_Index::Other && region_id != entry_region_id)
                    throw exception("Overlapping regions.");
                region_id = entry_region_id;
            }
        }
        const std::string sample_name = dataIds.begin()->first.Get<std::string>();
        //std::cout << "sample_name in AnaTuple: " << sample_name << std::endl;

        //std::cout << "sample_Id: " << static_cast<int>(GetSampleId(sample_name)) << std::endl;
//        std::cout << "weight, first: " << std::get<0>(dataIds.begin()->second) << ", second: " <<
//                     std::get<1>(dataIds.begin()->second) << std::endl;

        tuple().weight = std::get<0>(dataIds.begin()->second);
        tuple().sample_id = static_cast<int>(GetSampleId(sample_name));
        tuple().region_id = static_cast<int>(region_id);
        tuple().mva_score = def_val;
        tuple().has_2jets = event.HasBjetPair();
        tuple().m_sv = static_cast<float>(event.GetHiggsTTMomentum(true).M());
        tuple().pt_sv = static_cast<float>(event.GetHiggsTTMomentum(true).pt());
        tuple().eta_sv = static_cast<float>(event.GetHiggsTTMomentum(true).eta());
        tuple().phi_sv = static_cast<float>(event.GetHiggsTTMomentum(true).phi());

        if(event.HasBjetPair()) {
            tuple().m_ttbb = static_cast<float>(event.GetResonanceMomentum(true, false).M());
            const auto& kinfit = event.GetKinFitResults();
            tuple().m_ttbb_kinfit = kinfit.HasValidMass() ? static_cast<float>(kinfit.mass) : def_val;
            tuple().MT2 = static_cast<float>(event.GetMT2());
        } else {
            tuple().m_ttbb = def_val;
            tuple().m_ttbb_kinfit = def_val;
            tuple().MT2 = def_val;
        }

        const auto& Htt = event.GetHiggsTTMomentum(false);
        const auto& t1 = event.GetFirstLeg();
        const auto& t2 = event.GetSecondLeg();

        tuple().m_tt_vis = static_cast<float>(Htt.M());
        tuple().pt_H_tt = static_cast<float>(Htt.Pt());
        tuple().pt_H_tt_MET = static_cast<float>((Htt + event.GetMET().GetMomentum()).Pt());
        tuple().pt_1 = static_cast<float>(t1.GetMomentum().pt());
        tuple().eta_1 = static_cast<float>(t1.GetMomentum().eta());
        tuple().phi_1 = static_cast<float>(t1.GetMomentum().phi());
        tuple().E_1 = static_cast<float>(t1.GetMomentum().E());
        tuple().q_1 = event->q_1;
        tuple().iso_1 = static_cast<float>(t1.GetIsolation());
        tuple().iso_wp_1 = t1->iso_wp().GetResultBits();
        tuple().mt_1 = static_cast<float>(Calculate_MT(t1.GetMomentum(), event.GetMET().GetMomentum()));
        tuple().pt_2 = static_cast<float>(t2.GetMomentum().pt());
        tuple().eta_2 = static_cast<float>(t2.GetMomentum().eta());
        tuple().phi_2 = static_cast<float>(t2.GetMomentum().phi());
        tuple().E_2 = static_cast<float>(t2.GetMomentum().E());
        tuple().q_2 = event->q_2;
        tuple().iso_2 = static_cast<float>(t2.GetIsolation());
        tuple().iso_wp_2 = t2->iso_wp().GetResultBits();
        tuple().mt_2 = static_cast<float>(Calculate_MT(t2.GetMomentum(), event.GetMET().GetMomentum()));
        tuple().dR_l1l2 = static_cast<float>(DeltaR(t1.GetMomentum(),t2.GetMomentum()));
        tuple().abs_dphi_l1MET = static_cast<float>(std::abs(DeltaPhi(t1.GetMomentum(), event.GetMET().GetMomentum())));
        tuple().dphi_htautauMET = static_cast<float>(DeltaPhi(event.GetHiggsTTMomentum(true),
                                                              event.GetMET().GetMomentum()));
        tuple().dR_l1l2MET = static_cast<float>(DeltaR(event.GetHiggsTTMomentum(false), event.GetMET().GetMomentum()));
        tuple().dR_l1l2Pt_htautau = static_cast<float>(DeltaR(t1.GetMomentum(), t2.GetMomentum())
                                                       * event.GetHiggsTTMomentum(true).pt());
        tuple().mass_l1l2MET = static_cast<float>((event.GetHiggsTTMomentum(false) + event.GetMET().GetMomentum()).M());
        tuple().pt_l1l2MET = static_cast<float>((event.GetHiggsTTMomentum(false) + event.GetMET().GetMomentum()).pt());
        tuple().MT_htautau = static_cast<float>(Calculate_MT(event.GetHiggsTTMomentum(true),
                                                             event.GetMET().GetMomentum()));
        tuple().npv = event->npv;
        tuple().MET = static_cast<float>(event.GetMET().GetMomentum().Pt());
        tuple().phiMET = static_cast<float>(event.GetMET().GetMomentum().Phi());
        tuple().pt_MET = static_cast<float>(event.GetMET().GetMomentum().pt());
        tuple().mt_tot = static_cast<float>(Calculate_TotalMT(t1.GetMomentum(),t2.GetMomentum(),
                                                              event.GetMET().GetMomentum()));
        tuple().HT_otherjets = event->ht_other_jets;

        tuple().n_jets = event->n_jets;

        if(event.HasBjetPair()) {
            const auto& Hbb = event.GetHiggsBB();
            const auto& b1 = Hbb.GetFirstDaughter();
            const auto& b2 = Hbb.GetSecondDaughter();
            tuple().m_bb = static_cast<float>(Hbb.GetMomentum().M());
            tuple().pt_H_bb = static_cast<float>(Hbb.GetMomentum().Pt());
            tuple().pt_b1 = static_cast<float>(b1.GetMomentum().pt());
            tuple().eta_b1 = static_cast<float>(b1.GetMomentum().Eta());
            tuple().phi_b1 = static_cast<float>(b1.GetMomentum().phi());
            tuple().E_b1 = static_cast<float>(b1.GetMomentum().E());
            tuple().csv_b1 = b1->csv();
            tuple().pt_b2 = static_cast<float>(b2.GetMomentum().Pt());
            tuple().eta_b2 = static_cast<float>(b2.GetMomentum().Eta());
            tuple().phi_b2 = static_cast<float>(b2.GetMomentum().phi());
            tuple().E_b2 = static_cast<float>(b2.GetMomentum().E());
            tuple().csv_b2 = b2->csv();
            tuple().dphi_hbbhtautau = static_cast<float>(DeltaPhi(Hbb.GetMomentum(), event.GetHiggsTTMomentum(true)));
            tuple().deta_hbbhtautau = static_cast<float>((Hbb.GetMomentum()-event.GetHiggsTTMomentum(true)).Eta());
            tuple().costheta_METhbb = static_cast<float>(four_bodies::Calculate_cosTheta_2bodies(
                                                             event.GetMET().GetMomentum(), Hbb.GetMomentum()));
            tuple().dR_b1b2 = static_cast<float>(DeltaR(b1.GetMomentum(), b2.GetMomentum()));
            tuple().dR_b1b2_boosted = static_cast<float>(four_bodies::Calculate_dR_boosted(
                                                             b1.GetMomentum(), b2.GetMomentum(), Hbb.GetMomentum()));
            tuple().HT_total = static_cast<float>(b1.GetMomentum().pt() + b2.GetMomentum().Pt() + event->ht_other_jets);
        } else {
            tuple().m_bb = def_val;
            tuple().pt_H_bb = def_val;
            tuple().pt_b1 = def_val;
            tuple().eta_b1 = def_val;
            tuple().phi_b1 = def_val;
            tuple().E_b1 = def_val;
            tuple().csv_b1 = def_val;
            tuple().pt_b2 = def_val;
            tuple().eta_b2 = def_val;
            tuple().phi_b2 = def_val;
            tuple().E_b2 = def_val;
            tuple().csv_b2 = def_val;
            tuple().dphi_hbbhtautau = def_val;
            tuple().deta_hbbhtautau = def_val;
            tuple().costheta_METhbb = def_val;
            tuple().dR_b1b2 = def_val;
            tuple().dR_b1b2_boosted = def_val;
        }

        tuple.Fill();
    }

private:
    std::shared_ptr<TFile> file;
    AnaTuple tuple;
    AnaAuxTuple aux_tuple;
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

    AnaTupleReader(const std::string& file_name, Channel channel, NameSet& active_var_names) :
        file(root_ext::OpenRootFile(file_name))
    {
        static const NameSet essential_branches = { "dataIds", "all_weights", "has_2jets", "sample_id" , "weight" };
        static const NameSet other_branches = { "all_mva_scores" };
        NameSet enabled_branches;
        if(active_var_names.size()) {
            enabled_branches.insert(essential_branches.begin(), essential_branches.end());
            enabled_branches.insert(active_var_names.begin(), active_var_names.end());
            if(active_var_names.count("mva_score"))
                enabled_branches.insert("all_mva_scores");

        } else {
            enabled_branches = ExtractAllBranchNames(channel);
            for(const auto& name : enabled_branches) {
                if(!essential_branches.count(name) && !other_branches.count(name))
                    active_var_names.insert(name);
            }
        }
        tuple = std::make_shared<AnaTuple>(ToString(channel), file.get(), true, NameSet(), enabled_branches);

        AnaAuxTuple aux_tuple(file.get(), true);
        aux_tuple.GetEntry(0);
        ExtractDataIds(aux_tuple());
        ExtractMvaRanges(aux_tuple());
    }

    const DataId& GetDataIdByHash(Hash hash) const
    {
        const auto iter = known_data_ids.right.find(hash);
        if(iter == known_data_ids.right.end())
            throw exception("EventAnalyzerDataId not found for hash = %1%") % hash;
        return iter->second;
    }

    const DataId& GetDataIdByIndex(size_t dataId_index) const
    {
        if(dataId_index >= (*tuple)().dataIds.size())
            throw exception("DataId index is out of range.");
        return GetDataIdByHash((*tuple)().dataIds.at(dataId_index));
    }

    AnaTuple& GetAnaTuple() { return *tuple; }

    void UpdateSecondaryBranches(const DataId& dataId, size_t dataId_index)
    {
        if(dataId_index >= (*tuple)().dataIds.size())
            throw exception("DataId index is out of range.");
        if(dataId_index >= (*tuple)().all_weights.size())
            throw exception("Inconsistent AnaTuple entry.");
        (*tuple)().weight = (*tuple)().all_weights.at(dataId_index);
        if(dataId_index < (*tuple)().all_mva_scores.size()) {
            const auto raw_score = (*tuple)().all_mva_scores.at(dataId_index);
            (*tuple)().mva_score =  GetNormalizedMvaScore(dataId, raw_score);
        } else {
            (*tuple)().mva_score = 0;
        }
    }

private:
    void ExtractDataIds(const AnaAux& aux)
    {
        const size_t N = aux.dataIds.size();
        if(aux.dataId_names.size() != N)
            throw exception("Inconsistent dataId info in AnaAux tuple.");
        for(size_t n = 0; n < N; ++n) {
            const auto hash = aux.dataIds.at(n);
            const auto& dataId_name = aux.dataId_names.at(n);
            const auto dataId = DataId::Parse(dataId_name);
            if(known_data_ids.right.count(hash))
                throw exception("Duplicated hash = %1% in AnaAux tuple for dataId = %2%.\n"
                                "This hash is already assigned to %3%") % hash % dataId_name
                                % known_data_ids.right.at(hash);
            if(known_data_ids.left.count(dataId))
                throw exception("Duplicated dataId = '%1%' in AnaAux tuple.") % dataId_name;
            known_data_ids.insert({dataId, hash});
        }
    }

    void ExtractMvaRanges(const AnaAux& aux)
    {
        const size_t N = aux.mva_selections.size();
        if(aux.mva_min.size() != N || aux.mva_max.size() != N)
            throw exception("Inconsistent mva range info in AnaAux tuple.");
        for(size_t n = 0; n < N; ++n) {
            const SelectionCut sel = static_cast<SelectionCut>(aux.mva_selections.at(n));
            if(mva_ranges.count(sel))
                throw exception("Duplicated mva selection = %1% in AnaAux tuple.") % sel;
            mva_ranges[sel] = Range(aux.mva_min.at(n), aux.mva_max.at(n));
        }
    }

    NameSet ExtractAllBranchNames(Channel channel) const
    {
        AnaTuple tmpTuple(ToString(channel), file.get(), true);
        return tmpTuple.GetActiveBranches();
    }

    float GetNormalizedMvaScore(const DataId& dataId, float raw_score) const
    {
        SelectionCut sel;
        if(!dataId.Get<EventSubCategory>().TryGetLastMvaCut(sel))
            return raw_score;
        auto iter = mva_ranges.find(sel);
        if(iter == mva_ranges.end())
            throw exception("Mva range not found for %1%.") % sel;
        const Range& range = iter->second;
        const double result = mva_target_range.size() / range.size() * (raw_score - range.min())
                + mva_target_range.min();
        return static_cast<float>(result);
    }

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
