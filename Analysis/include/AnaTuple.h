/*! Definition of a containers that store analyzers output.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "EventAnalyzerDataId.h"

namespace analysis {

#define CREATE_VAR(r, type, name) VAR(type, name)
#define VAR_LIST(type, ...) BOOST_PP_SEQ_FOR_EACH(CREATE_VAR, type, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

#define ANA_EVENT_DATA() \
    VAR(std::vector<uint32_t>, dataIds) /* EventAnalyzerDataId */ \
    VAR(std::vector<double>, all_weights) /* all weight */ \
    VAR(std::vector<float>, all_mva_scores) /* all mva scores */ \
    VAR(bool, has_2jets) /* has 2 jets */ \
    VAR(double, weight) /* weight */ \
    VAR(float, mva_score) /* mva score */ \
    VAR_LIST(float, m_ttbb, m_ttbb_kinfit, m_sv, MT2, mt_tot, deta_hbbhtautau, dphi_hbbhtautau, m_tt_vis, pt_H_tt, \
             pt_H_tt_MET, pt_1, eta_1, iso_1, mt_1, pt_2, eta_2, iso_2, mt_2, dR_l1l2, abs_dphi_l1MET, \
             dphi_htautauMET, dR_l1l2MET, dR_l1l2Pt_htautau, mass_l1l2MET, pt_l1l2MET, MT_htautau, npv, MET, phiMET, \
             pt_MET, m_bb, pt_H_bb, pt_b1, eta_b1, csv_b1, pt_b2, eta_b2, csv_b2, costheta_METhbb, dR_b1b2, \
             dR_b1b2_boosted, HT_otherjets, mass_top1, mass_top2, p_zeta, p_zetavisible) \
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
    VAR(std::vector<uint32_t>, dataIds) /* EventAnalyzerDataId code */ \
    VAR(std::vector<std::string>, dataId_names) /* EventAnalyzerDataId name */ \
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
    using DataIdBiMap = boost::bimap<DataId, uint32_t>;
    using DataIdMap = std::map<DataId, std::tuple<double, double>>;

    AnaTupleWriter(const std::string& file_name, Channel channel) :
        file(root_ext::CreateRootFile(file_name)), tuple(ToString(channel), file.get(), false),
        aux_tuple(file.get(), false)
    {
    }

    ~AnaTupleWriter()
    {
        for(const auto& id : known_data_ids.left) {
            aux_tuple().dataIds.push_back(id.second);
            aux_tuple().dataId_names.push_back(id.first.GetName());
        }
        aux_tuple.Fill();
        aux_tuple.Write();
        tuple.Write();
    }

    template<typename EventInfo>
    void AddEvent(EventInfo& event, const DataIdMap& dataIds)
    {
        using namespace ROOT::Math::VectorUtil;
        static constexpr float def_val = std::numeric_limits<float>::lowest();

        if(!dataIds.size()) return;
        for(const auto& entry : dataIds) {
            if(!known_data_ids.left.count(entry.first)) {
                const uint32_t hash = tools::hash(entry.first.GetName());
                if(known_data_ids.right.count(hash))
                    throw exception("Duplicated hash for event id '%1%' and '%2%'.") % entry.first
                        %  known_data_ids.right.at(hash);
                known_data_ids.insert({entry.first, hash});
            }
            tuple().dataIds.push_back(known_data_ids.left.at(entry.first));
            tuple().all_weights.push_back(std::get<0>(entry.second));
            tuple().all_mva_scores.push_back(static_cast<float>(std::get<1>(entry.second)));
        }

        tuple().weight = def_val;
        tuple().mva_score = def_val;
        tuple().has_2jets = event.HasBjetPair();
        tuple().m_sv = event.GetHiggsTTMomentum(true).M();

        if(event.HasBjetPair()) {
            tuple().m_ttbb = event.GetResonanceMomentum(true, false).M();
            const auto& kinfit = event.GetKinFitResults();
            tuple().m_ttbb_kinfit = kinfit.HasValidMass() ? kinfit.mass : def_val;
            tuple().MT2 = event.GetMT2();
        } else {
            tuple().m_ttbb = def_val;
            tuple().m_ttbb_kinfit = def_val;
            tuple().MT2 = def_val;
        }

        const auto& Htt = event.GetHiggsTTMomentum(false);
        const auto& t1 = event.GetLeg(1);
        const auto& t2 = event.GetLeg(2);

        tuple().m_tt_vis = Htt.M();
        tuple().pt_H_tt = Htt.Pt();
        tuple().pt_H_tt_MET = (Htt + event.GetMET().GetMomentum()).Pt();
        tuple().pt_1 = t1.GetMomentum().pt();
        tuple().eta_1 = t1.GetMomentum().eta();
        tuple().iso_1 = t1.GetIsolation();
        tuple().mt_1 = Calculate_MT(t1.GetMomentum(), event.GetMET().GetMomentum());
        tuple().pt_2 = t2.GetMomentum().pt();
        tuple().eta_2 = t2.GetMomentum().eta();
        tuple().iso_2 = t2.GetIsolation();
        tuple().mt_2 = Calculate_MT(t2.GetMomentum(), event.GetMET().GetMomentum());
        tuple().dR_l1l2 = DeltaR(t1.GetMomentum(),t2.GetMomentum());
        tuple().abs_dphi_l1MET = std::abs(DeltaPhi(t1.GetMomentum(), event.GetMET().GetMomentum()));
        tuple().dphi_htautauMET = DeltaPhi(event.GetHiggsTTMomentum(true), event.GetMET().GetMomentum());
        tuple().dR_l1l2MET = DeltaR(event.GetHiggsTTMomentum(false), event.GetMET().GetMomentum());
        tuple().dR_l1l2Pt_htautau = DeltaR(t1.GetMomentum(), t2.GetMomentum()) * event.GetHiggsTTMomentum(true).pt();
        tuple().mass_l1l2MET = (event.GetHiggsTTMomentum(false) + event.GetMET().GetMomentum()).M();
        tuple().pt_l1l2MET = (event.GetHiggsTTMomentum(false) + event.GetMET().GetMomentum()).pt();
        tuple().MT_htautau = Calculate_MT(event.GetHiggsTTMomentum(true), event.GetMET().GetMomentum());
        tuple().npv = event->npv;
        tuple().MET = event.GetMET().GetMomentum().Pt();
        tuple().phiMET = event.GetMET().GetMomentum().Phi();
        tuple().pt_MET = event.GetMET().GetMomentum().pt();
        tuple().p_zeta = Calculate_Pzeta(t1.GetMomentum(), t2.GetMomentum(),  event.GetMET().GetMomentum());
        tuple().p_zetavisible = Calculate_visiblePzeta(t1.GetMomentum(), t2.GetMomentum());
        tuple().mt_tot = Calculate_TotalMT(t1.GetMomentum(),t2.GetMomentum(),event.GetMET().GetMomentum());
        tuple().HT_otherjets = event->ht_other_jets;

        if(event.HasBjetPair()) {
            const auto& Hbb = event.GetHiggsBB();
            const auto& b1 = Hbb.GetFirstDaughter();
            const auto& b2 = Hbb.GetSecondDaughter();
            tuple().m_bb = Hbb.GetMomentum().M();
            tuple().pt_H_bb = Hbb.GetMomentum().Pt();
            tuple().pt_b1 = b1.GetMomentum().pt();
            tuple().eta_b1 = b1.GetMomentum().Eta();
            tuple().csv_b1 = b1->csv();
            tuple().pt_b2 = b2.GetMomentum().Pt();
            tuple().eta_b2 = b2.GetMomentum().Eta();
            tuple().csv_b2 = b2->csv();
            tuple().dphi_hbbhtautau = DeltaPhi(Hbb.GetMomentum(), event.GetHiggsTTMomentum(true));
            tuple().deta_hbbhtautau = (Hbb.GetMomentum()-event.GetHiggsTTMomentum(true)).Eta();
            tuple().costheta_METhbb = four_bodies::Calculate_cosTheta_2bodies(event.GetMET().GetMomentum(),
                                                                              Hbb.GetMomentum());
            tuple().dR_b1b2 = DeltaR(b1.GetMomentum(), b2.GetMomentum());
            tuple().dR_b1b2_boosted = four_bodies::Calculate_dR_boosted(b1.GetMomentum(), b2.GetMomentum(),
                                                                        Hbb.GetMomentum());

            tuple().mass_top1 = four_bodies::Calculate_topPairMasses(t1.GetMomentum(), t2.GetMomentum(),
                                                                     b1.GetMomentum(), b2.GetMomentum(),
                                                                     event.GetMET().GetMomentum()).first;
            tuple().mass_top2 = four_bodies::Calculate_topPairMasses(t1.GetMomentum(), t2.GetMomentum(),
                                                                     b1.GetMomentum(), b2.GetMomentum(),
                                                                     event.GetMET().GetMomentum()).second;
        } else {
            tuple().m_bb = def_val;
            tuple().pt_H_bb = def_val;
            tuple().pt_b1 = def_val;
            tuple().eta_b1 = def_val;
            tuple().csv_b1 = def_val;
            tuple().pt_b2 = def_val;
            tuple().eta_b2 = def_val;
            tuple().csv_b2 = def_val;
            tuple().dphi_hbbhtautau = def_val;
            tuple().deta_hbbhtautau = def_val;
            tuple().costheta_METhbb = def_val;
            tuple().dR_b1b2 = def_val;
            tuple().dR_b1b2_boosted = def_val;
            tuple().mass_top1 = def_val;
            tuple().mass_top2 = def_val;
        }

        tuple.Fill();
    }

private:
    std::shared_ptr<TFile> file;
    AnaTuple tuple;
    AnaAuxTuple aux_tuple;
    DataIdBiMap known_data_ids;
};

class AnaTupleReader {
public:
    using DataId = EventAnalyzerDataId;
    using Hash = uint32_t;
    using DataIdBiMap = boost::bimap<DataId, Hash>;
    using DataIdMap = std::map<DataId, std::tuple<double, double>>;
    using NameSet = std::set<std::string>;

    AnaTupleReader(const std::string& file_name, Channel channel, NameSet& active_var_names) :
        file(root_ext::OpenRootFile(file_name))
    {
        static const NameSet essential_branches = { "dataIds", "all_weights", "has_2jets" };
        static const NameSet other_branches = { "all_mva_scores", "weight" };
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

    void UpdateSecondaryBranches(size_t dataId_index)
    {
        if(dataId_index >= (*tuple)().dataIds.size())
            throw exception("DataId index is out of range.");
        if(dataId_index >= (*tuple)().all_weights.size())
            throw exception("Inconsistent AnaTuple entry.");
        (*tuple)().weight = (*tuple)().all_weights.at(dataId_index);
        (*tuple)().mva_score = dataId_index < (*tuple)().all_mva_scores.size()
                             ? (*tuple)().all_mva_scores.at(dataId_index) : 0;
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
                throw exception("Duplicated hash = %1% for in AnaAux tuple.") % hash;
            if(known_data_ids.left.count(dataId))
                throw exception("Duplicated dataId = '%1%' for in AnaAux tuple.") % dataId_name;
            known_data_ids.insert({dataId, hash});
        }
    }

    NameSet ExtractAllBranchNames(Channel channel) const
    {
        AnaTuple tmpTuple(ToString(channel), file.get(), true);
        return tmpTuple.GetActiveBranches();
    }

private:
    std::shared_ptr<TFile> file;
    std::shared_ptr<AnaTuple> tuple;
    DataIdBiMap known_data_ids;
    NameSet var_branch_names;
};
} // namespace bbtautau
} // namespace analysis
