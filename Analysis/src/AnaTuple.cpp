/*! Definition of a containers that store analyzers output.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/AnaTuple.h"

namespace analysis {
namespace bbtautau {

AnaTupleWriter::AnaTupleWriter(const std::string& file_name, Channel channel, bool _runKinFit, bool _runSVfit,
                               bool _allow_calc_svFit) :
    file(root_ext::CreateRootFile(file_name)), tuple(ToString(channel), file.get(), false),
    aux_tuple(file.get(), false), runKinFit(_runKinFit), runSVfit(_runSVfit), allow_calc_svFit(_allow_calc_svFit)
{
}


AnaTupleWriter::~AnaTupleWriter()
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

void AnaTupleWriter::AddEvent(EventInfo& event, const AnaTupleWriter::DataIdMap& dataIds)
{
    using namespace ROOT::Math::VectorUtil;
    static constexpr float def_val = std::numeric_limits<float>::lowest();

    if(!dataIds.size()) return;
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
    }

    tuple().weight = def_val;
    tuple().mva_score = def_val;
    tuple().has_2jets = event.HasBjetPair();
    tuple().run = event->run;
    tuple().lumi = event->lumi;
    tuple().evt = event->evt;

    tuple().m_sv = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
                   static_cast<float>(event.GetHiggsTTMomentum(true,allow_calc_svFit).M()) : def_val;
     tuple().m_sv_error = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
                         static_cast<float>(event.GetSVFitResults(allow_calc_svFit).momentum_error.M()) : def_val;
     tuple().mt_sv = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
                  static_cast<float>(event.GetSVFitResults(allow_calc_svFit).transverseMass) : def_val;
     tuple().mt_sv_error = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
                        static_cast<float>(event.GetSVFitResults(allow_calc_svFit).transverseMass_error) : def_val;
    tuple().pt_sv = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
                    static_cast<float>(event.GetHiggsTTMomentum(true,allow_calc_svFit).Pt()) : def_val;
    tuple().pt_sv_error = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
                          static_cast<float>(event.GetSVFitResults(allow_calc_svFit).momentum_error.Pt()) : def_val;
    tuple().eta_sv = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
                     static_cast<float>(event.GetHiggsTTMomentum(true,allow_calc_svFit).Eta()) : def_val;
    tuple().eta_sv_error = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
                          static_cast<float>(event.GetSVFitResults(allow_calc_svFit).momentum_error.Eta()) : def_val;
    tuple().phi_sv = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
                     static_cast<float>(event.GetHiggsTTMomentum(true,allow_calc_svFit).Phi()) : def_val;
    tuple().phi_sv_error = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
                           static_cast<float>(event.GetSVFitResults(allow_calc_svFit).momentum_error.Phi()) : def_val;
    if(event.HasBjetPair()) {
        tuple().m_ttbb = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
            static_cast<float>(event.GetResonanceMomentum(true, false, allow_calc_svFit).M()) : def_val;
        if(runKinFit){
            const auto& kinfit = event.GetKinFitResults(allow_calc_svFit);
            tuple().m_ttbb_kinfit = kinfit.HasValidMass() ? static_cast<float>(kinfit.mass) : def_val;
            tuple().chi2_kinFit = kinfit.HasValidMass() ? static_cast<float>(kinfit.chi2) : def_val;
        }
        tuple().MT2 = static_cast<float>(event.GetMT2());
    } else {
        tuple().m_ttbb = def_val;
        tuple().m_ttbb_kinfit = def_val;
        tuple().MT2 = def_val;
    }
    const auto& Htt = event.GetHiggsTTMomentum(false);
    const auto& t1 = event.GetLeg(1);
    const auto& t2 = event.GetLeg(2);

    tuple().m_tt_vis = static_cast<float>(Htt.M());
    tuple().pt_H_tt = static_cast<float>(Htt.Pt());
    tuple().eta_H_tt = static_cast<float>(Htt.Eta());
    tuple().phi_H_tt = static_cast<float>(Htt.Phi());
    tuple().pt_H_tt_MET = static_cast<float>((Htt + event.GetMET().GetMomentum()).Pt());
    tuple().pt_1 = static_cast<float>(t1.GetMomentum().pt());
    tuple().eta_1 = static_cast<float>(t1.GetMomentum().eta());
    tuple().phi_1 = static_cast<float>(t1.GetMomentum().phi());
    tuple().m_1 = static_cast<float>(t1.GetMomentum().M());
    tuple().iso_1 = static_cast<float>(t1.GetIsolation());
    tuple().mt_1 = static_cast<float>(Calculate_MT(t1.GetMomentum(), event.GetMET().GetMomentum()));
    tuple().pt_2 = static_cast<float>(t2.GetMomentum().pt());
    tuple().eta_2 = static_cast<float>(t2.GetMomentum().eta());
    tuple().phi_2 = static_cast<float>(t2.GetMomentum().phi());
    tuple().m_2 = static_cast<float>(t2.GetMomentum().M());
    tuple().iso_2 = static_cast<float>(t2.GetIsolation());
    if(t2->leg_type() == LegType::tau){
      tuple().deepTau_vs_e_2 = static_cast<float>(t2->GetRawValue(TauIdDiscriminator::byDeepTau2017v2p1VSe));
      tuple().deepTau_vs_mu_2 = static_cast<float>(t2->GetRawValue(TauIdDiscriminator::byDeepTau2017v2p1VSmu));
      tuple().deepTau_vs_jet_2 = static_cast<float>(t2->GetRawValue(TauIdDiscriminator::byDeepTau2017v2p1VSjet));
      tuple().tauId_default = static_cast<float>(t2->GetRawValue(TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017));
    }
    tuple().mt_2 = static_cast<float>(Calculate_MT(t2.GetMomentum(), event.GetMET().GetMomentum()));
    tuple().dR_l1l2 = static_cast<float>(DeltaR(t1.GetMomentum(),t2.GetMomentum()));
    tuple().abs_dphi_l1MET = static_cast<float>(std::abs(DeltaPhi(t1.GetMomentum(), event.GetMET().GetMomentum())));
    tuple().dphi_htautauMET = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
                              static_cast<float>(DeltaPhi(event.GetHiggsTTMomentum(true,allow_calc_svFit), event.GetMET().GetMomentum())) : def_val;
    tuple().dR_l1l2MET = static_cast<float>(DeltaR(event.GetHiggsTTMomentum(false), event.GetMET().GetMomentum()));
    tuple().dR_l1l2Pt_htautau = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
                                static_cast<float>(DeltaR(t1.GetMomentum(), t2.GetMomentum())
                                * event.GetHiggsTTMomentum(true,allow_calc_svFit).pt()) : def_val;
    tuple().mass_l1l2MET = static_cast<float>((event.GetHiggsTTMomentum(false) + event.GetMET().GetMomentum()).M());
    tuple().pt_l1l2MET = static_cast<float>((event.GetHiggsTTMomentum(false) + event.GetMET().GetMomentum()).pt());
    tuple().MT_htautau = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
                         static_cast<float>(Calculate_MT(event.GetHiggsTTMomentum(true,allow_calc_svFit),event.GetMET().GetMomentum())) : def_val;
    tuple().npv = event->npv;
    tuple().MET = static_cast<float>(event.GetMET().GetMomentum().Pt());
    tuple().phiMET = static_cast<float>(event.GetMET().GetMomentum().Phi());
    tuple().pt_MET = static_cast<float>(event.GetMET().GetMomentum().pt());
    tuple().p_zeta = static_cast<float>(Calculate_Pzeta(t1.GetMomentum(), t2.GetMomentum(),
                                                        event.GetMET().GetMomentum()));
    tuple().p_zetavisible = static_cast<float>(Calculate_visiblePzeta(t1.GetMomentum(), t2.GetMomentum()));
    tuple().mt_tot = static_cast<float>(Calculate_TotalMT(t1.GetMomentum(),t2.GetMomentum(),
                                                          event.GetMET().GetMomentum()));
    tuple().HT_otherjets = event->ht_other_jets;
    tuple().HT_otherjets_gen = static_cast<float>(event.GetHT(false,true));
    tuple().HT_total_gen = static_cast<float>(event.GetHT(true,true));
    tuple().n_jets = event.EventInfo::SelectJets(20,5,false,false,JetOrdering::DeepCSV,
                                                     event.EventInfo::GetSelectedBjetIndicesSet()).size();
    tuple().n_jets_pu = event.EventInfo::SelectJets(20,5,true,false,JetOrdering::DeepCSV,
                                                        event.EventInfo::GetSelectedBjetIndicesSet()).size();
    tuple().n_jets_eta24_eta5 = event.EventInfo::SelectJets(20,5,false,false,JetOrdering::DeepCSV,
                                                                event.EventInfo::GetSelectedBjetIndicesSet(),2.4).size();
    tuple().n_jets_eta24_eta5_pu = event.EventInfo::SelectJets(20,5,true,false,JetOrdering::DeepCSV,
                                                                   event.EventInfo::GetSelectedBjetIndicesSet(),2.4).size();
    tuple().n_jets_eta24 = event.EventInfo::SelectJets(20,2.4,false,false,JetOrdering::DeepCSV,
                                                           event.EventInfo::GetSelectedBjetIndicesSet()).size();
    tuple().n_jets_eta24_pu = event.EventInfo::SelectJets(20,2.4,true,false,JetOrdering::DeepCSV,
                                                              event.EventInfo::GetSelectedBjetIndicesSet()).size();

    tuple().pt_VBF_1 = event.HasVBFjetPair() ? static_cast<float>(event.GetVBFJet(1).GetMomentum().Pt()) : def_val;
    tuple().eta_VBF_1 = event.HasVBFjetPair() ? static_cast<float>(event.GetVBFJet(1).GetMomentum().Eta()) : def_val;
    tuple().phi_VBF_1 = event.HasVBFjetPair() ? static_cast<float>(event.GetVBFJet(1).GetMomentum().Phi()) : def_val;
    tuple().m_VBF_1 = event.HasVBFjetPair() ? static_cast<float>(event.GetVBFJet(1).GetMomentum().M()) : def_val;
    tuple().deep_flavour_VBF_1 = event.HasVBFjetPair() ? static_cast<float>(event.GetVBFJet(1)->deepFlavour())
                                                       : def_val;
    tuple().pt_VBF_2 = event.HasVBFjetPair() ? static_cast<float>(event.GetVBFJet(2).GetMomentum().Pt()) : def_val;
    tuple().eta_VBF_2 = event.HasVBFjetPair() ? static_cast<float>(event.GetVBFJet(2).GetMomentum().Eta()) : def_val;
    tuple().phi_VBF_2 = event.HasVBFjetPair() ? static_cast<float>(event.GetVBFJet(2).GetMomentum().Phi()) : def_val;
    tuple().m_VBF_2 = event.HasVBFjetPair() ? static_cast<float>(event.GetVBFJet(2).GetMomentum().M()) : def_val;
    tuple().deep_flavour_VBF_2 = event.HasVBFjetPair() ? static_cast<float>(event.GetVBFJet(2)->deepFlavour())
                                                       : def_val;

    tuple().n_selected_gen_jets =  event->genJets_p4.size();
    int n_bflavour=0;
    int n_otherflavour=0;
    static constexpr double b_Flavour = 5;
    for(size_t i=0;i<event->genJets_hadronFlavour.size();i++){
        if(event->genJets_hadronFlavour.at(i)==b_Flavour) n_bflavour++;
        else n_otherflavour++;
    }
    tuple().n_selected_gen_bjets = n_bflavour;
    tuple().n_selected_gen_notbjets = n_otherflavour;
    tuple().genJets_nTotal = event->genJets_nTotal;
    tuple().jets_nTotal_hadronFlavour_b = event->jets_nTotal_hadronFlavour_b;
    tuple().jets_nTotal_hadronFlavour_c = event->jets_nTotal_hadronFlavour_c;

    tuple().gen_match_1 = static_cast<float>(event.GetLeg(1)->gen_match());
    tuple().gen_match_2 = static_cast<float>(event.GetLeg(2)->gen_match());

    if(event.HasBjetPair()) {
        const auto& Hbb = event.GetHiggsBB();
        const auto& b1 = Hbb.GetFirstDaughter();
        const auto& b2 = Hbb.GetSecondDaughter();
        tuple().m_bb = static_cast<float>(Hbb.GetMomentum().M());
        tuple().pt_H_bb = static_cast<float>(Hbb.GetMomentum().Pt());
        tuple().pt_b1 = static_cast<float>(b1.GetMomentum().pt());
        tuple().eta_b1 = static_cast<float>(b1.GetMomentum().Eta());
        tuple().phi_b1 = static_cast<float>(b1.GetMomentum().Phi());
        tuple().m_b1 = static_cast<float>(b1.GetMomentum().M());
        tuple().csv_b1 = b1->csv();
        tuple().deepcsv_b1 = b1->deepcsv();
        tuple().deep_flavour_b1 = b1->deepFlavour();
        tuple().pt_b2 = static_cast<float>(b2.GetMomentum().Pt());
        tuple().eta_b2 = static_cast<float>(b2.GetMomentum().Eta());
        tuple().phi_b2 = static_cast<float>(b2.GetMomentum().Phi());
        tuple().m_b2 = static_cast<float>(b2.GetMomentum().M());
        tuple().csv_b2 = b2->csv();
        tuple().deepcsv_b2 = b2->deepcsv();
        tuple().deep_flavour_b2 = b2->deepFlavour();
        tuple().dphi_hbbhtautau = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
                                  static_cast<float>(DeltaPhi(Hbb.GetMomentum(),
                                                    event.GetHiggsTTMomentum(true,allow_calc_svFit))) : def_val;
        tuple().deta_hbbhtautau = runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum ?
                                  static_cast<float>((Hbb.GetMomentum() -
                                  event.GetHiggsTTMomentum(true,allow_calc_svFit)).Eta()) : def_val;
        tuple().costheta_METhbb = static_cast<float>(four_bodies::Calculate_cosTheta_2bodies(
                                                         event.GetMET().GetMomentum(), Hbb.GetMomentum()));
        tuple().dR_b1b2 = static_cast<float>(DeltaR(b1.GetMomentum(), b2.GetMomentum()));
        tuple().dR_b1b2_boosted = static_cast<float>(four_bodies::Calculate_dR_boosted(
                                                         b1.GetMomentum(), b2.GetMomentum(), Hbb.GetMomentum()));
        tuple().mass_top1 = static_cast<float>(four_bodies::Calculate_topPairMasses(
                                                   t1.GetMomentum(), t2.GetMomentum(), b1.GetMomentum(),
                                                   b2.GetMomentum(), event.GetMET().GetMomentum()).first);
        tuple().mass_top2 = static_cast<float>(four_bodies::Calculate_topPairMasses(
                                                   t1.GetMomentum(), t2.GetMomentum(), b1.GetMomentum(),
                                                   b2.GetMomentum(), event.GetMET().GetMomentum()).second);
        tuple().HT_total = static_cast<float>(b1.GetMomentum().pt() + b2.GetMomentum().Pt() + event->ht_other_jets);
        tuple().dR_lj = static_cast<float>(four_bodies::Calculate_min_dR_lj(t1.GetMomentum(), t2.GetMomentum(),
                                           b1.GetMomentum(), b2.GetMomentum()));

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
        tuple().dR_lj = def_val;
    }
    tuple.Fill();
}

AnaTupleReader::AnaTupleReader(const std::string& file_name, Channel channel, NameSet& active_var_names) :
    file(root_ext::OpenRootFile(file_name))
{
    static const NameSet essential_branches = { "dataIds", "all_weights", "has_2jets" };
    static const NameSet other_branches = { "all_mva_scores", "weight", "evt", "run", "lumi" };
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

const AnaTupleReader::DataId& AnaTupleReader::GetDataIdByHash(Hash hash) const
{
    const auto iter = known_data_ids.right.find(hash);
    if(iter == known_data_ids.right.end())
        throw exception("EventAnalyzerDataId not found for hash = %1%") % hash;
    return iter->second;
}

const AnaTupleReader::DataId& AnaTupleReader::GetDataIdByIndex(size_t dataId_index) const
{
    if(dataId_index >= (*tuple)().dataIds.size())
        throw exception("DataId index is out of range.");
    return GetDataIdByHash((*tuple)().dataIds.at(dataId_index));
}

AnaTuple& AnaTupleReader::GetAnaTuple() { return *tuple; }

void AnaTupleReader::UpdateSecondaryBranches(const DataId& dataId, size_t dataId_index)
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

void AnaTupleReader::ExtractDataIds(const AnaAux& aux)
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

void AnaTupleReader::ExtractMvaRanges(const AnaAux& aux)
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

AnaTupleReader::NameSet AnaTupleReader::ExtractAllBranchNames(Channel channel) const
{
    AnaTuple tmpTuple(ToString(channel), file.get(), true);
    return tmpTuple.GetActiveBranches();
}

float AnaTupleReader::GetNormalizedMvaScore(const DataId& dataId, float raw_score) const
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

} // namespace bbtautau
} // namespace analysis
