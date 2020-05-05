/*! Definition of a containers that store analyzers output.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/AnaTuple.h"

namespace analysis {
namespace bbtautau {

void AnaEvent::UpdateSecondaryVariables()
{
    using namespace ROOT::Math::VectorUtil;
    static constexpr float def_val = std::numeric_limits<float>::lowest();

    #define FVAR(name) name##_float = name;
    #define CREATE_FVAR(r, x, name) FVAR(name)
    #define FVAR_LIST(x, ...) BOOST_PP_SEQ_FOR_EACH(CREATE_FVAR, x, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

    FVAR_LIST(x, tau1_q, tau1_gen_match, tau2_q, tau2_gen_match, b1_hadronFlavour, b2_hadronFlavour, VBF1_hadronFlavour,
              VBF2_hadronFlavour, SVfit_valid, kinFit_convergence)

    #undef FVAR
    #undef CREATE_FVAR
    #undef FVAR_LIST

    const LorentzVectorM t1(tau1_pt, tau1_eta, tau1_phi, tau1_m), t2(tau2_pt, tau2_eta, tau2_phi, tau2_m);
    const auto Htt = t1 + t2;
    const LorentzVectorM MET(MET_pt, 0, MET_phi, 0);

    boost::optional<LorentzVectorM> b1, b2, vbf1, vbf2, Hbb, SVfit;
    if(has_b_pair) {
        b1 = LorentzVectorM(b1_pt, b1_eta, b1_phi, b1_m);
        b2 = LorentzVectorM(b2_pt, b2_eta, b2_phi, b2_m);
        Hbb = *b1 + *b2;
    }
    if(has_VBF_pair) {
        vbf1 = LorentzVectorM(VBF1_pt, VBF1_eta, VBF1_phi, VBF1_m);
        vbf2 = LorentzVectorM(VBF2_pt, VBF2_eta, VBF2_phi, VBF2_m);
    }
    if(SVfit_valid)
        SVfit = LorentzVectorM(SVfit_pt, SVfit_eta, SVfit_phi, SVfit_m);

    m_tt_vis = static_cast<float>(Htt.M());
    pt_H_tt = static_cast<float>(Htt.Pt());
    eta_H_tt = static_cast<float>(Htt.Eta());
    phi_H_tt = static_cast<float>(Htt.Phi());
    pt_H_tt_MET = static_cast<float>((Htt + MET).Pt());
    mt_1 = static_cast<float>(Calculate_MT(t1, MET));
    mt_2 = static_cast<float>(Calculate_MT(t2, MET));
    dR_l1l2 = static_cast<float>(DeltaR(t1, t2));
    abs_dphi_l1MET = static_cast<float>(std::abs(DeltaPhi(t1, MET)));
    dphi_htautauMET = SVfit ? static_cast<float>(DeltaPhi(*SVfit, MET)) : def_val;
    dR_l1l2MET = static_cast<float>(DeltaR(Htt, MET));
    dR_l1l2Pt_htautau = SVfit ? static_cast<float>(DeltaR(t1, t2) * SVfit->pt()) : def_val;
    mass_l1l2MET = static_cast<float>((Htt + MET).M());
    pt_l1l2MET = static_cast<float>((Htt + MET).pt());
    MT_htautau = SVfit ? static_cast<float>(Calculate_MT(*SVfit, MET)) : def_val;
    p_zeta = static_cast<float>(Calculate_Pzeta(t1, t2, MET));
    p_zetavisible = static_cast<float>(Calculate_visiblePzeta(t1, t2));
    mt_tot = static_cast<float>(Calculate_TotalMT(t1, t2, MET));

    m_bb = Hbb ? static_cast<float>(Hbb->mass()) : def_val;
    pt_H_bb = Hbb ? static_cast<float>(Hbb->pt()) : def_val;

    dphi_hbbhtautau = Hbb && SVfit ? static_cast<float>(DeltaPhi(*Hbb, *SVfit)) : def_val;
    deta_hbbhtautau = Hbb && SVfit ? static_cast<float>(((*Hbb) - (*SVfit)).eta()) : def_val;
    costheta_METhbb = Hbb ? static_cast<float>(four_bodies::Calculate_cosTheta_2bodies(MET, *Hbb)) : def_val;
    dR_b1b2 = Hbb ? static_cast<float>(DeltaR(*b1, *b2)) : def_val;
    dR_b1b2_boosted = Hbb ? static_cast<float>(four_bodies::Calculate_dR_boosted(*b1, *b2, *Hbb)) : def_val;
    dR_lj = Hbb ? static_cast<float>(four_bodies::Calculate_min_dR_lj(t1, t2, *b1, *b2)) : def_val;

    boost::optional<std::pair<double, double>> topMasses;
    if(Hbb)
        topMasses = four_bodies::Calculate_topPairMasses(t1, t2, *b1, *b2, MET);
    mass_top1 = topMasses ? static_cast<float>(topMasses->first) : def_val;
    mass_top2 = topMasses ? static_cast<float>(topMasses->second) : def_val;

    hh_btag_b1b2 = Hbb ? b1_HHbtag + b2_HHbtag : def_val;
    hh_btag_b1_minus_b2 = Hbb ? b1_HHbtag - b2_HHbtag : def_val;
    hh_btag_VBF1VBF2 = has_VBF_pair ? VBF1_HHbtag + VBF2_HHbtag : def_val;
}

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
    for(const auto& id : known_sample_ids.left) {
        aux_tuple().sampleIds.push_back(id.second);
        aux_tuple().sampleId_names.push_back(id.first);
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
    static constexpr float def_val = std::numeric_limits<float>::lowest();
    static constexpr int def_val_int = std::numeric_limits<int>::lowest();

    if(!dataIds.size()) return;
    tuple().is_central_es = false;
    tuple().weight = def_val;
    tuple().mva_score = def_val;

    boost::optional<std::string> sample_id;
    for(const auto& entry : dataIds) {
        const DataId& data_id = entry.first;
        const double weight = std::get<0>(entry.second);
        const double mva_score = std::get<1>(entry.second);
        const UncertaintySource unc_source = data_id.Get<UncertaintySource>();

        if(!sample_id) {
            sample_id = data_id.Get<std::string>();
            if(!known_sample_ids.left.count(*sample_id)) {
                const size_t hash = std::hash<std::string>{}(*sample_id);
                if(known_sample_ids.right.count(hash))
                    throw exception("Duplicated hash for sample id '%1%' and '%2%'.") % (*sample_id)
                        %  known_sample_ids.right.at(hash);
                known_sample_ids.insert({*sample_id, hash});
            }
            tuple().sample_id = known_sample_ids.left.at(*sample_id);
        }
        if(*sample_id != data_id.Get<std::string>()) {
            throw exception("Single event has two distinct sample ids: %1% and %2%") % (*sample_id)
                    % data_id.Get<std::string>();
        }

        if(!known_data_ids.left.count(data_id)) {
            const size_t hash = std::hash<std::string>{}(data_id.GetName());
            if(known_data_ids.right.count(hash))
                throw exception("Duplicated hash for event id '%1%' and '%2%'.") % data_id
                    %  known_data_ids.right.at(hash);
            known_data_ids.insert({data_id, hash});
        }

        if(unc_source == UncertaintySource::None) {
            tuple().is_central_es = true;
            tuple().weight = weight;
            tuple().mva_score = static_cast<float>(mva_score);
        }

        SelectionCut mva_cut;
        if(data_id.Get<EventSubCategory>().TryGetLastMvaCut(mva_cut)) {
            mva_ranges[mva_cut] = mva_ranges[mva_cut].Extend(mva_score);
        }

        tuple().dataIds.push_back(known_data_ids.left.at(data_id));
        tuple().all_weights.push_back(weight);
        tuple().all_mva_scores.push_back(static_cast<float>(mva_score));
    }

    tuple().has_b_pair = event.HasBjetPair();
    tuple().has_VBF_pair = event.HasVBFjetPair();
    tuple().run = event->run;
    tuple().lumi = event->lumi;
    tuple().evt = event->evt;

    #define TAU_DATA(name, obj) \
        tuple().name##_pt = static_cast<float>(obj.GetMomentum().pt()); \
        tuple().name##_eta = static_cast<float>(obj.GetMomentum().eta()); \
        tuple().name##_phi = static_cast<float>(obj.GetMomentum().phi()); \
        tuple().name##_m = static_cast<float>(obj.GetMomentum().M()); \
        tuple().name##_iso = obj->leg_type() != LegType::tau ? static_cast<float>(obj.GetIsolation()) : def_val; \
        tuple().name##_DeepTauVSe = obj->leg_type() == LegType::tau \
                                  ? obj->GetRawValue(TauIdDiscriminator::byDeepTau2017v2p1VSe) : def_val; \
        tuple().name##_DeepTauVSmu = obj->leg_type() == LegType::tau \
                                   ? obj->GetRawValue(TauIdDiscriminator::byDeepTau2017v2p1VSmu) : def_val; \
        tuple().name##_DeepTauVSjet = obj->leg_type() == LegType::tau \
                                    ? obj->GetRawValue(TauIdDiscriminator::byDeepTau2017v2p1VSjet) : def_val; \
        tuple().name##_q = obj->charge(); \
        tuple().name##_gen_match = static_cast<int>(obj->gen_match()); \
        /**/

    const auto& t1 = event.GetLeg(1);
    const auto& t2 = event.GetLeg(2);
    TAU_DATA(tau1, t1)
    TAU_DATA(tau2, t2)
    #undef TAU_DATA

    #define JET_DATA(name, obj) \
        tuple().name##_pt = obj ? static_cast<float>(obj->GetMomentum().pt()) : def_val; \
        tuple().name##_eta = obj ? static_cast<float>(obj->GetMomentum().eta()) : def_val; \
        tuple().name##_phi = obj ? static_cast<float>(obj->GetMomentum().phi()) : def_val; \
        tuple().name##_m = obj ? static_cast<float>(obj->GetMomentum().M()) : def_val; \
        tuple().name##_CSV = obj ? (*obj)->csv() : def_val; \
        tuple().name##_DeepCSV = obj ? (*obj)->deepcsv() : def_val; \
        tuple().name##_DeepFlavour = obj ? (*obj)->deepFlavour() : def_val; \
        tuple().name##_DeepFlavour_CvsL = obj ? (*obj)->deepFlavour_CvsL() : def_val; \
        tuple().name##_DeepFlavour_CvsB = obj ? (*obj)->deepFlavour_CvsB() : def_val; \
        tuple().name##_HHbtag = obj ? (*obj)->hh_btag() : def_val; \
        tuple().name##_valid = obj != nullptr; \
        tuple().name##_hadronFlavour = obj ? (*obj)->hadronFlavour() : def_val_int; \
        /**/

    const JetCandidate *b1 = nullptr, *b2 = nullptr, *vbf1 = nullptr, *vbf2 = nullptr;
    if(event.HasBjetPair()) {
        b1 = &event.GetBJet(1);
        b2 = &event.GetBJet(2);
    }
    if(event.HasVBFjetPair()) {
        vbf1 = &event.GetVBFJet(1);
        vbf2 = &event.GetVBFJet(2);
    }
    JET_DATA(b1, b1)
    JET_DATA(b2, b2)
    JET_DATA(VBF1, vbf1)
    JET_DATA(VBF2, vbf2)

    #undef JET_DATA

    tuple().MET_pt = static_cast<float>(event.GetMET().GetMomentum().pt());
    tuple().MET_phi = static_cast<float>(event.GetMET().GetMomentum().phi());


    const sv_fit_ana::FitResults* SVfit = nullptr;
    if(runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum){
        SVfit = &event.GetSVFitResults(allow_calc_svFit);}
    tuple().SVfit_valid = SVfit != nullptr;
    tuple().SVfit_pt = SVfit ? static_cast<float>(SVfit->momentum.pt()) : def_val;
    tuple().SVfit_eta = SVfit ? static_cast<float>(SVfit->momentum.eta()) : def_val;
    tuple().SVfit_phi = SVfit ? static_cast<float>(SVfit->momentum.phi()) : def_val;
    tuple().SVfit_m = SVfit ? static_cast<float>(SVfit->momentum.mass()) : def_val;
    tuple().SVfit_mt = SVfit ? static_cast<float>(SVfit->transverseMass) : def_val;
    tuple().SVfit_pt_error = SVfit ? static_cast<float>(SVfit->momentum_error.pt()) : def_val;
    tuple().SVfit_eta_error = SVfit ? static_cast<float>(SVfit->momentum_error.eta()) : def_val;
    tuple().SVfit_phi_error = SVfit ? static_cast<float>(SVfit->momentum_error.phi()) : def_val;
    tuple().SVfit_m_error = SVfit ? static_cast<float>(SVfit->momentum_error.mass()) : def_val;
    tuple().SVfit_mt_error = SVfit ? static_cast<float>(SVfit->transverseMass_error) : def_val;

    const kin_fit::FitResults* kinFit = nullptr;
    if(runKinFit && event.HasBjetPair())
        kinFit = &event.GetKinFitResults(allow_calc_svFit);
    tuple().kinFit_convergence = kinFit ? kinFit->convergence : def_val_int;
    tuple().kinFit_m = kinFit && kinFit->HasValidMass() ? static_cast<float>(kinFit->mass) : def_val;
    tuple().kinFit_chi2 = kinFit && kinFit->HasValidMass() ? static_cast<float>(kinFit->chi2) : def_val;

    tuple().MT2 = event.HasBjetPair() ? static_cast<float>(event.GetMT2()) : def_val;

    tuple().npv = event->npv;
    tuple().HT_total = static_cast<float>(event.GetHT(true, true));
    tuple().HT_otherjets = static_cast<float>(event.GetHT(false, true));
    tuple().lhe_HT = event->lhe_HT;
    tuple().n_jets = event.GetAllJets().size();
    tuple().n_jets_eta24 = event.GetCentralJets().size();
    tuple().n_jets_eta24_eta5 = event.GetForwardJets().size();

    tuple().n_selected_gen_jets =  event->genJets_p4.size();
    int n_bflavour=0;
    static constexpr double b_Flavour = 5;
    for(size_t i = 0; i < event->genJets_hadronFlavour.size(); ++i) {
        if(event->genJets_hadronFlavour.at(i) == b_Flavour) ++n_bflavour;
    }
    tuple().n_selected_gen_bjets = n_bflavour;
    tuple().genJets_nTotal = event->genJets_nTotal;
    tuple().jets_nTotal_hadronFlavour_b = event->jets_nTotal_hadronFlavour_b;
    tuple().jets_nTotal_hadronFlavour_c = event->jets_nTotal_hadronFlavour_c;

    tuple.Fill();
}

AnaTupleReader::AnaTupleReader(const std::string& file_name, Channel channel, NameSet& active_var_names) :
    file(root_ext::OpenRootFile(file_name))
{
    static const NameSet essential_branches = { "dataIds", "all_weights", "has_b_pair", "has_VBF_pair" };
    static const NameSet other_branches = {
        "is_central_es", "sample_id", "all_mva_scores", "weight", "evt", "run", "lumi", "tau1_q", "tau1_gen_match",
        "tau2_q", "tau2_gen_match", "b1_valid", "b1_hadronFlavour", "b2_valid", "b2_hadronFlavour", "VBF1_valid",
        "VBF1_hadronFlavour", "VBF2_valid", "VBF2_hadronFlavour", "SVfit_valid", "kinFit_convergence"
    };
    tuple = std::make_shared<AnaTuple>(ToString(channel), file.get(), true);
    if(active_var_names.empty()) {
        std::vector<const std::set<std::string>*> names = {
            &tuple->GetActiveBranches(), &tuple->GetSecondaryVariables()
        };
        for(auto name_set : names) {
            for(const auto& name : *name_set) {
                if(!essential_branches.count(name) && !other_branches.count(name))
                    active_var_names.insert(name);
            }
        }
    }

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
