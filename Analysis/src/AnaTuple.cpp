/*! Definition of a containers that store analyzers output.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/AnaTuple.h"

namespace analysis {
namespace bbtautau {

AnaTupleWriter::AnaTupleWriter(const std::string& file_name, Channel channel, bool _runKinFit, bool _runSVfit,
                               bool _allow_calc_svFit) :
    file(root_ext::CreateRootFile(file_name, ROOT::kLZ4, 4)), tuple(ToString(channel), file.get(), false),
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

void AnaTupleWriter::AddEvent(EventInfo& event, const DataIdMap& dataIds, const bool pass_VBF_trigger,
                              const CategoriesFlags& categories_flags,
                              const std::map<DiscriminatorWP, std::map<UncertaintyScale, float>>& btag_weights,
                              const std::map<UncertaintySource, std::map<UncertaintyScale, float>>& uncs_weight_map)
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
        // if(*sample_id != data_id.Get<std::string>()) {
        //     throw exception("Single event has two distinct sample ids: %1% and %2%") % (*sample_id)
        //             % data_id.Get<std::string>();
        // }

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

    auto fill_unc_weight_vec = [](const std::map<UncertaintyScale, float>& weights_in, std::vector<float>& weights_out,
                                  bool store_relative = true) {
        static const std::vector<UncertaintyScale> scales = { UncertaintyScale::Up, UncertaintyScale::Down };

        weights_out.clear();
        const float central_w = weights_in.at(UncertaintyScale::Central);
        if(!store_relative)
            weights_out.push_back(central_w);
        for(auto scale : scales) {
            if(!weights_in.count(scale)) break;
            float w = weights_in.at(scale);
            if(store_relative)
                w /= central_w;
            weights_out.push_back(w);
        }
    };

    if(!event->isData){
        fill_unc_weight_vec(btag_weights.at(DiscriminatorWP::Loose), tuple().btag_weight_Loose, false);
        fill_unc_weight_vec(btag_weights.at(DiscriminatorWP::Medium), tuple().btag_weight_Medium, false);
        fill_unc_weight_vec(btag_weights.at(DiscriminatorWP::Tight), tuple().btag_weight_Tight, false);
        if(event.GetEventCandidate().GetUncSource() == UncertaintySource::None){
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::EleTriggerUnc), tuple().unc_EleTriggerUnc);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::MuonTriggerUnc), tuple().unc_MuonTriggerUnc);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauTriggerUnc_DM0), tuple().unc_TauTriggerUnc_DM0);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauTriggerUnc_DM1), tuple().unc_TauTriggerUnc_DM1);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauTriggerUnc_DM10), tuple().unc_TauTriggerUnc_DM10);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauTriggerUnc_DM11), tuple().unc_TauTriggerUnc_DM11);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauVSjetSF_DM0), tuple().unc_TauVSjetSF_DM0);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauVSjetSF_DM1), tuple().unc_TauVSjetSF_DM1);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauVSjetSF_3prong), tuple().unc_TauVSjetSF_3prong);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauVSjetSF_pt20to25), tuple().unc_TauVSjetSF_pt20to25);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauVSjetSF_pt25to30), tuple().unc_TauVSjetSF_pt25to30);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauVSjetSF_pt30to35), tuple().unc_TauVSjetSF_pt30to35);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauVSjetSF_pt35to40), tuple().unc_TauVSjetSF_pt35to40);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauVSjetSF_ptgt40), tuple().unc_TauVSjetSF_ptgt40);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauVSeSF_barrel), tuple().unc_TauVSeSF_barrel);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauVSeSF_endcap), tuple().unc_TauVSeSF_endcap);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauVSmuSF_etaLt0p4), tuple().unc_TauVSmuSF_etaLt0p4);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauVSmuSF_eta0p4to0p8), tuple().unc_TauVSmuSF_eta0p4to0p8);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauVSmuSF_eta1p2to1p7), tuple().unc_TauVSmuSF_eta1p2to1p7);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauVSmuSF_eta1p2to1p7), tuple().unc_TauVSmuSF_eta1p2to1p7);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TauVSmuSF_etaGt1p7), tuple().unc_TauVSmuSF_etaGt1p7);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::EleIdIsoUnc), tuple().unc_EleIdIsoUnc);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::MuonIdIsoUnc), tuple().unc_MuonIdIsoUnc);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::TopPt), tuple().unc_TopPt, true);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::L1_prefiring), tuple().unc_L1_prefiring, true);
            fill_unc_weight_vec(uncs_weight_map.at(UncertaintySource::PileUp), tuple().unc_PileUp, true);
        }
    }
    tuple().has_b_pair = event.HasBjetPair();
    tuple().has_VBF_pair = event.HasVBFjetPair();
    tuple().pass_VBF_trigger = pass_VBF_trigger;
    tuple().run = event->run;
    tuple().lumi = event->lumi;
    tuple().evt = event->evt;
    tuple().channelId = event->channelId;
    tuple().num_central_jets =  static_cast<Int_t>(categories_flags.num_jets);
    tuple().num_btag_loose = static_cast<Int_t>(categories_flags.num_btag_loose);
    tuple().num_btag_medium = static_cast<Int_t>(categories_flags.num_btag_medium);
    tuple().num_btag_tight = static_cast<Int_t>(categories_flags.num_btag_tight);
    tuple().is_vbf = categories_flags.is_vbf;
    tuple().is_boosted = categories_flags.is_boosted;

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
        tuple().name##_decay_mode = obj->decayMode();
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

    auto centralJets = event.GetCentralJets();
    std::sort(centralJets.begin(), centralJets.end(), [](auto jet1, auto jet2) {
        return jet1->GetMomentum().pt() > jet2->GetMomentum().pt();
    });

    auto forwardJets = event.GetForwardJets();
    std::sort(forwardJets.begin(), forwardJets.end(), [](auto jet1, auto jet2) {return jet1->GetMomentum().pt() > jet2->GetMomentum().pt(); });

    //remove forward jets with pt < 30
    auto low_pt_iter = std::find_if(forwardJets.begin(), forwardJets.end(), [](auto jet) {
        return jet->GetMomentum().Pt() <= cuts::hh_bbtautau_Run2::jetID::vbf_pt; });
    forwardJets.erase(low_pt_iter, forwardJets.end());

    //remove from forward and central signal jets
    const std::set<const JetCandidate*> signal_jets = { b1, b2, vbf1, vbf2 };
    auto is_signal = [&signal_jets](const JetCandidate* jet) -> bool { return signal_jets.count(jet); };
    centralJets.erase(std::remove_if(centralJets.begin(), centralJets.end(), is_signal), centralJets.end());
    forwardJets.erase(std::remove_if(forwardJets.begin(), forwardJets.end(), is_signal), forwardJets.end());

    centralJets.resize(5, nullptr);
    forwardJets.resize(5, nullptr);

    JET_DATA(central_jet1, centralJets.at(0))
    JET_DATA(central_jet2, centralJets.at(1))
    JET_DATA(central_jet3, centralJets.at(2))
    JET_DATA(central_jet4, centralJets.at(3))
    JET_DATA(central_jet5, centralJets.at(4))


    JET_DATA(forward_jet1, forwardJets.at(0))
    JET_DATA(forward_jet2, forwardJets.at(1))
    JET_DATA(forward_jet3, forwardJets.at(2))
    JET_DATA(forward_jet4, forwardJets.at(3))
    JET_DATA(forward_jet5, forwardJets.at(4))

    #undef JET_DATA

    #define FAT_JET_DATA(name, obj) \
        tuple().name##_pt = obj ? static_cast<float>(obj->GetMomentum().pt()) : def_val; \
        tuple().name##_eta = obj ? static_cast<float>(obj->GetMomentum().eta()) : def_val; \
        tuple().name##_phi = obj ? static_cast<float>(obj->GetMomentum().phi()) : def_val; \
        tuple().name##_m = obj ? static_cast<float>(obj->GetMomentum().M()) : def_val; \
        tuple().name##_valid = obj != nullptr; \
        tuple().name##_m_softDrop = obj ? static_cast<float>((*obj)->m(ntuple::TupleFatJet::MassType::SoftDrop)) : def_val; \
        /**/

    FAT_JET_DATA(fat_jet, categories_flags.fat_jet_cand)

    #undef FAT_JET_DATA

    tuple().MET_pt = static_cast<float>(event.GetMET().GetMomentum().pt());
    tuple().MET_phi = static_cast<float>(event.GetMET().GetMomentum().phi());


    const sv_fit_ana::FitResults* SVfit = nullptr;
    if(runSVfit && event.GetSVFitResults(allow_calc_svFit).has_valid_momentum)
        SVfit = &event.GetSVFitResults(allow_calc_svFit);
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

const AnaTupleReader::NameSet AnaTupleReader::BoolBranches = {
    "has_b_pair", "has_VBF_pair", "pass_VBF_trigger"
};

const AnaTupleReader::NameSet AnaTupleReader::IntBranches = {
    "tau1_q", "tau1_gen_match", "tau1_decay_mode","tau2_q", "tau2_gen_match", "tau2_decay_mode",
    "b1_valid", "b1_hadronFlavour", "b2_valid", "b2_hadronFlavour",
    "VBF1_valid", "VBF1_hadronFlavour", "VBF2_valid", "VBF2_hadronFlavour",
    "SVfit_valid", "kinFit_convergence", "central_jet1_valid", "central_jet1_hadronFlavour", "central_jet2_valid",
    "central_jet2_hadronFlavour", "central_jet3_valid", "central_jet3_hadronFlavour", "central_jet4_valid",
    "central_jet4_hadronFlavour", "central_jet5_valid", "central_jet5_hadronFlavour"
};

AnaTupleReader::AnaTupleReader(const std::string& file_name, Channel channel, NameSet& active_var_names) :
    file(root_ext::OpenRootFile(file_name)), tree(root_ext::ReadObject<TTree>(*file, ToString(channel))),
    dataFrame(*tree), df(dataFrame)
{
    static const NameSet support_branches = {
        "dataIds", "all_weights", "is_central_es", "sample_id", "all_mva_scores", "weight", "btag_weight",
        "evt", "run", "lumi", "tau1_p4", "tau2_p4", "b1_valid", "b1_p4", "b2_valid", "b2_p4", "MET_p4",
        "Hbb_p4", "Htt_p4", "HttMET_p4", "VBF1_valid", "VBF1_p4", "VBF2_valid", "VBF2_p4", "SVfit_p4", "mass_top_pair",
        "is_boosted", "channelId", "central_jet1_valid", "central_jet1_p4", "central_jet2_valid", "central_jet2_p4",
        "central_jet3_valid", "central_jet3_p4", "central_jet4_valid", "central_jet4_p4", "central_jet5_valid",
         "central_jet5_p4"
    };

    DefineBranches(active_var_names, active_var_names.empty());
    if(active_var_names.empty()) {
        std::vector<std::vector<std::string>> names = {
            df.GetColumnNames(),
        };
        for(auto& other_df : skimmed_df)
            names.push_back(other_df.GetColumnNames());
        for(const auto& name_set : names) {
            for(const auto& name : name_set) {
                if(!support_branches.count(name))
                    active_var_names.insert(name);
            }
        }
    }

    AnaAuxTuple aux_tuple(file.get(), true);
    aux_tuple.GetEntry(0);
    ExtractDataIds(aux_tuple());
    ExtractMvaRanges(aux_tuple());
}

void AnaTupleReader::DefineBranches(const NameSet& active_var_names, bool all)
{
    // const auto Define = [&](RDF& target_df, const std::string& var, const std::string& expr) {
    //     if(all || active_var_names.count(var))
    //         target_df = target_df.Define(var, boost::str(boost::format("static_cast<double>(%1%)") % expr));
    // };

    const auto Define = [&](RDF& target_df, const std::string& var, auto expr,
                              const std::vector<std::string>& columns, bool force = false) {
        if(force || all || active_var_names.count(var))
            target_df = target_df.Define(var, expr, columns);
    };

    const auto Filter = [](RDF& target_df, const std::string& var) -> RDF {
        return target_df.Filter([](bool flag){ return flag; }, {var});
    };

    const auto FilterInt = [](RDF& target_df, const std::string& var) -> RDF {
        return target_df.Filter([](int flag) -> bool { return flag; }, {var});
    };

    const auto Sum = [](float a, float b) -> double { return a + b; };
    const auto Delta = [](float a, float b) -> double { return a - b; };

    const auto ReturnP4 = [](float pt, float eta, float phi, float mass) {
        return LorentzVectorM(pt, eta, phi, mass);
    };

    const auto ReturnMETP4 = [](float pt, float phi) {
        return LorentzVectorM(pt, 0, phi, 0);
    };


    const auto SumP4 = [](const LorentzVectorM& p4_1, const LorentzVectorM& p4_2) {
        return p4_1 + p4_2;
    };

    const auto GetPt = [](const LorentzVectorM& p4) { return p4.pt(); };
    const auto GetEta = [](const LorentzVectorM& p4) { return p4.eta(); };
    const auto GetPhi = [](const LorentzVectorM& p4) { return p4.phi(); };
    const auto GetMass = [](const LorentzVectorM& p4) { return p4.mass(); };

    const auto DeltaPhi = [](const LorentzVectorM& p4_1, const LorentzVectorM& p4_2) {
        return ROOT::Math::VectorUtil::DeltaPhi(p4_1, p4_2);
    };
    const auto AbsDeltaPhi = [](const LorentzVectorM& p4_1, const LorentzVectorM& p4_2) {
        return std::abs(ROOT::Math::VectorUtil::DeltaPhi(p4_1, p4_2));
    };
    const auto DeltaEta = [](const LorentzVectorM& p4_1, const LorentzVectorM& p4_2) {
        return p4_1.eta() - p4_2.eta();
    };
    const auto DeltaR = [](const LorentzVectorM& p4_1, const LorentzVectorM& p4_2) {
        return ROOT::Math::VectorUtil::DeltaR(p4_1, p4_2);
    };

    const auto DefineP4 = [&](RDF& target_df, const std::string& prefix) {
        Define(target_df, prefix + "_p4", ReturnP4,
               { prefix + "_pt", prefix + "_eta", prefix + "_phi", prefix + "_m" }, true);
    };

    const auto _Calculate_MT = [](const LorentzVectorM& p4, const LorentzVectorM& MET_p4) {
        return Calculate_MT(p4, MET_p4);
    };

    DefineP4(df, "tau1");
    DefineP4(df, "tau2");
    Define(df, "Htt_p4", SumP4, { "tau1_p4", "tau2_p4" }, true);
    Define(df, "MET_p4", ReturnMETP4, {"MET_pt", "MET_phi"}, true);
    Define(df, "HttMET_p4", SumP4, { "Htt_p4", "MET_p4" }, true);

    auto df_bb = Filter(df, "has_b_pair");
    DefineP4(df_bb, "b1");
    DefineP4(df_bb, "b2");
    Define(df_bb, "Hbb_p4", SumP4, { "b1_p4", "b2_p4" }, true);

    auto df_vbf = Filter(df_bb, "has_VBF_pair");
    DefineP4(df_vbf, "VBF1");
    DefineP4(df_vbf, "VBF2");

    DefineP4(df, "central_jet1");
    DefineP4(df, "central_jet2");
    DefineP4(df, "central_jet3");
    DefineP4(df, "central_jet4");
    DefineP4(df, "central_jet5");

    auto df_sv = FilterInt(df, "SVfit_valid");
    DefineP4(df_sv, "SVfit");

    auto df_bb_sv = FilterInt(df_bb, "SVfit_valid");
    DefineP4(df_bb_sv, "SVfit");

    Define(df, "m_tt_vis", GetMass, {"Htt_p4"});
    Define(df, "pt_H_tt", GetPt, {"Htt_p4"});
    Define(df, "eta_H_tt", GetEta, {"Htt_p4"});
    Define(df, "phi_H_tt", GetPhi, {"Htt_p4"});
    Define(df, "mt_1", _Calculate_MT, {"tau1_p4", "MET_p4"});
    Define(df, "mt_2", _Calculate_MT, {"tau2_p4", "MET_p4"});
    Define(df_sv, "MT_htautau", _Calculate_MT, {"SVfit_p4", "MET_p4"});

    Define(df, "dR_l1l2", DeltaR, {"tau1_p4", "tau2_p4"});
    Define(df, "abs_dphi_l1MET", AbsDeltaPhi, {"tau1_p4", "MET_p4"});
    Define(df_sv, "dphi_htautauMET", DeltaPhi, {"SVfit_p4", "MET_p4"});
    Define(df, "dR_l1l2MET", DeltaR, {"Htt_p4", "MET_p4"});
    Define(df_sv, "dR_l1l2Pt_htautau", [](const LorentzVectorM& p4_1, const LorentzVectorM& p4_2, float pt)
        { return ROOT::Math::VectorUtil::DeltaR(p4_1, p4_2) * pt; }, {"tau1_p4", "tau2_p4", "SVfit_pt"});
    Define(df, "mass_l1l2MET", GetMass, {"HttMET_p4"});
    Define(df, "pt_l1l2MET", GetPt, {"HttMET_p4"});
    Define(df, "p_zeta", [](const LorentzVectorM& tau1_p4, const LorentzVectorM& tau2_p4, const LorentzVectorM& MET_p4)
        { return Calculate_Pzeta(tau1_p4, tau2_p4, MET_p4); }, {"tau1_p4", "tau2_p4", "MET_p4"});
    Define(df, "p_zetavisible", [](const LorentzVectorM& tau1_p4, const LorentzVectorM& tau2_p4)
        { return Calculate_visiblePzeta(tau1_p4, tau2_p4); }, {"tau1_p4", "tau2_p4"});
    Define(df, "mt_tot", [](const LorentzVectorM& tau1_p4, const LorentzVectorM& tau2_p4, const LorentzVectorM& MET_p4)
        { return Calculate_TotalMT(tau1_p4, tau2_p4, MET_p4); }, {"tau1_p4", "tau2_p4", "MET_p4"});

    Define(df_bb, "m_bb", GetMass, {"Hbb_p4"});
    Define(df_bb, "pt_H_bb", GetPt, {"Hbb_p4"});
    Define(df_bb_sv, "dphi_hbbhtautau", DeltaPhi, {"Hbb_p4", "SVfit_p4"});
    Define(df_bb_sv, "deta_hbbhtautau", DeltaEta, {"Hbb_p4", "SVfit_p4"});
    Define(df_bb, "costheta_METhbb", [](const LorentzVectorM& MET_p4, const LorentzVectorM& Hbb_p4)
        { return four_bodies::Calculate_cosTheta_2bodies(MET_p4, Hbb_p4); }, {"MET_p4", "Hbb_p4"});
    Define(df_bb, "dR_b1b2", DeltaR, {"b1_p4", "b2_p4"});
    Define(df_bb, "p_zeta", [](const LorentzVectorM& b1_p4, const LorentzVectorM& b2_p4, const LorentzVectorM& Hbb_p4)
        { return four_bodies::Calculate_dR_boosted(b1_p4, b2_p4, Hbb_p4); }, {"b1_p4", "b2_p4", "Hbb_p4"});
    Define(df_bb, "dR_lj", [](const LorentzVectorM& tau1_p4, const LorentzVectorM& tau2_p4,
                              const LorentzVectorM& b1_p4, const LorentzVectorM& b2_p4)
        { return four_bodies::Calculate_min_dR_lj(tau1_p4, tau2_p4, b1_p4, b2_p4); },
        { "tau1_p4", "tau2_p4", "b1_p4", "b2_p4"});

    Define(df_bb, "mass_top_pair", [](const LorentzVectorM& tau1_p4, const LorentzVectorM& tau2_p4,
                                      const LorentzVectorM& b1_p4, const LorentzVectorM& b2_p4,
                                      const LorentzVectorM& MET_p4)
        { return four_bodies::Calculate_topPairMasses(tau1_p4, tau2_p4, b1_p4, b2_p4, MET_p4); },
        { "tau1_p4", "tau2_p4", "b1_p4", "b2_p4", "MET_p4"});

    Define(df_bb, "mass_top1", [](const std::pair<double, double>& m_top) { return m_top.first; },
        {"mass_top_pair"});
    Define(df_bb, "mass_top2", [](const std::pair<double, double>& m_top) { return m_top.second; },
        {"mass_top_pair"});
    Define(df_bb, "hh_btag_b1b2", Sum, {"b1_HHbtag", "b2_HHbtag"});
    Define(df_bb, "hh_btag_b1_minus_b2", Delta, {"b1_HHbtag", "b2_HHbtag"});
    Define(df_vbf, "hh_btag_VBF1VBF2", Sum, {"VBF1_HHbtag", "VBF2_HHbtag"});

    struct category_storage{
        size_t num_jets;
        size_t num_btag_loose;
        size_t num_btag_medium;
        size_t num_btag_tight;
        bool is_vbf;
        bool is_boosted;
        //const FatJetCandidate* fat_jet_cand;

        category_storage(size_t _num_jets, size_t _num_btag_loose, size_t _num_btag_medium, size_t _num_btag_tight,
                         bool _is_vbf, bool _is_boosted):
             num_jets(_num_jets), num_btag_loose(_num_btag_loose), num_btag_medium(_num_btag_medium),
             num_btag_tight(_num_btag_tight), is_vbf(_is_vbf), is_boosted(_is_boosted) {}

    };

    const auto return_category_storage = [] (size_t num_jets, size_t num_btag_loose, size_t num_btag_medium,
            size_t num_btag_tight, bool is_vbf, bool is_boosted) {
        return category_storage(num_jets, num_btag_loose,num_btag_medium, num_btag_tight,is_vbf, is_boosted);
    };

    Define(df, "category_storage",return_category_storage,{"num_jets", "num_btag_loose", "num_btag_medium", "num_btag_tight",
           "is_vbf", "is_boosted"});


    skimmed_df.push_back(df_bb);
    skimmed_df.push_back(df_vbf);
    skimmed_df.push_back(df_sv);
    skimmed_df.push_back(df_bb_sv);
}

const AnaTupleReader::DataId& AnaTupleReader::GetDataIdByHash(Hash hash) const
{
    const auto iter = known_data_ids.right.find(hash);
    if(iter == known_data_ids.right.end())
        throw exception("EventAnalyzerDataId not found for hash = %1%") % hash;
    return iter->second;
}

size_t AnaTupleReader::GetNumberOfEntries() const { return static_cast<size_t>(tree->GetEntries()); }
const AnaTupleReader::RDF& AnaTupleReader::GetDataFrame() const { return df; }
const std::list<AnaTupleReader::RDF>& AnaTupleReader::GetSkimmedDataFrames() const { return skimmed_df; }

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

std::string HyperPoint::ToString()
{
    std::vector<std::string> points;
    if(spin)
        points.push_back("s" + ::analysis::ToString(*spin));
    if(mass)
        points.push_back("m" + ::analysis::ToString(*mass));
    if(kl)
        points.push_back("kl" + ::analysis::ToString(*kl));
    std::ostringstream ss;
    if(!points.size())
        throw exception("Empty hyperparameter.");
    ss << points.at(0);
    for(size_t n = 1; n < points.size(); ++n)
        ss << "_" << points.at(n);
    return ss.str();
}

} // namespace bbtautau
} // namespace analysis
