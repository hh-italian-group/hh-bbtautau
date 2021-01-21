/*! Definition of a containers that store analyzers output.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/AnaTuple.h"

namespace analysis {
namespace bbtautau {

AnaTupleWriter::AnaTupleWriter(const std::string& file_name, Channel channel, bool _runKinFit, bool _runSVfit,
                               bool _allow_calc_svFit) :
    file(root_ext::CreateRootFile(file_name, ROOT::kLZMA, 9)), tuple(ToString(channel), file.get(), false),
    aux_tuple(file.get(), false), runKinFit(_runKinFit), runSVfit(_runSVfit), allow_calc_svFit(_allow_calc_svFit)
{
}

AnaTupleWriter::~AnaTupleWriter()
{
    for(const auto& [name, hash] : known_datasets) {
        aux_tuple().dataset_hashes.push_back(hash);
        aux_tuple().dataset_names.push_back(name);
    }

    for(const auto& [name, hash] : known_regions) {
        aux_tuple().region_hashes.push_back(hash);
        aux_tuple().region_names.push_back(name);
    }

    aux_tuple.Fill();
    aux_tuple.Write();
    tuple.Write();
}

void AnaTupleWriter::AddEvent(EventInfo& event, const DataIdMap& dataIds, const CategoriesFlags& categories_flags,
                              const BTagWeights& btag_weights,
                              const std::map<UncertaintySource, std::map<UncertaintyScale, float>>& uncs_weight_map)
{
    static constexpr float def_val = std::numeric_limits<float>::lowest();
    static constexpr int def_val_int = std::numeric_limits<int>::lowest();

    if(!dataIds.size()) return;

    tuple().run = event->run;
    tuple().lumi = event->lumi;
    tuple().evt = event->evt;
    tuple().channelId = event->channelId;
    tuple().is_data = event->isData;

    tuple().has_b_pair = event.HasBjetPair();
    tuple().has_VBF_pair = event.HasVBFjetPair();
    tuple().pass_VBF_trigger = categories_flags.pass_vbf_trigger;
    tuple().pass_only_VBF_trigger = categories_flags.pass_only_vbf_trigger;
    tuple().num_central_jets =  static_cast<Int_t>(categories_flags.num_jets);

    #define SET_NBTAG(col, wp) \
        tuple().col##_##wp = static_cast<Int_t>(categories_flags.col.at(DiscriminatorWP::wp))

    SET_NBTAG(num_btag, Loose);
    SET_NBTAG(num_btag, Medium);
    SET_NBTAG(num_btag, Tight);
    SET_NBTAG(nbtag, Loose);
    SET_NBTAG(nbtag, Medium);
    SET_NBTAG(nbtag, Tight);

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
        tuple().name##_resolution = obj ? (*obj)->resolution() : def_val; \
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

    static const auto order_jets = [] (const auto& jet1, const auto& jet2) {
        return jet1->GetMomentum().pt() > jet2->GetMomentum().pt();
    };

    auto centralJets = event.GetCentralJets();
    std::sort(centralJets.begin(), centralJets.end(), order_jets);

    auto forwardJets = event.GetForwardJets();
    std::sort(forwardJets.begin(), forwardJets.end(), order_jets);

    //remove forward jets with pt < 30
    const auto low_pt_iter = std::find_if(forwardJets.begin(), forwardJets.end(), [](const auto& jet) {
        return jet->GetMomentum().Pt() <= cuts::hh_bbtautau_Run2::jetID::vbf_pt;
    });
    forwardJets.erase(low_pt_iter, forwardJets.end());

    //remove from forward and central signal jets
    const std::set<const JetCandidate*> signal_jets = { b1, b2, vbf1, vbf2 };
    const auto is_signal = [&signal_jets](const JetCandidate* jet) -> bool { return signal_jets.count(jet); };
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
    tuple().MET_cov_00 = static_cast<float>(event.GetMET().GetCovMatrix()(0, 0));
    tuple().MET_cov_01 = static_cast<float>(event.GetMET().GetCovMatrix()(0, 1));
    tuple().MET_cov_11 = static_cast<float>(event.GetMET().GetCovMatrix()(1, 1));


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

    if(!event->isData) {
        #define FILL_BTAG(wp, suffix) \
            FillUncWeightVec(btag_weights.weights.at(DiscriminatorWP::wp), &tuple().weight_btag_##suffix, \
                             &tuple().unc_btag_##suffix##_Up, &tuple().unc_btag_##suffix##_Down)
        FILL_BTAG(Loose, Loose);
        FILL_BTAG(Medium, Medium);
        FILL_BTAG(Tight, Tight);

        #undef FILL_BTAG

        tuple().weight_btag_IterativeFit = btag_weights.iter_weight;
        tuple().weight_btag_IterativeFit_withJES = btag_weights.iter_weight_with_jes ?
                                                   *btag_weights.iter_weight_with_jes : def_val;
        tuple().is_TuneCP5 = btag_weights.is_TuneCP5;
        #define FILL_UNC(r, _, name) \
           FillUncWeightVec(uncs_weight_map.at(UncertaintySource::name), nullptr, \
               &tuple().BOOST_PP_CAT(BOOST_PP_CAT(unc_, name), _Up), \
               &tuple().BOOST_PP_CAT(BOOST_PP_CAT(unc_, name), _Down) );
       #define FILL_UNC_LIST(_, ...) BOOST_PP_SEQ_FOR_EACH(FILL_UNC, _, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

        if(event.GetEventCandidate().GetUncSource() == UncertaintySource::None){
            FILL_UNC_LIST(float,
                EleTriggerUnc, MuonTriggerUnc,
                TauTriggerUnc_DM0, TauTriggerUnc_DM1, TauTriggerUnc_DM10, TauTriggerUnc_DM11,
                TauVSjetSF_DM0, TauVSjetSF_DM1, TauVSjetSF_3prong,
                TauVSjetSF_pt20to25, TauVSjetSF_pt25to30, TauVSjetSF_pt30to35, TauVSjetSF_pt35to40, TauVSjetSF_ptgt40,
                TauVSeSF_barrel, TauVSeSF_endcap,
                TauVSmuSF_etaLt0p4, TauVSmuSF_eta0p4to0p8, TauVSmuSF_eta0p8to1p2, TauVSmuSF_eta1p2to1p7,
                TauVSmuSF_etaGt1p7,
                EleIdIsoUnc, MuonIdIsoUnc,
                TopPt, L1_prefiring, PileUp, PileUpJetId_eff, PileUpJetId_mistag,
                TauCustomSF_DM0, TauCustomSF_DM1, TauCustomSF_DM10, TauCustomSF_DM11, VBFTriggerUnc_jets,
                VBFTauTriggerUnc_DM0, VBFTauTriggerUnc_DM1, VBFTauTriggerUnc_3prong, btag_lf, btag_hf, btag_hfstats1,
                btag_hfstats2, btag_lfstats1, btag_lfstats2, btag_cferr1, btag_cferr2
            )
        }
        #undef FILL_UNC
        #undef FILL_UNC_LIST
    }

    for(const auto& [data_id, weight] : dataIds) {
        const std::string dataset_name = data_id.Get<std::string>();
        const std::string event_region = data_id.Get<EventRegion>().ToString();

        if(!known_datasets.count(dataset_name))
            known_datasets[dataset_name] = static_cast<unsigned>(known_datasets.size()); 
        tuple().dataset = known_datasets.at(dataset_name);

        if(!known_regions.count(event_region))
            known_regions[event_region] = static_cast<unsigned>(known_regions.size());
        tuple().event_region = known_regions.at(event_region);

        tuple().unc_source = static_cast<int>(data_id.Get<UncertaintySource>());
        tuple().unc_scale = static_cast<int>(data_id.Get<UncertaintyScale>());

        tuple().weight = static_cast<float>(weight);

        tuple.Fill();
    }
}

void AnaTupleWriter::FillUncWeightVec(const std::map<UncertaintyScale, float>& weights_in,
                                      float* weight, float* unc_up, float* unc_down)
{
    static constexpr float def_val = std::numeric_limits<float>::lowest();

    const float central_w = weights_in.at(UncertaintyScale::Central);
    if(weight)
        *weight = central_w;

    const auto set_unc = [&](UncertaintyScale scale, float* output) {
        if(output) {
            float w = def_val;
            auto iter = weights_in.find(scale);
            if(iter != weights_in.end())
                w = iter->second / central_w;
            *output = w;
        }
    };

    set_unc(UncertaintyScale::Up, unc_up);
    set_unc(UncertaintyScale::Down, unc_down);
}

} // namespace bbtautau
} // namespace analysis
