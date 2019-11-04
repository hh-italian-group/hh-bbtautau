/*! Definition of functions for calulating most used quantities
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#pragma once

using LorentzVectorM = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >;
using LorentzVectorE = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> >;
using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >;

LorentzVectorM getHTTp4 (const ROOT::VecOps::RVec<LorentzVectorM>& lep_p4, const ROOT::VecOps::RVec<int>& signal_tau_idx)
{
    return lep_p4.at(signal_tau_idx.at(0)) + lep_p4.at(signal_tau_idx.at(1));
}

float getHTTScalarPt (const ROOT::VecOps::RVec<LorentzVectorM>& lep_p4, const ROOT::VecOps::RVec<int>& signal_tau_idx)
{
    return lep_p4.at(signal_tau_idx.at(0)).Pt() + lep_p4.at(signal_tau_idx.at(1)).Pt();
}

ROOT::VecOps::RVec<int> getSignalTauIndices_Gen(const ROOT::VecOps::RVec<LorentzVectorM>& lep_p4,
                                                const ROOT::VecOps::RVec<int>& lep_genTauIndex)
{
    ROOT::VecOps::RVec<int> indexes;
    int n_tau = 0;
    for(size_t n = 0; n < lep_p4.size(); ++n) {
        if(lep_genTauIndex.at(n) >= 0) {
            indexes.push_back(n);
            n_tau++;
        }
    }
    if(n_tau != 2)
        throw std::runtime_error("too few taus");
    return indexes;
}

ROOT::VecOps::RVec<int> getSignalTauIndicesDeep_Tau(const ROOT::VecOps::RVec<LorentzVectorM>& lep_p4,
                                                    const ROOT::VecOps::RVec<float>& byDeepTau2017v2p1VSjet)
{
    ROOT::VecOps::RVec<size_t> ordered_index;
    ROOT::VecOps::RVec<size_t> selected_taus_indices;
    for(size_t lep_index = 0; lep_index < lep_p4.size(); ++lep_index)
        ordered_index.push_back(lep_index);

    std::sort(ordered_index.begin(), ordered_index.end(), [&](size_t a, size_t b){
        return byDeepTau2017v2p1VSjet.at(a) > byDeepTau2017v2p1VSjet.at(b);
    });

    for(size_t n = 0; n < 2; ++n)
        selected_taus_indices.push_back(ordered_index.at(n));

    return selected_taus_indices;

}


float HTTScalarPt (const ROOT::VecOps::RVec<LorentzVectorM>& lep_p4, const ROOT::VecOps::RVec<int>& lep_genTauIndex)
{
    size_t n_tau = 0;
    float h_tautau_pt_scalar = 0.f;
    for(size_t n = 0; n < lep_p4.size(); ++n) {
        if(lep_genTauIndex.at(n) >= 0) {
            h_tautau_pt_scalar += lep_p4.at(n).Pt();
            n_tau++;
        }
    }
    if(n_tau != 2)
        throw std::runtime_error("too few taus");

    return h_tautau_pt_scalar;
}

ROOT::VecOps::RVec<float> MakeDeepFlavour_bVSall (const ROOT::VecOps::RVec<float>& jets_deepFlavour_b,
                                                  const ROOT::VecOps::RVec<float>& jets_deepFlavour_bb,
                                                  const ROOT::VecOps::RVec<float>& jets_deepFlavour_lepb)
{
    ROOT::VecOps::RVec<float> b_vs_all(jets_deepFlavour_b.size());
    for(size_t n = 0; n < jets_deepFlavour_b.size(); ++n)
        b_vs_all.at(n) = jets_deepFlavour_b.at(n) + jets_deepFlavour_bb.at(n) + jets_deepFlavour_lepb.at(n);
    return b_vs_all;
}

float jets_deepFlavour (const ROOT::VecOps::RVec<float>& deepFlavour_bVSall,
                        const ROOT::VecOps::RVec<size_t>& df_index, const size_t n)
{
    if(n < df_index.size()){
        size_t selected_index = df_index.at(n);
        return deepFlavour_bVSall.at(selected_index);
    }
    else
        return 0;
}

float jets_deepCSV (const ROOT::VecOps::RVec<float>& deepCsv_BvsAll,
                    const ROOT::VecOps::RVec<size_t>& df_index, const size_t n)
{
    if(n < df_index.size()){
        size_t selected_index = df_index.at(n);
        return deepCsv_BvsAll.at(selected_index);
    }
    else
        return 0;
}


float jet_p4_pt (const ROOT::VecOps::RVec<size_t>& df_index,
                 const ROOT::VecOps::RVec<LorentzVectorE>& jets_p4, const size_t n)
{   if(n < df_index.size())
        return jets_p4.at(df_index.at(n)).pt();
    else
        return 0;
}

float jet_p4_eta (const ROOT::VecOps::RVec<size_t>& df_index,
                  const ROOT::VecOps::RVec<LorentzVectorE>& jets_p4, const size_t n)
{   if(n < df_index.size())
        return jets_p4.at(df_index.at(n)).eta();
    else
        return 0;
}

float jet_p4_E (const ROOT::VecOps::RVec<size_t>& df_index,
                const ROOT::VecOps::RVec<LorentzVectorE>& jets_p4, const size_t n)
{   if(n < df_index.size())
        return jets_p4.at(df_index.at(n)).E();
    else
        return 0;
}

float jet_p4_M (const ROOT::VecOps::RVec<size_t>& df_index,
                const ROOT::VecOps::RVec<LorentzVectorE>& jets_p4, const size_t n)
{   if(n < df_index.size())
        return jets_p4.at(df_index.at(n)).M();
    else
        return 0;
}

float rel_jet_M_pt (const ROOT::VecOps::RVec<size_t>& df_index,
                    const ROOT::VecOps::RVec<LorentzVectorE>& jets_p4, const size_t n)
{   if(n < df_index.size())
        return jets_p4.at(df_index.at(n)).M() / jets_p4.at(df_index.at(n)).Pt();
    else
        return 0;
}

float rel_jet_E_pt (const ROOT::VecOps::RVec<size_t>& df_index,
                    const ROOT::VecOps::RVec<LorentzVectorE>& jets_p4, const size_t n)
{   if(n < df_index.size())
        return jets_p4.at(df_index.at(n)).E() / jets_p4.at(df_index.at(n)).Pt();
    else
        return 0;
}


int jet_genbJet (const ROOT::VecOps::RVec<int>& jet_genJetIndex,
                 const ROOT::VecOps::RVec<size_t>& ordered_jet_indexes, const size_t n,
                 const ROOT::VecOps::RVec<LorentzVectorE>& jets_p4)
{
    if(n < ordered_jet_indexes.size()){
        size_t index = ordered_jet_indexes.at(n);
        return jet_genJetIndex.at(index) >= 0;
    }
    else
        return 0;
}

ROOT::VecOps::RVec<int> MakeGenbJet (const ROOT::VecOps::RVec<int>& jet_genJetIndex,
                                     const ROOT::VecOps::RVec<size_t>& ordered_jet_indexes)
{
    ROOT::VecOps::RVec<int> jets_genbJet(ordered_jet_indexes.size());
    for(size_t n = 0; n < ordered_jet_indexes.size(); ++n){
        const size_t index = ordered_jet_indexes.at(n);
        jets_genbJet.at(n) = jet_genJetIndex.at(index) >= 0;
    }
    return jets_genbJet;
}

ROOT::VecOps::RVec<size_t> CreateOrderedIndex (const ROOT::VecOps::RVec<LorentzVectorE>& jets_p4,
                                               const ROOT::VecOps::RVec<float>& jets_deepFlavour, bool apply_acceptance,
                                               const ROOT::VecOps::RVec<int>& tau_indeces,
                                               const ROOT::VecOps::RVec<LorentzVectorM>& lep_p4,
                                               size_t max_jet = std::numeric_limits<size_t>::max())
{
    ROOT::VecOps::RVec<size_t> ordered_index;
    for(size_t jet_index = 0; jet_index < jets_p4.size() && ordered_index.size() < max_jet; ++jet_index) {
        if(apply_acceptance && (jets_p4.at(jet_index).pt() < 20 || abs(jets_p4.at(jet_index).eta()) > 2.4)) continue;
        bool has_overlap_with_taus = false;
        for(size_t lep_index = 0; lep_index < tau_indeces.size() && !has_overlap_with_taus; ++lep_index) {
            if(ROOT::Math::VectorUtil::DeltaR(lep_p4.at(tau_indeces.at(lep_index)), jets_p4.at(jet_index)) < 0.5)
                has_overlap_with_taus = true;
        }
        if(!has_overlap_with_taus)
            ordered_index.push_back(jet_index);
    }

    std::sort(ordered_index.begin(), ordered_index.end(), [&](size_t a, size_t b){
        return jets_deepFlavour.at(a) > jets_deepFlavour.at(b);
    });

    return ordered_index;
}

float httDeltaPhi_jet (const LorentzVectorM& htt_p4,const ROOT::VecOps::RVec<size_t>& df_index,
                       const ROOT::VecOps::RVec<LorentzVectorE>& jets_p4, const size_t n)
{
    if(n < df_index.size())
        return ROOT::Math::VectorUtil::DeltaPhi(htt_p4, jets_p4.at(df_index.at(n)));
    else
        return 0;
}

float httDeltaEta_jet (const LorentzVectorM& htt_p4, const ROOT::VecOps::RVec<size_t>& df_index,
                       const ROOT::VecOps::RVec<LorentzVectorE>& jets_p4, size_t n)
{
    if(n < df_index.size())
        return (htt_p4.eta() - jets_p4.at(df_index.at(n)).eta());
    else
        return 0;
}
