import ROOT
from getDataIdHash import LoadIdFrames, 

def select(df):
    df = df.Filter('''genLepton_index >= 0 && genLepton_kind == 5 && genLepton_lastMotherIndex == 0
                      && (genParticle_pdgId.at(genLepton_lastMotherIndex) == 9900012
                      || genParticle_pdgId.at(genLepton_lastMotherIndex) == 9990012)''')
    df = df.Define('genLepton', '''reco_tau::gen_truth::GenLepton::fromRootTuple(
        genLepton_lastMotherIndex, genParticle_pdgId, genParticle_mother, genParticle_charge,
        genParticle_isFirstCopy, genParticle_isLastCopy,
        genParticle_pt, genParticle_eta, genParticle_phi, genParticle_mass,
        genParticle_vtx_x, genParticle_vtx_y, genParticle_vtx_z)
    ''')
    # df = df.Filter('''genLepton.kind() == reco_tau::gen_truth::GenLepton::Kind::TauDecayedToHadrons
    #                && genLepton.mothers().size() == 1 && (*genLepton.mothers().begin())->pdgId == 9900012
    # ''')
    df = df.Define('genLepton_mother_prod_vtx_x', '(*genLepton.mothers().begin())->vertex.x()')
    df = df.Define('genLepton_mother_prod_vtx_y', '(*genLepton.mothers().begin())->vertex.y()')
    df = df.Define('genLepton_mother_prod_vtx_z', '(*genLepton.mothers().begin())->vertex.z()')
    df = df.Define('genLepton_prod_vtx_x', 'genLepton.lastCopy().vertex.x()')
    df = df.Define('genLepton_prod_vtx_y', 'genLepton.lastCopy().vertex.y()')
    df = df.Define('genLepton_prod_vtx_z', 'genLepton.lastCopy().vertex.z()')
    df = df.Define('genLepton_decay_vtx_x', '(*genLepton.lastCopy().daughters.begin())->vertex.x()')
    df = df.Define('genLepton_decay_vtx_y', '(*genLepton.lastCopy().daughters.begin())->vertex.y()')
    df = df.Define('genLepton_decay_vtx_z', '(*genLepton.lastCopy().daughters.begin())->vertex.z()')
    df = df.Define('genLepton_nChargedHadrons', 'static_cast<int>(genLepton.nChargedHadrons())')
    df = df.Define('genLepton_nNeutralHadrons', 'static_cast<int>(genLepton.nNeutralHadrons())')
    return df
