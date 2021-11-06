import ROOT

ROOT.gInterpreter.Declare('''

using IdType = unsigned long long;
using LorentzVectorXYZ = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>;
using LorentzVectorM = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;

struct JetEntry {
    float eta, phi, bRegCorr, bRegRes;
};

static std::map<ULong64_t, std::vector<JetEntry>> regCache;

void AddEventId(IdType event, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_phi,
                const ROOT::VecOps::RVec<float>& Jet_bRegCorr, const ROOT::VecOps::RVec<float>& Jet_bRegRes)
{
    static std::mutex m;
    std::lock_guard<std::mutex> lock(m);
    for(size_t n = 0; n < Jet_eta.size(); ++n) {
        JetEntry entry;
        entry.eta = Jet_eta.at(n);
        entry.phi = Jet_phi.at(n);
        entry.bRegCorr = Jet_bRegCorr.at(n);
        entry.bRegRes = Jet_bRegRes.at(n);
        regCache[event].push_back(entry);
    }
}

void LoadRegressionCache(const std::string& file_name)
{
    ROOT::RDataFrame df("Events", file_name);
    df.Foreach(AddEventId, {"event", "Jet_eta", "Jet_phi", "Jet_bRegCorr", "Jet_bRegRes"});
}

JetEntry GetRegressedInfo(IdType evt, const LorentzVectorM& orig_p4)
{
    static constexpr double dR_thr = 0.01;
    static constexpr double dR2_thr = dR_thr * dR_thr;
    const auto iter = regCache.find(evt);
    if(iter != regCache.end()) {
        const auto& jets = iter->second;
        std::set<size_t> matched;
        for(size_t n = 0; n < jets.size(); ++n) {
            const double dEta = jets.at(n).eta - orig_p4.eta();
            const double dPhi = ROOT::Math::VectorUtil::Phi_mpi_pi(jets.at(n).phi - orig_p4.phi());
            const double dR2 = std::pow(dEta, 2) + std::pow(dPhi, 2);
            if(dR2 < dR2_thr)
                matched.insert(n);
        }
        if(matched.size() == 1)
            return jets.at(*matched.begin());
        if(matched.size() > 1)
            throw std::runtime_error("Multiple matching.");
    }

    throw std::runtime_error("Jet entry not found.");
}

''')

columns_to_exclude = [ 'b1_p4_orig', 'b2_p4_orig', 'b1_reg_info', 'b2_reg_info', 'b1_p4_reg', 'b2_p4_reg',
                       'MET_px', 'MET_py', 'MET_p4' ]


def add(df, nanoFile):
    ROOT.LoadRegressionCache(nanoFile)
    for idx in range(1, 3):
        df = df.Define('b{}_p4_orig'.format(idx),
                       'LorentzVectorM(b{0}_pt_orig, b{0}_eta, b{0}_phi, b{0}_m_orig)'.format(idx))
        df = df.Define('b{}_reg_info'.format(idx), 'GetRegressedInfo(evt, b{}_p4_orig)'.format(idx))
        df = df.Define('b{}_p4_reg'.format(idx), 'b{0}_p4_orig * b{0}_reg_info.bRegCorr'.format(idx))
        df = df.Define('b{}_pt'.format(idx), 'float(b{}_p4_reg.pt())'.format(idx))
        df = df.Define('b{}_m'.format(idx), 'float(b{}_p4_reg.mass())'.format(idx))
        df = df.Define('b{}_resolution'.format(idx), 'float(b{0}_p4_reg.pt() * b{0}_reg_info.bRegRes)'.format(idx))

    df = df.Define('MET_px', '(b1_p4_reg.Px() - b1_p4_orig.Px()) + (b2_p4_reg.Px() - b2_p4_orig.Px())')
    df = df.Define('MET_py', '(b1_p4_reg.Py() - b1_p4_orig.Py()) + (b2_p4_reg.Py() - b2_p4_orig.Py())')
    df = df.Define('MET_p4', 'LorentzVectorXYZ(MET_px, MET_py, 0, std::hypot(MET_px, MET_py))')
    df = df.Define('MET_pt', 'MET_p4.pt()')
    df = df.Define('MET_phi', 'MET_p4.phi()')

    return df
