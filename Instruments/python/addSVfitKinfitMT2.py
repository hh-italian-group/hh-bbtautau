import ROOT

ROOT.gROOT.ProcessLine(".include .")

ROOT.gROOT.ProcessLine('#include "h-tautau/Instruments/src/Lester_mt2_bisect.cpp"')
ROOT.gROOT.ProcessLine('#include "HHKinFit2/HHKinFit2/interface/HHKinFitMasterHeavyHiggs.h"')
ROOT.gROOT.ProcessLine('#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"')
ROOT.gROOT.ProcessLine('#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"')
ROOT.gROOT.ProcessLine('#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"')

ROOT.gSystem.Load('../lib/slc7_amd64_gcc820/libHHKinFit2HHKinFit2.so')
ROOT.gSystem.Load('../lib/slc7_amd64_gcc820/libTauAnalysisClassicSVfit.so')


ROOT.gInterpreter.Declare('''

using LorentzVectorXYZ = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>;
using LorentzVectorM = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;

template<typename LVector1, typename LVector2, typename LVector3, typename LVector4, typename LVector5 >
double Calculate_MT2(const LVector1& lepton1_p4, const LVector2& lepton2_p4, const LVector3& bjet_1, const LVector4& bjet_2, const LVector5& met_p4)
{
    asymm_mt2_lester_bisect::disableCopyrightMessage();
    const double mVisA = bjet_1.mass();
    const double pxA = bjet_1.px();
    const double pyA = bjet_1.py();
    const double mVisB = bjet_2.mass();
    const double pxB = bjet_2.px();
    const double pyB = bjet_2.py();
    const double pxMiss = lepton1_p4.px() + lepton2_p4.px() + met_p4.px();
    const double pyMiss = lepton1_p4.py() + lepton2_p4.py() + met_p4.py();
    double chiA = lepton1_p4.mass(); // hypothesised mass of invisible on side A.  Must be >=0.
    double chiB = lepton2_p4.mass(); // hypothesised mass of invisible on side B.  Must be >=0.
    double MT2 =  asymm_mt2_lester_bisect::get_mT2(mVisA, pxA, pyA,mVisB, pxB, pyB,pxMiss, pyMiss,chiA, chiB,0);
    return MT2;
}

template<typename LVector>
TLorentzVector ConvertVector(const LVector& v)
{
    return TLorentzVector(v.Px(), v.Py(), v.Pz(), v.E());
}

namespace kin_fit {

struct FitResults {
    double mass, chi2, probability;
    int convergence;
    bool HasValidMass() const { return convergence > 0; }

    FitResults() : convergence(std::numeric_limits<int>::lowest()) {}
    FitResults(double _mass, double _chi2, double _probability, int _convergence) :
        mass(_mass), chi2(_chi2), probability(_probability), convergence(_convergence) {}
};

class FitProducer {
public:
    template<typename LVector1, typename LVector2, typename LVector3, typename LVector4, typename LVector5>
    static FitResults Fit(const LVector1& lepton1_p4, const LVector2& lepton2_p4,
                          const LVector3& jet1_p4, const LVector4& jet2_p4,
                          const LVector5& met_p4, double MET_cov_00, double MET_cov_01, double MET_cov_11,
                          double resolution_1, double resolution_2, int verbosity = 0)
    {
        TMatrixD met_cov(2, 2);
        met_cov(0, 0) = MET_cov_00;
        met_cov(0, 1) = MET_cov_01;
        met_cov(1, 0) = MET_cov_01;
        met_cov(1, 1) = MET_cov_11;

        return FitImpl(ConvertVector(lepton1_p4), ConvertVector(lepton2_p4), ConvertVector(jet1_p4),
                       ConvertVector(jet2_p4), TVector2(met_p4.Px(), met_p4.Py()), met_cov,
                       resolution_1, resolution_2, verbosity);
    }

    static FitResults FitImpl(const TLorentzVector& lepton1_p4, const TLorentzVector& lepton2_p4,
                              const TLorentzVector& jet1_p4, const TLorentzVector& jet2_p4, const TVector2& met,
                              const TMatrixD& met_cov, double resolution_1, double resolution_2, int verbosity)
    {
      FitResults result;
      try {

          HHKinFit2::HHKinFitMasterHeavyHiggs hh_kin_fit(jet1_p4, jet2_p4, lepton1_p4, lepton2_p4, met, met_cov,
                                                         resolution_1, resolution_2);
          hh_kin_fit.verbosity = verbosity;
          hh_kin_fit.fit();

          result.convergence = hh_kin_fit.getConvergence();
          if(result.HasValidMass()) {
              result.mass = hh_kin_fit.getMH();
              result.chi2 = hh_kin_fit.getChi2();
              result.probability = hh_kin_fit.getFitProb();
          }

          if(verbosity > 0) {
              std::cout << "Convergence = " << result.convergence << std::endl;
              if(result.HasValidMass()) {
                  std::cout << std::setprecision(6);
                  std::cout << "Mass = " << result.mass << std::endl;
                  std::cout << "chi2 = " << result.chi2 << std::endl;
                  std::cout << "probability = " << result.probability << std::endl;
              }
          }
      } catch(std::exception&) {}

      return result;
    }
};
}

enum class LegType { e = 0, mu = 1, tau = 2, jet = 3 };

namespace sv_fit {

struct FitResults {
    bool has_valid_momentum;
    LorentzVectorM momentum;
    LorentzVectorM momentum_error;
    double transverseMass;
    double transverseMass_error;

    FitResults() :
        has_valid_momentum(false), transverseMass(std::numeric_limits<double>::lowest()),
        transverseMass_error(std::numeric_limits<double>::lowest()) {}
    FitResults(bool _has_valid_momentum, LorentzVectorM _momentum, LorentzVectorM _momentum_error,
               double _transverseMass, double _transverseMass_error) :
        has_valid_momentum(_has_valid_momentum), momentum(_momentum), momentum_error(_momentum_error),
        transverseMass(_transverseMass), transverseMass_error(_transverseMass_error) {}
};

class FitProducer {
public:
    static classic_svFit::MeasuredTauLepton CreateMeasuredLepton(const LorentzVectorM& momentum, LegType leg_type,
                                                                 int decay_mode)
    {
        double preciseVisMass = momentum.mass();
        classic_svFit::MeasuredTauLepton::kDecayType decay_type;
        if(leg_type == LegType::e) {
            static const double minVisMass = classic_svFit::electronMass, maxVisMass = minVisMass;
            preciseVisMass = std::clamp(preciseVisMass, minVisMass, maxVisMass);
            decay_type = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
            decay_mode = -1;
        } else if(leg_type == LegType::mu) {
            decay_type = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
            decay_mode = -1;
        } else if(leg_type == LegType::tau){
            const double minVisMass = decay_mode == 0 ? classic_svFit::chargedPionMass : 0.3;
            const double maxVisMass = decay_mode == 0 ? classic_svFit::chargedPionMass : 1.5;
            preciseVisMass = std::clamp(preciseVisMass, minVisMass, maxVisMass);
            decay_type = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
        } else {
            throw std::runtime_error("Leg Type not supported for SVFitAnaInterface.");
        }

        return classic_svFit::MeasuredTauLepton(decay_type, momentum.Pt(), momentum.Eta(), momentum.Phi(),
                                                preciseVisMass, decay_mode);
    }

    static FitResults Fit(const LorentzVectorM& tau1_p4, LegType tau1_leg_type, int tau1_decay_mode,
                          const LorentzVectorM& tau2_p4, LegType tau2_leg_type, int tau2_decay_mode,
                          const LorentzVectorM& met_p4, double MET_cov_00, double MET_cov_01, double MET_cov_11,
                          int verbosity = 0)
    {
        static const auto init = []() { TH1::AddDirectory(false); return true; };
        static const bool initialized = init();
        (void) initialized;

        const std::vector<classic_svFit::MeasuredTauLepton> measured_leptons = {
            CreateMeasuredLepton(tau1_p4, tau1_leg_type, tau1_decay_mode),
            CreateMeasuredLepton(tau2_p4, tau2_leg_type, tau2_decay_mode)
        };

        TMatrixD met_cov_t(2, 2);
        met_cov_t(0, 0) = MET_cov_00;
        met_cov_t(0, 1) = MET_cov_01;
        met_cov_t(1, 0) = MET_cov_01;
        met_cov_t(1, 1) = MET_cov_11;

        ClassicSVfit algo(verbosity);
        algo.addLogM_fixed(false);
        algo.addLogM_dynamic(false);
        // algo.setDiTauMassConstraint(-1.0);
        algo.integrate(measured_leptons, met_p4.Px(), met_p4.Py(), met_cov_t);

        FitResults result;
        if(algo.isValidSolution()) {
            auto histoAdapter = dynamic_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter());
            result.momentum = LorentzVectorM(histoAdapter->getPt(), histoAdapter->getEta(), histoAdapter->getPhi(), histoAdapter->getMass());
            result.momentum_error = LorentzVectorM(histoAdapter->getPtErr(), histoAdapter->getEtaErr(), histoAdapter->getPhiErr(), histoAdapter->getMassErr());
            result.transverseMass = histoAdapter->getTransverseMass();
            result.transverseMass_error = histoAdapter->getTransverseMassErr();
            result.has_valid_momentum = true;
        }
        return result;
    }
};

} // namespace sv_fit_ana


''')

columns_to_exclude = [ 'tau1_p4', 'tau2_p4', 'b1_p4', 'b2_p4', 'MET_p4', 'kinFit_result', 'SVfit_result' ]


def add(df):
    for leg in [ 'tau', 'b' ]:
        for idx in range(1, 3):
            df = df.Define('{}{}_p4'.format(leg, idx),
                           'LorentzVectorM({0}{1}_pt, {0}{1}_eta, {0}{1}_phi, {0}{1}_m)'.format(leg, idx))
    df = df.Define('MET_p4', 'LorentzVectorM(MET_pt, 0, MET_phi, 0)')

    df = df.Define('MT2', 'float(Calculate_MT2(tau1_p4, tau2_p4, b1_p4, b2_p4, MET_p4))')

    df = df.Define('kinFit_result',
                   '''kin_fit::FitProducer::Fit(tau1_p4, tau2_p4, b1_p4, b2_p4, MET_p4, MET_cov_00, MET_cov_01,
                                                MET_cov_11, b1_resolution, b2_resolution)''')
    df = df.Define('kinFit_convergence', 'kinFit_result.convergence')
    df = df.Define('kinFit_m', 'float(kinFit_result.mass)')
    df = df.Define('kinFit_chi2', 'float(kinFit_result.chi2)')

    df = df.Define('SVfit_result',
                   '''sv_fit::FitProducer::Fit(tau1_p4, LegType::tau, tau1_decay_mode,
                                               tau2_p4, LegType::tau, tau2_decay_mode,
                                               MET_p4, MET_cov_00, MET_cov_01, MET_cov_11)''')

    df = df.Define('SVfit_valid', 'int(SVfit_result.has_valid_momentum)')
    df = df.Define('SVfit_pt', 'float(SVfit_result.momentum.pt())')
    df = df.Define('SVfit_eta', 'float(SVfit_result.momentum.eta())')
    df = df.Define('SVfit_phi', 'float(SVfit_result.momentum.phi())')
    df = df.Define('SVfit_m', 'float(SVfit_result.momentum.mass())')
    df = df.Define('SVfit_pt_error', 'float(SVfit_result.momentum_error.pt())')
    df = df.Define('SVfit_eta_error', 'float(SVfit_result.momentum_error.eta())')
    df = df.Define('SVfit_phi_error', 'float(SVfit_result.momentum_error.phi())')
    df = df.Define('SVfit_m_error', 'float(SVfit_result.momentum_error.mass())')
    df = df.Define('SVfit_mt', 'float(SVfit_result.transverseMass)')
    df = df.Define('SVfit_mt_error', 'float(SVfit_result.transverseMass_error)')
    return df
