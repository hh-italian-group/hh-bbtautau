/*! Definition of wrappers for CMSSW KinFitter code.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "h-tautau/Analysis/include/KinFit.h"

#include "FWCore/Utilities/src/EDMException.cc"
#include "FWCore/Utilities/src/Exception.cc"

#include "FWCore/MessageLogger/src/MessageLogger.cc"
#include "FWCore/MessageLogger/src/MessageLoggerQ.cc"
#include "FWCore/MessageLogger/src/MessageDrop.cc"
#include "FWCore/MessageLogger/src/ELseverityLevel.cc"
#include "FWCore/MessageLogger/src/MessageSender.cc"
#include "FWCore/MessageLogger/src/AbstractMLscribe.cc"
#include "FWCore/MessageLogger/src/ErrorObj.cc"
#include "FWCore/MessageLogger/src/ELextendedID.cc"
#include "FWCore/MessageLogger/src/ELstring.cc"

using std::setiosflags;
#include "PhysicsTools/KinFitter/src/TAbsFitConstraint.cc"
#include "PhysicsTools/KinFitter/src/TAbsFitParticle.cc"
#include "PhysicsTools/KinFitter/src/TFitConstraintEp.cc"
#include "PhysicsTools/KinFitter/src/TFitConstraintM.cc"
#include "PhysicsTools/KinFitter/src/TFitConstraintMGaus.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleCart.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleECart.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleEMomDev.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleEScaledMomDev.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleESpher.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleEtEtaPhi.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleEtThetaPhi.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleMCCart.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleMCMomDev.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleMCPInvSpher.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleMCSpher.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleMomDev.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleSpher.cc"
#include "PhysicsTools/KinFitter/src/TKinFitter.cc"
#include "PhysicsTools/KinFitter/src/TSLToyGen.cc"

namespace analysis {
namespace kinematic_fit {
namespace two_body {

namespace detail {
//cca from https://github.com/cvernier/kinfit/blob/master/test/kinFit4b.C

inline Double_t ErrEt(Float_t Et, Float_t Eta)
{
Double_t InvPerr2, a, b, c;
const int NYBINS = 5;
//double YBND[NYBINS+1] = {0,0.5,1.0,1.5,2.0,2.5};
double PAR[NYBINS][3] = {
{3.51, 0.826, 0.0364},{3.23, 0.851, 0.0367},{4.36, 0.871, 0.0415},
{5.22, 0.713, 0.0229},{5.07, 0.610, 0.0207}};
int iy =0;
if (fabs(Eta) < 2.5){
if (fabs(Eta) < 0.5 && fabs(Eta) >0.) iy =0;
if (fabs(Eta) < 1. && fabs(Eta) >0.5) iy =1;
if (fabs(Eta) < 1.5 && fabs(Eta) >1.) iy =2;
if (fabs(Eta) < 2. && fabs(Eta) >1.5) iy =3;
if (fabs(Eta) < 2. && fabs(Eta) >2.5) iy =4;
///std::cout<< iy << " kin iy "<<std::endl;
InvPerr2 = sqrt(pow(PAR[iy][0]/Et,2) + pow(PAR[iy][1],2)/Et + pow(PAR[iy][2],2));
//std::cout<< InvPerr2 << " kin invPerr "<<std::endl;
}
else {
a = 4.8;
b = 0.89;
c = 0.043;
InvPerr2 = (a * a) + (b * b) * Et + (c * c) * Et * Et;
}
return InvPerr2;
}
//cca
inline Double_t ErrEta(Float_t Et, Float_t Eta)
{
    Double_t InvPerr2, a, b, c;
    if(fabs(Eta) < 1.4){
        a = 1.215;
        b = 0.037;
        c = 7.941 * 0.0001;
    } else {
        a = 1.773;
        b = 0.034;
        c = 3.56 * 0.0001;
    }
    InvPerr2 = a/(Et * Et) + b/Et + c;
    return InvPerr2;
}
/*cca
inline Double_t ErrEt(Float_t Et, Float_t Eta)
{
    Double_t InvPerr2, a, b, c;
      if(fabs(Eta) < 1.4){
        a = 5.6;
    b = 1.25;
        c = 0.033;
      }else{
        a = 4.8;
        b = 0.89;
        c = 0.043;
    }
    InvPerr2 = (a * a) + (b * b) * Et + (c * c) * Et * Et;
    return InvPerr2;
}
*/

inline Double_t ErrPhi(Float_t Et, Float_t Eta)
{
    Double_t InvPerr2, a, b, c;
      if(fabs(Eta) < 1.4){
        a = 6.65;
        b = 0.04;
        c = 8.49 * 0.00001;
      }else{
          a = 2.908;
          b = 0.021;
          c = 2.59 * 0.0001;
      }
    InvPerr2 = a/(Et * Et) + b/Et + c;
      return InvPerr2;
}

} // namespace detail

inline FitResults Fit_KinFitter(const FitInput& input)
{
    FitResults result;
    TMatrixD m1(3, 3);
    TMatrixD m2(3, 3);
    m1.Zero();
    m2.Zero();

    TLorentzVector mom1(input.bjet_momentums.at(0)), mom2(input.bjet_momentums.at(1));

    m1(0, 0) = detail::ErrEt(mom1.Et(), mom1.Eta());
    m1(1, 1) = detail::ErrEta(mom1.Et(), mom1.Eta());
    m1(2, 2) = detail::ErrPhi(mom1.Et(), mom1.Eta());

    m2(0, 0) = detail::ErrEt(mom2.Et(), mom2.Eta());
    m2(1, 1) = detail::ErrEta(mom2.Et(), mom2.Eta());
    m2(2, 2) = detail::ErrPhi(mom2.Et(), mom2.Eta());


    TFitParticleEtEtaPhi _jet1(&mom1, &m1);
    TFitParticleEtEtaPhi _jet2(&mom2, &m2);

    TFitConstraintM m_bb;
    m_bb.addParticle1(&_jet1);
    m_bb.addParticle1(&_jet2);
    m_bb.setMassConstraint(125.6);
//cca Higgs mass 125.6 and not 125.0

    TKinFitter _fitter;
    _fitter.addMeasParticle(&_jet1);
    _fitter.addMeasParticle(&_jet2);
    _fitter.addConstraint(&m_bb);

    _fitter.setMaxNbIter(30);
    _fitter.setMaxDeltaS(1e-2);
    _fitter.setMaxF(1e-1);
    _fitter.setVerbosity(0);
    const Int_t fit_result = _fitter.fit();

/*    std::cout << "test fit result: " << fit_result << std::endl;
    if(fit_result == 0) {
        std::cout << "Old 1: " << input.bjet_momentums.at(0) << std::endl;
        std::cout << "Old 2: " << input.bjet_momentums.at(1) << std::endl;
        std::cout << "Old mass 1: " << input.bjet_momentums.at(0).M() << std::endl;
        std::cout << "Old mass 2: " << input.bjet_momentums.at(1).M() << std::endl;
        std::cout << "New 1: " << *_jet1.getCurr4Vec() << std::endl;
        std::cout << "New 2: " << *_jet2.getCurr4Vec() << std::endl;
        std::cout << "New mass 1: " << _jet1.getCurr4Vec()->M() << std::endl;
        std::cout << "New mass 2: " << _jet2.getCurr4Vec()->M() << std::endl;
        const double m_old = (input.bjet_momentums.at(0) + input.bjet_momentums.at(1)).M();
        const double m_new = ((*_jet1.getCurr4Vec()) + (*_jet2.getCurr4Vec())).M();
        std::cout << "M old: " << m_old << std::endl;
        std::cout << "M new: " << m_new << std::endl;
   }
*/
    result.convergence = fit_result;
    result.bjet_momentums.push_back(*_jet1.getCurr4Vec());
    result.bjet_momentums.push_back(*_jet2.getCurr4Vec());
    return result;
}

} // namespace two_body
} // namespace kinematic_fit
} // namespace analysis
