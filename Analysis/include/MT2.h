#include "Lester_mt2_bisect.h"

template<typename LVector1, typename LVector2, typename LVector3, typename LVector4, typename LVector5 >
double Calculate_MT2(const LVector1& lepton1_p4, const LVector2& lepton2_p4, const LVector3& bjet_1, const LVector4& bjet_2, const LVector5& met_p4)
{
    asymm_mt2_lester_bisect::disableCopyrightMessage();
    const double mVisA = (lepton1_p4 + bjet_1).mass();
    const double pxA = (lepton1_p4 + bjet_1).px();
    const double pyA = (lepton1_p4 + bjet_1).py();
    const double mVisB = (lepton2_p4 + bjet_2).mass();
    const double pxB = (lepton2_p4 + bjet_2).px();
    const double pyB = (lepton2_p4 + bjet_2).py();
    const double pxMet = met_p4.px();
    const double pyMet = met_p4.py();
    double chiA = 0.; // hypothesised mass of invisible on side A.  Must be >=0.
    double chiB = 0.; // hypothesised mass of invisible on side B.  Must be >=0.
    double MT2 =  asymm_mt2_lester_bisect::get_mT2(mVisA, pxA, pyA,mVisB, pxB, pyB,pxMet, pyMet,chiA, chiB,0);
    return MT2;
}
