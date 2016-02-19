/*! Definition of Config class.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/SimpleConfigBase.h"

namespace analysis {

class Config : public BaseConfig {
public:
    Config(const std::string& fileName) { Read(fileName); }

    ANA_CONFIG_PARAMETER(unsigned, ReportInterval, 10)
    ANA_CONFIG_PARAMETER(bool, RunSingleEvent, false)
    ANA_CONFIG_PARAMETER(unsigned, SingleEventId, 0)
    ANA_CONFIG_PARAMETER(unsigned, MaxTreeVersion, 1)

    ANA_CONFIG_PARAMETER(bool, isMC, false)
    ANA_CONFIG_PARAMETER(bool, ApplyTauESCorrection, false)
    ANA_CONFIG_PARAMETER(bool, ApplyRecoilCorrection, false)
    ANA_CONFIG_PARAMETER(bool, ApplyRecoilCorrectionForW, false)
    ANA_CONFIG_PARAMETER(bool, ApplyRecoilCorrectionForZ, false)
    ANA_CONFIG_PARAMETER(bool, ExpectedOneNonSMResonance, false)
    ANA_CONFIG_PARAMETER(bool, ExpectedAtLeastOneSMResonanceToTauTauOrToBB, false)
    ANA_CONFIG_PARAMETER(bool, RequireSpecificFinalState, false)
    ANA_CONFIG_PARAMETER(bool, DoZEventCategorization, false)
    ANA_CONFIG_PARAMETER(bool, isDYEmbeddedSample, false)
    ANA_CONFIG_PARAMETER(bool, isTTEmbeddedSample, false)
    ANA_CONFIG_PARAMETER(bool, ApplyEtoTauFakeRate, false)
    ANA_CONFIG_PARAMETER(bool, ApplyJetToTauFakeRate, false)
    ANA_CONFIG_PARAMETER(bool, ApplyDMweight, false)

    ANA_CONFIG_PARAMETER(bool, EstimateTauEnergyUncertainties, false)
    ANA_CONFIG_PARAMETER(bool, EstimateJetEnergyUncertainties, false)
    ANA_CONFIG_PARAMETER(bool, EstimateBtagEfficiencyUncertainties, false)
    ANA_CONFIG_PARAMETER(bool, EstimateBtagFakeUncertainties, false)
    ANA_CONFIG_PARAMETER(std::string, JetEnergyUncertainties_inputFile, "")
    ANA_CONFIG_PARAMETER(std::string, JetEnergyUncertainties_inputSection, "")

    ANA_CONFIG_PARAMETER(bool, ApplyPUreweight, false)
    ANA_CONFIG_PARAMETER(std::string, PUreweight_fileName, "")
    ANA_CONFIG_PARAMETER(double, PUreweight_maxAvailablePU, 60.0)
    ANA_CONFIG_PARAMETER(double, PUreweight_defaultWeight, 0.0)

    ANA_CONFIG_PARAMETER(double, MvaMet_dZcut, 0.1)
    ANA_CONFIG_PARAMETER(std::string, MvaMet_inputFileNameU, "")
    ANA_CONFIG_PARAMETER(std::string, MvaMet_inputFileNameDPhi, "")
    ANA_CONFIG_PARAMETER(std::string, MvaMet_inputFileNameCovU1, "")
    ANA_CONFIG_PARAMETER(std::string, MvaMet_inputFileNameCovU2, "")

    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileCorrectTo_MuTau, "")
    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileZmmData_MuTau, "")
    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileZmmMC_MuTau, "")

    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileCorrectTo_ETau, "")
    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileZmmData_ETau, "")
    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileZmmMC_ETau, "")

    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileCorrectTo_TauTau, "")
    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileZmmData_TauTau, "")
    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileZmmMC_TauTau, "")

    bool extractMCtruth()
    {
        return ApplyTauESCorrection() || ApplyRecoilCorrection() || RequireSpecificFinalState()
                || ExpectedOneNonSMResonance() || ExpectedAtLeastOneSMResonanceToTauTauOrToBB()
                || DoZEventCategorization() || ApplyEtoTauFakeRate() || IsEmbeddedSample() || ApplyDMweight();
    }

    bool IsEmbeddedSample() { return isDYEmbeddedSample() || isTTEmbeddedSample(); }

};

} // analysis
