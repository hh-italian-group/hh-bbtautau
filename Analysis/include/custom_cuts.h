/*! Cuts which are customized for HHbbtautau analysis
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

namespace cuts {
const double DeltaR_MC_Match = 0.3; // <

namespace massWindow{
    const double m_tautau_low = 90;
    const double m_tautau_high = 150;
    const double m_bb_low = 70;
    const double m_bb_high = 150;
}

namespace WjetsBackgroundEstimation {
    const double HighMtRegion = 70; // > For W-jets data driven estimation
    const double HighMtRegion_low = 60; // > For W-jets data driven estimation in 2jet2tag for ltau channels
    const double HighMtRegion_high = 120; // < For W-jets data driven estimation in 2jet2tag for ltau channels
}

namespace IsolationRegionForLeptonicChannel {
    const double pfRelIso = 0.1;
    const double isolation_low = 0.2; // > For QCD data driven estimation in 2jet*tag for ltau channels
    const double isolation_high = 0.5; // < For QCD data driven estimation in 2jet*tag for ltau channels
}

namespace jetCorrections {
    const double energyUncertainty = 0.05;
}

namespace skim {
    namespace TauTau {
        namespace tauID {
            const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 10; // <
        }
    }

    namespace ETau {
        const double pFRelIso = 0.5; // <
        namespace tauID {
            const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 10; // <
        }
    }

    namespace MuTau {
        const double pFRelIso = 0.5; // <
        namespace tauID {
            const double againstMuonLoose = 0.5; // >
            const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 10; // <
        }
    }
} // namespace skim
} // namespace cuts
