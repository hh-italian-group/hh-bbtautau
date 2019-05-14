#pragma once

#include <cmath>
#include <string>
#include <map>
#include <set>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <iostream>
#include <limits>
#include "AnalysisTools/Core/include/Tools.h"

namespace particles {
class ParticleCode {

public:
    enum pdg {e= 11, mu =13, tau = 15 , b = 5, higgs = 25, nu_e = 12, nu_mu = 14, nu_tau = 16};

};

static const std::set<int> neutrinos = { ParticleCode::nu_e, ParticleCode::nu_mu, ParticleCode::nu_tau };

static const std::set<int> light_leptons = { ParticleCode::e, ParticleCode::mu };

static const std::set<int> light_and_invisible = analysis::tools::union_sets({light_leptons, neutrinos});

} // particles
