/*! Definiton of analysis::FlatEventInfo class.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/GenParticle.h"

namespace analysis {

GenParticle::GenParticle(const ntuple::Event& events, size_t n)
    : index(n)
{
    if(n >= events.genParticles_p4.size())
        throw std::runtime_error("events index is out of range.");

    pdg = events.genParticles_pdg.at(n);
    status = events.genParticles_status.at(n);
    genStatusFlags = GenStatusFlags(events.genParticles_statusFlags.at(n));
    momentum = events.genParticles_p4.at(n);
}


GenEvent::GenEvent(const ntuple::Event& event)
{
    genParticles.reserve(event.genParticles_p4.size());
    for(size_t n = 0; n < event.genParticles_p4.size(); ++n)
        genParticles.emplace_back(event,n);

    for (size_t relIndex = 0; relIndex < event.genParticles_rel_mIndex.size(); relIndex++) {
        int motherIndex = event.genParticles_rel_mIndex.at(relIndex);
        int daughterIndex = event.genParticles_rel_pIndex.at(relIndex);

        GenParticle* mother = &genParticles.at(static_cast<size_t>(motherIndex));
        GenParticle* daughter = &genParticles.at(static_cast<size_t>(daughterIndex));

        mother->daughters.push_back(daughter);
        daughter->mothers.push_back(mother);
    }

    for(const GenParticle& genParticle : genParticles ){
        if (!genParticle.mothers.size())
            primaryParticles.insert(&genParticle);


         if (genParticle.genStatusFlags.isPrompt())
            particleCodeMap[std::abs(genParticle.pdg)].insert(&genParticle);
    }
}

GenParticleSet GenEvent::GetParticles(int particle_pgd) const
{
    GenParticleSet results;
    const ParticleCodeMap::const_iterator code_iter = particleCodeMap.find(particle_pgd);

    if (code_iter != particleCodeMap.end()){
        for (const GenParticle* particle : code_iter->second){
            if(!particle->genStatusFlags.isLastCopy()) continue;
            results.insert(particle);
        }
    }
    return results;
}

std::vector<const GenParticle*> GenEvent::GetTypesParticles(std::vector<std::string> type_names, const GenParticle* possible_mother ) const
{
    std::vector<int> selected_particles_type_pgd;
    for(auto pdg_type : particle_types){
        for(size_t type_index = 0; type_index < type_names.size(); ++type_index){
            if(pdg_type.second ==  type_names.at(type_index))
                selected_particles_type_pgd.push_back(pdg_type.first);
        }
    }
    std::vector<const GenParticle*> results;
    results.clear();
    for (const auto& code : selected_particles_type_pgd){
        const ParticleCodeMap::const_iterator code_iter = particleCodeMap.find(code);
        if (code_iter == particleCodeMap.end())
            continue;
        for (const GenParticle* particle : code_iter->second){
           if(!particle->genStatusFlags.isLastCopy()) continue;
           if(!particle->genStatusFlags.isPrompt()) continue;
           if(!areParented(particle,possible_mother)) continue;
            //to check for shared particles in the daughters of the possible mother
//               if(!areParented(particle,possible_mother->daughters.at(0)))continue;
//               if(!areParented(particle,possible_mother->daughters.at(1)))continue;
           results.push_back(particle);
        }
    }
    return results;
}

bool GenEvent::areParented(const GenParticle* daughter, const GenParticle* possible_mother) const
{
    for (size_t i = 0; i < daughter->mothers.size(); i++) {
        const GenParticle* mother = daughter->mothers.at(i);
        if(mother == possible_mother)
            return true;
        if(areParented(mother,possible_mother))
            return true;
    }
    return false;
}

static const std::string& GenEvent::GetParticleName(int pdgId)
{
    auto iter = particle_names.find(pdgId);
    if(iter == particle_names.end()) throw exception("Name not found for particle with pdgId = %1%") % pdgId;
    return iter->second;
}

static void GenEvent::intializeNames(const std::string& fileName)
{
    particle_names.clear();
    particle_types.clear();
    std::ifstream f(fileName.c_str());
    if(!f.is_open() )
        throw analysis::exception("Unable to read the configuration file ");

    while(f.good()) {
        std::string line;
        std::getline(f, line);
        if(!line.length() || line[0] == '#')
            continue;
        auto pdg_name_type = SplitValueList(line, false, ",");
        if(pdg_name_type.size() < 2)
            throw exception("Invalid particle definition (1): '%1%'") % line;

        const int pdgId = Parse<int>(pdg_name_type.at(0));
        const std::string& name = pdg_name_type.at(1);
        if(particle_names.count(pdgId))
            throw exception("Duplicated definition of particle with pdgId = %1%") % pdgId;
        particle_names[pdgId] = name;

        if(pdg_name_type.size() == 3){
            const std::string& type = pdg_name_type.at(2);
            if(particle_types.count(pdgId))
                throw exception("Duplicated definition of particle with pdgId = %1%") % pdgId;
            particle_types[pdgId] = type;
        }
    }
}

void GenEvent::Print() const
{
    for (const GenParticle* particle : primaryParticles) {
        PrintChain(particle, "");
    }
}

void GenEvent::PrintChain(const GenParticle* particle, const std::string& pre, unsigned iteration = 0) const
{
    const int pdgParticle = particle->pdg;
    const auto particleName = GetParticleName(pdgParticle);
    const int particleStatus = particle->status;
    const LorentzVectorM_Float genParticle_momentum = particle->momentum;
    const std::bitset<15> genStatusFlags_ = particle->genStatusFlags.flags ;
    for (unsigned n = 0; n < iteration; ++n)
        std::cout << "  ";
    auto mother_index = particle->mothers.size() > 0 ?  particle->mothers.at(0)->index : 0;
    std::cout << particleName                              << " <" << pdgParticle
              << "> pt=" << genParticle_momentum.Pt()      << " eta=" << genParticle_momentum.Eta()
              << " phi=" << genParticle_momentum.Phi()     << " E=" << genParticle_momentum.E()
              << " m=" << genParticle_momentum.M()         << " index=" << particle->index
              << " mother_index=" << mother_index          << " status=" << particleStatus
              << " statusFlags=" << genStatusFlags_        << std::endl;

    for(unsigned n = 0; n < particle->daughters.size(); ++n) {
        const GenParticle* daughter = particle->daughters.at(n);
        std::cout << pre << "+-> ";
        const char pre_first = particle->daughters.at(n) == particle->daughters.at(particle->daughters.size() -1) ? ' ' : '|';
        const std::string pre_d = pre + pre_first + "   ";
        PrintChain(daughter, pre_d, iteration);
    }
}

void FindFinalStateDaughters(const GenParticle& particle, std::set<const GenParticle*>& daughters,
                             const std::set<int>& pdg_to_exclude)
{
    if(!particle.daughters.size()) {
        const int abs_pdg = std::abs(particle.pdg);
        if(!pdg_to_exclude.count(abs_pdg))
            daughters.insert(&particle);
    } else {
        for(const auto& daughter : particle.daughters)
            FindFinalStateDaughters(*daughter, daughters, pdg_to_exclude);
    }
}


LorentzVectorM_Float GetFinalStateMomentum(const GenParticle& particle, std::vector<const GenParticle*>& visible_daughters,
                                   bool excludeInvisible, bool excludeLightLeptons)
{
    using pair = std::pair<bool, bool>;
    static const std::set<int> empty = {};

    static const std::map<pair, const std::set<int>*> to_exclude {
        { pair(false, false), &empty }, { pair(true, false), &particles::neutrinos },
        { pair(false, true), &particles::light_leptons }, { pair(true, true), &particles::light_and_invisible },
    };

    std::set<const GenParticle*> daughters_set;

    FindFinalStateDaughters(particle, daughters_set, *to_exclude.at(pair(excludeInvisible, false)));
    visible_daughters.clear();
    visible_daughters.insert(visible_daughters.begin(), daughters_set.begin(), daughters_set.end());

    LorentzVectorM_Float p4;
    for(auto daughter : visible_daughters) {
        if(excludeLightLeptons && particles::light_leptons.count(std::abs(daughter->pdg))
            && daughter->genStatusFlags.isDirectTauDecayProduct()) continue;
        p4 += daughter->momentum;
    }
    return p4;
}
} //analysis
