#pragma once

#include <set>
#include <TLorentzVector.h>

#include "h-tautau/Core/include/EventTuple.h"
#include "GenStatusFlags.h"
#include "hh-bbtautau/Analysis/include/Particle.h"
#include "TextIO.h"


// #include "Particle.h"
#include "exception.h"
#include <iostream>
#include <fstream>


namespace analysis {

class GenParticle;

using GenParticleVector = std::vector<GenParticle>;
using GenParticleSet = std::set<const GenParticle*>;
using GenParticlePtrVector = std::vector<const GenParticle*>;
using ParticleCodeMap = std::map<int,GenParticleSet>;
using GenParticleVector2D = std::vector<GenParticlePtrVector>;


class GenParticle {
public:
    size_t index;
    int pdg;
    int status;
    GenStatusFlags genStatusFlags;
    LorentzVectorM_Float momentum;
    GenParticlePtrVector mothers;
    GenParticlePtrVector daughters;

public:
    GenParticle(const ntuple::Event& events, size_t n)
        : index(n)
    {
        if(n >= events.genParticles_p4.size())
            throw std::runtime_error("events index is out of range.");

        pdg = events.genParticles_pdg.at(n);
        status = events.genParticles_status.at(n);
        genStatusFlags = GenStatusFlags(events.genParticles_statusFlags.at(n));
        momentum = events.genParticles_p4.at(n);
    }
};

class GenEvent {
public:
    GenParticleVector genParticles;
    ParticleCodeMap particleCodeMap;
    ParticleCodeMap hardParticleCodeMap;
    GenParticleSet primaryParticles;

public:
    GenEvent(const ntuple::Event& event)
    {
        genParticles.reserve(event.genParticles_p4.size());
        for(size_t n = 0; n < event.genParticles_p4.size(); ++n)
            genParticles.emplace_back(event,n);

        for (size_t relIndex = 0; relIndex < event.genParticles_rel_mIndex.size(); relIndex++) {
            size_t motherIndex = event.genParticles_rel_mIndex.at(relIndex);
            size_t daughterIndex = event.genParticles_rel_pIndex.at(relIndex);

            GenParticle* mother = &genParticles.at(motherIndex);
            GenParticle* daughter = &genParticles.at(daughterIndex);

            mother->daughters.push_back(daughter);
            daughter->mothers.push_back(mother);
        }

        for(const GenParticle& genParticle : genParticles ){
            if (!genParticle.mothers.size())
                primaryParticles.insert(&genParticle);


            // if (genParticle.genStatusFlags.isPrompt()) // particle.genStatusFlags.isLastCopy()){
                particleCodeMap[std::abs(genParticle.pdg)].insert(&genParticle);

        // else if(genParticle.genStatusFlags.fromHardProcess())
        //     hardParticleCodeMap[genParticle.pdg].insert(&genParticle);
        }
    }

    GenParticleSet GetParticles(int particle_pgd/*, double pt = 0*/) const {
        GenParticleSet results;
        const ParticleCodeMap::const_iterator code_iter = particleCodeMap.find(particle_pgd);

        if (code_iter != particleCodeMap.end()){
            for (const GenParticle* particle : code_iter->second){
                // if (particle->momentum.Pt() <= pt) continue;
                if(!particle->genStatusFlags.isLastCopy()) continue;
                results.insert(particle);
            }
        }
        return results;
    }

    bool areParented(std::set<const GenParticle*>& daughters, const GenParticle& possible_mother){
        bool family = false;
        for(const auto& daughter: daughters){
            if(!daughter->daughters.size()) {
                for (size_t i = 0; i < daughter->mothers.size(); i++) {
                    if((daughter->mothers.at(i)) == &possible_mother){
                        std::cout << "Nutella" << '\n';
                        family = true;
                    }
                }
            }
        }
        for(const auto& nonna : possible_mother.mothers)
            areParented(daughters,*nonna);
        return family;
    }

    GenParticleSet GetBaryons(int particle_pgd) const {
        GenParticleSet results;
        const ParticleCodeMap::const_iterator code_iter = particleCodeMap.find(particle_pgd);

        if (code_iter != particleCodeMap.end()){
            for (const GenParticle* particle : code_iter->second){
                if(!particle->genStatusFlags.isLastCopy() && !particle->genStatusFlags.isPrompt()
                    && particle_types[particle_pgd] != "baryon" ) continue;
                results.insert(particle);
            }
        }
        return results;
    }

    GenParticleSet GetMeson(int particle_pgd) const {
        GenParticleSet results;
        const ParticleCodeMap::const_iterator code_iter = particleCodeMap.find(particle_pgd);

        if (code_iter != particleCodeMap.end()){
            for (const GenParticle* particle : code_iter->second){
                if(!particle->genStatusFlags.isLastCopy() && !particle->genStatusFlags.isPrompt()
                    && particle_types[particle_pgd] != "meson" ) continue;
                results.insert(particle);
            }
        }
        return results;
    }

    static const std::string& GetParticleName(int pdgId)
    {
        auto iter = particle_names.find(pdgId);
        if(iter == particle_names.end()) throw exception("Name not found for particle with pdgId = %1%") % pdgId;
        return iter->second;
    }

    static void intializeNames(const std::string& fileName, const std::string& typesFileName)
    {
        particle_names.clear();
        particle_types.clear();
        std::ifstream f(fileName.c_str());
        std::ifstream g(typesFileName.c_str());
        if(!f.is_open() || !g.is_open())
            throw analysis::exception("Unable to read the configuration file ");
        while(f.good()) {
            std::string line;
            std::getline(f, line);
            if(!line.length() || line[0] == '#')
                continue;
            auto value_name = SplitValueList(line, false, ",");
            if(value_name.size() != 2)
                throw exception("Invalid particle definition: '%1%'") % line;
            const int value = Parse<int>(value_name.at(0)); // cambiare nome a PdgID
            const std::string& name = value_name.at(1);
            if(particle_names.count(value))
                throw exception("Duplicated definition of particle with pdgId = %1%") % value;
            particle_names[value] = name;
        }
        while(g.good()) {
            std::string line;
            std::getline(f, line);
            if(!line.length() || line[0] == '#')
                continue;
            auto value_name = SplitValueList(line, false, ",");
            if(value_name.size() != 2)
                throw exception("Invalid particle definition: '%1%'") % line;
            const int value = Parse<int>(value_name.at(0)); // cambiare nome a PdgID
            const std::string& type = value_name.at(1);
            if(particle_types.count(value))
                throw exception("Duplicated definition of particle with pdgId = %1%") % value;
            particle_types[value] = type;
        }
    }

    void Print() const
    {
        for (const GenParticle* particle : primaryParticles) {
            PrintChain(particle);
        }
    }

    void PrintChain(const GenParticle* particle, unsigned iteration = 0) const
    {
        const int pdgParticle = particle->pdg;
        const auto particleName = GetParticleName(pdgParticle);
        const int particleStatus = particle->status;
        const LorentzVectorM_Float genParticle_momentum = particle->momentum;
        const std::bitset<15> genStatusFlags_ = particle->genStatusFlags.flags ;
        for (unsigned n = 0; n < iteration; ++n)
            std::cout << "  ";
        auto mother_index = particle->mothers.size() > 0 ?  particle->mothers.at(0)->index : -1;
        std::cout <<particleName  << " <" << pdgParticle << ">" <<  " pt= " << genParticle_momentum.Pt()
            << " eta= " << genParticle_momentum.Eta() << " phi= " << genParticle_momentum.Phi()
            << " E= " << genParticle_momentum.E()  << " m= " << genParticle_momentum.M() << " index=" << particle->index
            << " mother_index=" << mother_index << " status=" << particleStatus  << " statusFlags=" << genStatusFlags_ << std::endl;

        for(unsigned n = 0; n < particle->daughters.size(); ++n) {
            const GenParticle* daughter = particle->daughters.at(n);
                PrintChain(daughter ,iteration+1);
        }
    }



private:
    static std::map<int, std::string> particle_names;
    static std::map<int, std::string> particle_types;
};

std::map<int, std::string> GenEvent::particle_names;
std::map<int, std::string> GenEvent::particle_types;

class VisibleGenObject {
public:
    const GenParticle* origin;

    GenParticleSet finalStateChargedLeptons;
    GenParticleSet finalStateChargedHadrons;
    GenParticleSet finalStateNeutralHadrons;
    GenParticleSet finalStateBquarks;


    LorentzVectorM_Float chargedLeptonsMomentum;
    LorentzVectorM_Float chargedHadronsMomentum;
    LorentzVectorM_Float finalStateBquarkMomentum;

    LorentzVectorM_Float neutralHadronsMomentum;
    LorentzVectorM_Float visibleMomentum;
    LorentzVectorM_Float invisibleMomentum;

    GenParticleSet particlesProcessed;

public:

    explicit VisibleGenObject() : origin(nullptr) {}

    explicit VisibleGenObject(const GenParticle *_origin) : origin(_origin)
    {
        const GenParticleSet particlesToIgnore;
        CollectInfo(origin, particlesProcessed, particlesToIgnore);
    }

    VisibleGenObject(const GenParticle *_origin, const GenParticleSet& particlesToIgnore) : origin(_origin)
    {
        CollectInfo(origin, particlesProcessed, particlesToIgnore);
    }

    bool operator < (const VisibleGenObject& other) const
    {
        return origin < other.origin;
    }


private:
    void CollectInfo(const GenParticle* particle, GenParticleSet& particlesProcessed,
                     const GenParticleSet& particlesToIgnore)
    {
        // if(particle->status == 1 && particle->daughters.size() != 0)
        // if(particle->daughters.size() != 0) //bool isLastCopy()
        //     throw exception("Invalid gen particle");

        if(particlesProcessed.count(particle) || particlesToIgnore.count(particle)) return;
        particlesProcessed.insert(particle);
        if(particle->genStatusFlags.isLastCopy()) { //prima era status == FinalStateParticle
            if(particles::neutrinos.count(particle->pdg)) {
                invisibleMomentum += particle->momentum;
                return;
            }

            visibleMomentum += particle->momentum;
            if(std::abs(particle->pdg) == particles::ParticleCode::e || std::abs(particle->pdg) == particles::ParticleCode::mu) {
                finalStateChargedLeptons.insert(particle);
                chargedLeptonsMomentum += particle->momentum;
            } else /*if(particle->pdg > 0) */ {
                finalStateChargedHadrons.insert(particle);
                chargedHadronsMomentum += particle->momentum;
            } //else {
            //     finalStateNeutralHadrons.insert(particle);
            //     neutralHadronsMomentum += particle->momentum;
            // }

            // else if(particle->pdg == particles::ParticleCode::b) {
            //     finalStateBquarks.insert(particle);
            //     finalStateBquarkMomentum += particle->momentum;
            // }
        }
        //
        for(const GenParticle* daughter : particle->daughters){
            CollectInfo(daughter, particlesProcessed, particlesToIgnore);
        }
    }
};

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

typedef std::vector<VisibleGenObject> VisibleGenObjectVector;

inline bool HasMatchWithMCObject(const GenParticle& particle, std::vector<const GenParticle*>& visible_daughters,
                                 const LorentzVectorM_Float& momentum_reco, double deltaR, bool useVisibleMomentum = false)
{
    auto visibleMomentum = GetFinalStateMomentum(particle, visible_daughters, true, false);
    const LorentzVectorM_Float& momentum = useVisibleMomentum ? visibleMomentum : particle.momentum;
    return ROOT::Math::VectorUtil::DeltaR(momentum, momentum_reco) < deltaR;
}
} //analysis
