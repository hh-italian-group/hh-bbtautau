#pragma once

#include <set>
#include <TLorentzVector.h>

#include "h-tautau/Core/include/EventTuple.h"
#include "GenStatusFlags.h"
#include "hh-bbtautau/Analysis/include/Particle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

// #include "Particle.h"
#include "exception.h"

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


            if (genParticle.genStatusFlags.isPrompt()) // particle.genStatusFlags.isLastCopy()){
                particleCodeMap[std::abs(genParticle.pdg)].insert(&genParticle);

        // else if(genParticle.genStatusFlags.fromHardProcess())
        //     hardParticleCodeMap[genParticle.pdg].insert(&genParticle);
        }
    }

    GenParticleSet GetParticles(int particle_pgd, double pt = 0) const {
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

    void Print() const
    {
        for (const GenParticle* particle : primaryParticles) {
            PrintChain(particle);
        }
    }
    
    void PrintChain(const GenParticle* particle, unsigned iteration = 0) const
    {
        const int pdgParticle = particle->pdg;
        const int particleStatus = particle->status;
        const LorentzVectorM_Float genParticle_momentum = particle->momentum;
        // const GenStatusFlags genStatusFlags = particle->genStatusFlags;
        for (unsigned n = 0; n < iteration; ++n)
            std::cout << "  ";
        std::cout << "index=" << particle->index << " name=" << pdgParticle << " status=" << particleStatus
                  <<  " pt= " << genParticle_momentum.Pt() <<"\n";
        for(unsigned n = 0; n < particle->daughters.size(); ++n) {
            const GenParticle* daughter = particle->daughters.at(n);
                PrintChain(daughter,iteration+1);
        }
    }
};
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

// inline bool HasMatchWithMCObject(const LorentzVectorM_Float& candidateMomentum, const LorentzVectorM_Float& momentum_reco,
//                                  const VisibleGenObject* genObject, double deltaR, bool useVisibleMomentum = false)
inline bool HasMatchWithMCObject(const GenParticle& particle, std::vector<const GenParticle*>& visible_daughters,
                                 const LorentzVectorM_Float& momentum_reco, double deltaR, bool useVisibleMomentum = false)
{
    auto visibleMomentum = GetFinalStateMomentum(particle, visible_daughters, true, false);
    const LorentzVectorM_Float& momentum = useVisibleMomentum ? visibleMomentum : particle.momentum;
    return ROOT::Math::VectorUtil::DeltaR(momentum, momentum_reco) < deltaR;
}

// inline bool HasMatchWithMCObject(const LorentzVectorM_Float& candidateMomentum_1,
//                                  const LorentzVectorM_Float& candidateMomentum_2,
//                                  double deltaR)
// {
//     return ROOT::Math::VectorUtil::DeltaR(candidateMomentum_1, candidateMomentum_2) < deltaR;
// }

// template<typename Container>
// inline VisibleGenObjectVector FindMatchedObjects(const LorentzVectorM_Float& candidateMomentum_1,
//                                                  const LorentzVectorM_Float& candidateMomentum_2,
//                                                  const Container& genObjects, double deltaR)
// {
//     VisibleGenObjectVector matchedGenObjects;
//     for (const VisibleGenObject& genObject : genObjects){
//         if (HasMatchWithMCObject(candidateMomentum_1, candidateMomentum_2, &genObjects, deltaR))
//             matchedGenObjects.push_back(genObject);
//     }
//     return matchedGenObjects;
// }



} //analysis
