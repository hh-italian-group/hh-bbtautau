/*! Definition of data and event categories used in HH->bbTauTau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/GenParticle.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Analysis/src/EventInfo.cpp"
#include <functional>   // std::greater
#include <iostream>     // std::cout
#include <algorithm>

struct Arguments {
    REQ_ARG(std::string, inputPath);
    REQ_ARG(analysis::Channel, channel);
    REQ_ARG(std::string, outputFile);
    REQ_ARG(std::string, particleNameTypeFile);

};

namespace analysis {

using Event = ntuple::Event;
using EventPtr = std::shared_ptr<Event>;

class FileMerger : public root_ext::AnalyzerData {
 public:
     using AnalyzerData::AnalyzerData;

    TH1D_ENTRY(h_tautau_matches, 10, -0.5, 9.5)
    TH1D_ENTRY(baryons_mesons_leg_1_matches, 10, -0.5, 9.5)
    TH1D_ENTRY(baryons_mesons_leg_2_matches, 10, -0.5, 9.5)
    TH1D_ENTRY(baryons_mesons_matches, 10, -0.5, 9.5)
    TH1D_ENTRY(baryons_mesons_matches_visible_mass, 10, -0.5, 9.5)
    TH1D_ENTRY(reco_0jet_matched_gen_particles, 10, -0.5, 9.5)
    TH1D_ENTRY(reco_1jet_matched_gen_particles, 10, -0.5, 9.5)
    TH2D_ENTRY_EX(n_baryons_vs_n_mesons, 10, -0.5, 9.5, 10, -0.5, 9.5, "number of baryons", "number of mesons", false, 1, true)

    TH1D_ENTRY_EX(rel_energy, 80, 0, 1.5, "Rel Energy", "events", true, 1, false, true)
    TH1D_ENTRY_EX(rel_momentum, 80, 0, 1.5, "Rel Pt", "events", true, 1, false, true)
    TH1D_ENTRY_EX(pt_tau_1, 50, 0, 500, "Pt [GeV]", "events", true, 1, false, true)

    TH1D_ENTRY_EX(pt_charged_vis_tau, 50, 0, 500, "Pt [GeV]", "events", true, 1, false, true)

    TH1D_ENTRY_EX(pt_tau_2, 50, 0, 500, "Pt [GeV]", "events", true, 1, false, true)
    TH1D_ENTRY_EX(eta_tau_1, 20, -2.5, 2.5,"#eta ", "events", true, 1, false, true)
    TH1D_ENTRY_EX(eta_tau_2, 20, -2.5, 2.5, "#eta ", "events", true, 1, false, true)

//    TH2D_ENTRY_EX(sum_shared_sum_normal, 10, -0.5, 9.5, 10, -0.5, 9.5, "sum of baryons plus mesons", "sum of baryons plus mesons shared", false, 1, true)
//    TH1D_ENTRY_EX(sum_bm, 10, -0.5, 9.5, "Sum of baryons and mesons", "events", true, 1, false, true)
//    TH1D_ENTRY_EX(sum_bm_shared, 10, -0.5, 9.5, "Sum of baryons and mesons", "events", true, 1, false, true)
 };

enum class GenDecayMode { Electron, Muon, Hadrons };

struct HHGenEvent {
    const GenParticle *h_tautau, *h_bb;
    std::array<const GenParticle*, 2> tau, b;
    std::array<LorentzVectorM, 2> vis_tau, vis_b;
    std::array<LorentzVectorXYZ, 2> vis_charged_tau;
    std::array<GenDecayMode, 2> tau_decay;
};

class GenStudy {
public:
    GenStudy(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.outputFile())), anaData(output)
        {
            GenEvent::intializeNames(args.particleNameTypeFile());
        }

    bool isInsideAcceptance(const HHGenEvent& hh_event, Channel channel)
    {
        for(size_t tau_index = 0; tau_index < hh_event.h_tautau->daughters.size(); ++tau_index){
            if(std::abs(hh_event.vis_tau.at(tau_index).eta()) > 2.3 ||
                hh_event.vis_tau.at(tau_index).Pt() < 20 )
                return false;
        }

        if(channel == Channel::ETau){
            if((hh_event.tau_decay.at(0) == GenDecayMode::Electron &&  hh_event.tau_decay.at(1)  == GenDecayMode::Hadrons) ||
                    (hh_event.tau_decay.at(1) == GenDecayMode::Electron &&  hh_event.tau_decay.at(0)  == GenDecayMode::Hadrons))
                return true;
        }
        else if(channel == Channel::MuTau){
            if((hh_event.tau_decay.at(0) == GenDecayMode::Muon &&  hh_event.tau_decay.at(1)  == GenDecayMode::Hadrons) ||
                    (hh_event.tau_decay.at(1) == GenDecayMode::Muon &&  hh_event.tau_decay.at(0)  == GenDecayMode::Hadrons))
                return true;
        }
        else if(channel == Channel::TauTau){
            if((hh_event.tau_decay.at(0) == GenDecayMode::Hadrons &&  hh_event.tau_decay.at(1)  == GenDecayMode::Hadrons))
                return true;
        }
//        else
            return false;
    }

    struct greater
    {
        template<class T>
        bool operator()(T const &a, T const &b) const { return a > b; }
    };

    void Run()
    {
        auto file = root_ext::OpenRootFile(args.inputPath());
        auto tuple = ntuple::CreateEventTuple(ToString(args.channel()), file.get(), true, ntuple::TreeState::Full);
        const double deltaR_value = 0.2;

        for(const auto& event : *tuple) {
           if(event.eventEnergyScale != 0) continue;

//           const EventIdentifier EventId(event.run, event.lumi, event.evt);
//           const EventIdentifier EventIdTest(1,386,385640);
//           if(!(EventId == EventIdTest)) continue;

            GenEvent genEvent(event);
//            genEvent.Print();
            const GenParticleSet higgsPair = genEvent.GetParticles(particles::ParticleCode::higgs, true);
            if(higgsPair.size() != 2)
                throw analysis::exception("Higgs pair size must be 2, insteed has a size of '%1%'") %higgsPair.size(); //Stampare eventId e albero

            HHGenEvent HH_Gen_Event;

            HH_Gen_Event.h_tautau = nullptr;
            HH_Gen_Event.h_bb = nullptr;

            for(const auto& higgs : higgsPair){
                if(higgs->daughters.size() != 2)
                    throw analysis::exception("Higgs decay products must be 2, insteed has a size of '%1%'") %higgs->daughters.size();

                if(std::abs(higgs->daughters.at(0)->pdg) == particles::ParticleCode::tau &&
                    std::abs(higgs->daughters.at(1)->pdg) == particles::ParticleCode::tau)
                    HH_Gen_Event.h_tautau = higgs;
                else if(std::abs(higgs->daughters.at(0)->pdg) == particles::ParticleCode::b &&
                    std::abs(higgs->daughters.at(1)->pdg) == particles::ParticleCode::b)
                    HH_Gen_Event.h_bb = higgs;
                 else
                     throw analysis::exception("Higgs pair doesn't decay in a pair of taus or b jets h1 -> '%1%', h2 -> '%2%'.") %higgs->daughters.at(0)->pdg %higgs->daughters.at(1)->pdg;
            }

            if(HH_Gen_Event.h_tautau == nullptr || HH_Gen_Event.h_bb == nullptr)
                throw analysis::exception("The Higgs pairs are not filled");
            if(HH_Gen_Event.h_tautau->daughters.size() != 2 || HH_Gen_Event.h_bb->daughters.size() != 2)
                throw analysis::exception("Each of the Higgs pairs must decay in two particles.");

            for (size_t tau_index = 0; tau_index < HH_Gen_Event.h_tautau->daughters.size(); tau_index++) {
                std::vector<const GenParticle*> visible_daughters;

                HH_Gen_Event.tau[tau_index] = HH_Gen_Event.h_tautau->daughters.at(tau_index);
                HH_Gen_Event.vis_tau[tau_index] = genEvent.GetFinalStateMomentum(*HH_Gen_Event.h_tautau->daughters.at(tau_index),
                                                                                 visible_daughters, true, false);

               auto particle_charge = GenEvent::particleCharge();
               for(size_t vis_daughter_index = 0 ; vis_daughter_index <visible_daughters.size(); ++vis_daughter_index ){
                    if (std::abs(visible_daughters.at(vis_daughter_index)->pdg) == particles::ParticleCode::e )
                        HH_Gen_Event.tau_decay[tau_index] = GenDecayMode::Electron;
                    else if (std::abs(visible_daughters.at(vis_daughter_index)->pdg) == particles::ParticleCode::mu )
                        HH_Gen_Event.tau_decay[tau_index] = GenDecayMode::Muon;
                    else
                        HH_Gen_Event.tau_decay[tau_index] = GenDecayMode::Hadrons;

                    if(particle_charge[visible_daughters.at(vis_daughter_index)->pdg] != 0)
                        HH_Gen_Event.vis_charged_tau[tau_index] += visible_daughters.at(vis_daughter_index)->momentum;

                }
            }
            if(!isInsideAcceptance(HH_Gen_Event, args.channel())) continue;

            size_t n_tau_matches_total = 0;
            std::map<size_t, std::set<size_t>> tau_matches;
            for (size_t tau_index = 0; tau_index < HH_Gen_Event.h_tautau->daughters.size(); tau_index++) {
                for (size_t reco_tau_index = 0; reco_tau_index < event.lep_p4.size(); reco_tau_index++) {
                    const auto& visibleMomentum = HH_Gen_Event.vis_tau[tau_index];
                    if(ROOT::Math::VectorUtil::DeltaR(visibleMomentum, event.lep_p4.at(reco_tau_index)) < deltaR_value) {
                        ++n_tau_matches_total;
                        tau_matches[tau_index].insert(reco_tau_index);
                    }
//                    std::cout <<"tau_index=" << tau_index<< " reco_tau_index=" << reco_tau_index << " n_tau_matches="
//                             << n_tau_matches_total << " Pt=" << visibleMomentum.Pt() << " Eta=" << visibleMomentum.Eta()
//                             << " Phi=" << visibleMomentum.Phi() << " Pt_reco=" << event.lep_p4.at(reco_tau_index).Pt() << " Eta_reco=" << event.lep_p4.at(reco_tau_index).Eta()
//                             << " Phi_reco=" << event.lep_p4.at(reco_tau_index).Phi()  << " gen_match="<<
//                             event.lep_gen_match.at(reco_tau_index)<< std::endl;

                }
            }

            GenParticleSet baryons_plus_mesons_set;
            std::vector<const GenParticle*> baryons_plus_mesons;

            genEvent.GetTypesParticles({particles::ParticleCode::ParticleType::baryon, particles::ParticleCode::ParticleType::meson},
                                         HH_Gen_Event.h_bb, baryons_plus_mesons_set);

            for (auto bm : baryons_plus_mesons_set)
                baryons_plus_mesons.push_back((bm));

            if(baryons_plus_mesons.size() < 2)
                throw analysis::exception("There're less than two barions or mesons");

//             if(baryons_plus_mesons.size() == 3)
//                 std::cout << event.run <<"," << event.lumi << "," << event.evt << std::endl;

            std::vector<std::pair<int,int>> energy_index;
            for (size_t bm_index = 0; bm_index < baryons_plus_mesons.size(); bm_index++)
                energy_index.push_back(std::make_pair(baryons_plus_mesons.at(bm_index)->momentum.E(), bm_index));
            std::sort(energy_index.begin(), energy_index.end(), greater());

            std::vector<const GenParticle*> most_energetic_pair;
            most_energetic_pair.push_back(baryons_plus_mesons.at(static_cast<size_t>(energy_index.at(0).second)));
            most_energetic_pair.push_back(baryons_plus_mesons.at(static_cast<size_t>(energy_index.at(1).second)));

            double energy_of_the_pair =  most_energetic_pair.at(0)->momentum.E() + most_energetic_pair.at(1)->momentum.E();

            double total_energy = 0;

            double momentum_of_the_pair =  most_energetic_pair.at(0)->momentum.Pt() + most_energetic_pair.at(1)->momentum.Pt();
            double total_momentum = 0;

            for (size_t i = 0; i < baryons_plus_mesons.size(); i++){
                total_energy += baryons_plus_mesons.at(i)->momentum.E();
                total_momentum += baryons_plus_mesons.at(i)->momentum.Pt();
            }

            double relative_energy = 0;
            double relative_momentum = 0;
            if(baryons_plus_mesons.size() > 2)  {
                 relative_energy = energy_of_the_pair /total_energy;
                 relative_momentum = momentum_of_the_pair /total_momentum;
            }

            std::vector<int> variables (10);
            for (size_t i = 0; i < most_energetic_pair.size(); i++) {
//                if(!isInsideAcceptance(HH_Gen_Event, args.channel())) continue;

                std::vector<const GenParticle*> visible_daughters;
                auto visibleMomentum = genEvent.GetFinalStateMomentum(*most_energetic_pair.at(i), visible_daughters, false, false);

                for (size_t reco_b_index = 0; reco_b_index < event.jets_p4.size(); reco_b_index++) {
                    if(ROOT::Math::VectorUtil::DeltaR(visibleMomentum, event.jets_p4.at(reco_b_index)) < deltaR_value)
                        ++variables.at(0);// ++mb_matches_VM;

                    if(ROOT::Math::VectorUtil::DeltaR(most_energetic_pair.at(i)->momentum, event.jets_p4.at(reco_b_index)) < deltaR_value)
                        ++variables.at(1);//++mb_matches;

                    if(ROOT::Math::VectorUtil::DeltaR(most_energetic_pair.at(0)->momentum, event.jets_p4.at(reco_b_index)) < deltaR_value)
                        ++variables.at(2);// ++mb_matches_1;

                    if(ROOT::Math::VectorUtil::DeltaR(most_energetic_pair.at(1)->momentum, event.jets_p4.at(reco_b_index)) < deltaR_value)
                        ++variables.at(3);//++mb_matches_2;

                }
                if(ROOT::Math::VectorUtil::DeltaR(most_energetic_pair.at(i)->momentum, event.jets_p4.at(0)) < deltaR_value)
                    ++variables.at(4);//++test_reco_0;
                if(event.jets_p4.size() >1){
                    if(ROOT::Math::VectorUtil::DeltaR(most_energetic_pair.at(i)->momentum, event.jets_p4.at(1)) < deltaR_value)
                        ++variables.at(5);//++LorentzVectorM;
                }

            }
            anaData.h_tautau_matches().Fill(n_tau_matches_total);
//            if(n_tau_matches == 1 )
//                std::cout << event.run << "," << event.lumi << "," << event.evt << std::endl;

            anaData.baryons_mesons_leg_1_matches().Fill(variables.at(2));
            anaData.baryons_mesons_leg_2_matches().Fill(variables.at(3));
            anaData.baryons_mesons_matches().Fill(variables.at(1));
            anaData.baryons_mesons_matches_visible_mass().Fill(variables.at(0));
            anaData.reco_0jet_matched_gen_particles().Fill(variables.at(4));
            anaData.reco_1jet_matched_gen_particles().Fill(variables.at(5));

            anaData.rel_energy().Fill(relative_energy);
            anaData.rel_momentum().Fill(relative_momentum);

            anaData.pt_tau_1().Fill(HH_Gen_Event.h_tautau->daughters.at(0)->momentum.Pt());
            anaData.pt_tau_2().Fill(HH_Gen_Event.h_tautau->daughters.at(1)->momentum.Pt());
            anaData.eta_tau_1().Fill(HH_Gen_Event.h_tautau->daughters.at(0)->momentum.Eta());
            anaData.eta_tau_2().Fill(HH_Gen_Event.h_tautau->daughters.at(1)->momentum.Eta());

            if(n_tau_matches_total == 1) {
                for(size_t tau_index = 0; tau_index < 2; ++tau_index) {
                    if(tau_matches[tau_index].empty()){
                        anaData.pt_charged_vis_tau("missing").Fill(HH_Gen_Event.vis_charged_tau[tau_index].Pt());
//                        if(HH_Gen_Event.vis_charged_tau[tau_index].Pt() > 50 )
//                            std::cout << event.run << "," << event.lumi << "," << event.evt << std::endl;
                    }
                }

            }

//            anaData.sum_bm().Fill(baryons_plus_mesons.size());
//            anaData.sum_bm_shared().Fill(baryons_plus_mesons_shared.size());
//            if(baryons.size() + mesons.size() == 4)
//              std::cout << event.run << "," << event.lumi << "," << event.evt << std::endl;


        }
    }

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    FileMerger anaData;
};

} // namespace analysis
PROGRAM_MAIN(analysis::GenStudy, Arguments)
