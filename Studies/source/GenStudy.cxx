/*! Definition of data and event categories used in HH->bbTauTau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/GenParticle.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Analysis/src/EventInfo.cpp"


struct Arguments {
    REQ_ARG(std::string, inputPath);
    REQ_ARG(std::string, treeName);
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
     TH1D_ENTRY(bm_1, 10, -0.5, 9.5)
     TH1D_ENTRY(bm_2, 10, -0.5, 9.5)
     TH1D_ENTRY(bm, 10, -0.5, 9.5)
     TH1D_ENTRY(bm_VM, 10, -0.5, 9.5)
     TH1D_ENTRY(reco_0jet_matched_gen_particles, 10, -0.5, 9.5)
     TH1D_ENTRY(reco_1jet_matched_gen_particles, 10, -0.5, 9.5)
 };

class GenStudy {
public:
    GenStudy(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.outputFile())), anaData(output)
        {
            GenEvent::intializeNames(args.particleNameTypeFile());
        }

    size_t getIndexMaxPt(std::vector<const GenParticle*> visible_daughters){
        size_t max_pt_index = 0;
        for (size_t i = 1; i < visible_daughters.size(); i++)
            max_pt_index = visible_daughters.at(0)->momentum.pt() > visible_daughters.at(i)->momentum.pt() ? 0 : (i);
        return max_pt_index;
    }

    void Run()
    {
        auto file = root_ext::OpenRootFile(args.inputPath());
        auto tuple = ntuple::CreateEventTuple(args.treeName(), file.get(), true, ntuple::TreeState::Full);
        const double deltaR_value = 0.2;

        for(const auto& event : *tuple) {
           if(event.eventEnergyScale != 0) continue;

            GenEvent genEvent(event);
//            genEvent.Print();
            const GenParticleSet higgsPair = genEvent.GetParticles(particles::ParticleCode::higgs);
            if(higgsPair.size() != 2)
                throw analysis::exception("Higgs pair size must be 2, insteed has a size of '%1%'") %higgsPair.size(); //Stampare eventId e albero

            const GenParticle* h_tautau = nullptr;
            const GenParticle* h_bb = nullptr;

            for(const auto& higgs : higgsPair){
                if(higgs->daughters.size() != 2)
                throw analysis::exception("Higgs decay products must be 2, insteed has a size of '%1%'") %higgs->daughters.size();

                if(std::abs(higgs->daughters.at(0)->pdg) == particles::ParticleCode::tau &&
                    std::abs(higgs->daughters.at(1)->pdg) == particles::ParticleCode::tau)
                    h_tautau = higgs;
                else if(std::abs(higgs->daughters.at(0)->pdg) == particles::ParticleCode::b &&
                    std::abs(higgs->daughters.at(1)->pdg) == particles::ParticleCode::b)
                     h_bb = higgs;
                 else
                     throw analysis::exception("Higgs pair doesn't decay in a pair of taus or b jets h1 -> '%1%', h2 -> '%2%'.") %higgs->daughters.at(0)->pdg %higgs->daughters.at(1)->pdg;
            }
            if(h_tautau == nullptr || h_bb == nullptr)
                throw analysis::exception("The Higgs pairs are not filled");
            if(h_tautau->daughters.size() != 2 || h_bb->daughters.size() != 2)
                throw analysis::exception("Each of the Higgs pairs must decay in two particles.");

            int n_tau_matches = 0;
            for (size_t tau_index = 0; tau_index < h_tautau->daughters.size(); tau_index++) {
                std::vector<const GenParticle*> visible_daughters;
                if(std::abs(h_tautau->daughters.at(tau_index)->momentum.eta()) > 2.3 ) continue;
                auto visibleMomentum = genEvent.GetFinalStateMomentum(*h_tautau->daughters.at(tau_index), visible_daughters, true, false);
                for (size_t reco_tau_index = 0; reco_tau_index < 2; reco_tau_index++) { 
                    auto reco_tau_momentum = reco_tau_index == 0 ? event.p4_1 : event.p4_2;
                    if(ROOT::Math::VectorUtil::DeltaR(visibleMomentum, reco_tau_momentum) < deltaR_value)
                        ++n_tau_matches;
                }
            }

            const auto baryons_plus_mesons = genEvent.GetTypesParticles({"baryon", "meson"}, h_bb);
            if(baryons_plus_mesons.size() < 2)
                throw analysis::exception("There're less than two barions or mesons");

            size_t max_i = 0;
            std::vector<const GenParticle*> most_energetic_pair;
            for (size_t bm_index = 0; bm_index < baryons_plus_mesons.size(); bm_index++) {
                if(baryons_plus_mesons.at(bm_index)->momentum.E() > baryons_plus_mesons.at(max_i)->momentum.E())
                    max_i = bm_index;
            }

            size_t second_max_i = max_i == 0 ? 1 : 0;
            for (size_t bm_index = 0; bm_index < baryons_plus_mesons.size(); bm_index++) {
                if(bm_index == max_i) continue;
                if(baryons_plus_mesons.at(bm_index)->momentum.E() > baryons_plus_mesons.at(second_max_i)->momentum.E())
                    second_max_i = bm_index;
            }

            most_energetic_pair.push_back(baryons_plus_mesons.at(max_i));
            most_energetic_pair.push_back(baryons_plus_mesons.at(second_max_i));

            int mb_matches_1 = 0;
            int mb_matches_2 = 0;
            int mb_matches = 0;
            int mb_matches_VM = 0;
            int test_reco_0 = 0;
            int test_reco_1 = 0;

            for (size_t i = 0; i < most_energetic_pair.size(); i++) {
                if(std::abs(most_energetic_pair.at(i)->momentum.eta()) > 2.3) continue;
                std::vector<const GenParticle*> visible_daughters;
                auto visibleMomentum = genEvent.GetFinalStateMomentum(*most_energetic_pair.at(i), visible_daughters, false, false);

                for (size_t reco_b_index = 0; reco_b_index < event.jets_p4.size(); reco_b_index++) {
                    if(ROOT::Math::VectorUtil::DeltaR(visibleMomentum, event.jets_p4.at(reco_b_index)) < deltaR_value)
                        ++mb_matches_VM;

                    if(ROOT::Math::VectorUtil::DeltaR(most_energetic_pair.at(i)->momentum, event.jets_p4.at(reco_b_index)) < deltaR_value)
                        ++mb_matches;

                    if(ROOT::Math::VectorUtil::DeltaR(most_energetic_pair.at(0)->momentum, event.jets_p4.at(reco_b_index)) < deltaR_value)
                        ++mb_matches_1;

                    if(ROOT::Math::VectorUtil::DeltaR(most_energetic_pair.at(1)->momentum, event.jets_p4.at(reco_b_index)) < deltaR_value)
                        ++mb_matches_2;

                }
                if(ROOT::Math::VectorUtil::DeltaR(most_energetic_pair.at(i)->momentum, event.jets_p4.at(0)) < deltaR_value)
                    ++test_reco_0;
                if(event.jets_p4.size() >1){
                    if(ROOT::Math::VectorUtil::DeltaR(most_energetic_pair.at(i)->momentum, event.jets_p4.at(1)) < deltaR_value)
                        ++test_reco_1;
                }
            }

            anaData.h_tautau_matches().Fill(n_tau_matches);
            anaData.bm_1().Fill(mb_matches_1);
            anaData.bm_2().Fill(mb_matches_2);
            anaData.bm().Fill(mb_matches);
            anaData.bm_VM().Fill(mb_matches_VM);
            anaData.reco_0jet_matched_gen_particles().Fill(test_reco_0);
            anaData.reco_1jet_matched_gen_particles().Fill(test_reco_1);
        }
    }

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    FileMerger anaData;
};

} // namespace analysis
PROGRAM_MAIN(analysis::GenStudy, Arguments)
