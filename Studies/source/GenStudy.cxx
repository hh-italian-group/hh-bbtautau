#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "hh-bbtautau/Analysis/include/GenParticle.h"
#include "hh-bbtautau/Analysis/include/Particle.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"

#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "h-tautau/Core/include/AnalysisTypes.h"


struct Arguments {
    REQ_ARG(std::string, inputPath);
    REQ_ARG(std::string, treeName);
    REQ_ARG(std::string, outputFile);
    REQ_ARG(std::string, particleNameFile);
    REQ_ARG(std::string, particleTypeFile);

};

namespace analysis {

using Event = ntuple::Event;
using EventPtr = std::shared_ptr<Event>;

class FileMerger : public root_ext::AnalyzerData {
 public:
     using AnalyzerData::AnalyzerData;

     TH1D_ENTRY(h_tautau_matches, 10, -0.5, 9.5)
     TH1D_ENTRY(h_bb_matchesV1, 10, -0.5, 9.5)
     TH1D_ENTRY(h_bb_matchesV2, 10, -0.5, 9.5)
     TH1D_ENTRY(B_star_zero, 10, -0.5, 9.5)
     TH1D_ENTRY(num, 10, -0.5, 9.5)
     TH1D_ENTRY(number_of_b, 10, -0.5, 9.5)
     TH1D_ENTRY(dR_1, 50, 0, 1)
     TH1D_ENTRY(dR_orignial_b_vs_VM, 50, 0, 1)
     TH1D_ENTRY(dR_2, 50, 0, 1)
     TH1D_ENTRY(dR_max_pt_vs_original_momentum, 50, 0, 1)
     TH1D_ENTRY(pt, 50, 0, 50)
     TH1D_ENTRY(pt_max, 50, 0, 80)
     TH1D_ENTRY(pt_b, 50, 0, 50)
     TH2D_ENTRY_EX(h_2D, 50, 0, 50, 50, 0, 1, "P_{T} [GeV]", "dR", false, 1, true )
 };

class GenStudy {
public:
    GenStudy(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.outputFile())), anaData(output)
        {
            GenEvent::intializeNames(args.particleNameFile(), args.particleTypeFile());
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

            const EventIdentifier EventId(event.run, event.lumi, event.evt);
            const EventIdentifier EventIdTest(1, 385, 384944);
            if(!(EventId == EventIdTest)) continue;

            if(event.eventEnergyScale != 0) continue;

            GenEvent genEvent(event);
            // genEvent.Print();
            const GenParticleSet higgsPair = genEvent.GetParticles(particles::ParticleCode::higgs);
            const GenParticleSet B_star_zero = genEvent.GetParticles(513);
            const GenParticleSet Lambda_b0 = genEvent.GetParticles(5122);
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

            int n_tau_matches = 0;
            int n_b_matches_v1 = 0;
            int n_b_matches_v2 = 0;
            for (size_t tau_index = 0; tau_index < h_tautau->daughters.size(); tau_index++) {
                std::vector<const GenParticle*> visible_daughters;
                auto visibleMomentum = GetFinalStateMomentum(*h_tautau->daughters.at(tau_index), visible_daughters, true, false);
                for (size_t reco_tau_index = 0; reco_tau_index < 2; reco_tau_index++) {
                    auto reco_tau_momentum = reco_tau_index == 0 ? event.p4_1 : event.p4_2;
                    if(ROOT::Math::VectorUtil::DeltaR(visibleMomentum, reco_tau_momentum) < deltaR_value)
                        ++n_tau_matches;
                }
              //   std::cout << "reco_tau_1: pt= " << event.p4_1.pt() << " eta= " << event.p4_1.eta() << " phi= " << event.p4_1.phi()
              //             << event.p4_1.mass() << '\n';
              // std::cout << "reco_tau_2: " << "pt= " << event.p4_2.pt() << "eta= " << event.p4_2.eta() << "phi= " << event.p4_2.phi()
              //           << event.p4_2.mass() << '\n';
            }
            for (size_t b_index = 0; b_index < h_bb->daughters.size(); b_index++) {
                std::vector<const GenParticle*> visible_daughters;
                auto visibleMomentum = GetFinalStateMomentum(*h_bb->daughters.at(b_index), visible_daughters, true, false);
                auto dR_1 = ROOT::Math::VectorUtil::DeltaR(h_bb->daughters.at(b_index)->momentum,visibleMomentum);
                anaData.dR_1().Fill(dR_1);
                anaData.pt_b().Fill(h_bb->daughters.at(b_index)->momentum.pt());

                // std::cout << "h_bb->daughters_" << b_index << ": pt= " << h_bb->daughters.at(b_index)->momentum.pt() << " eta= "
                //           << h_bb->daughters.at(b_index)->momentum.eta() << " phi= " << h_bb->daughters.at(b_index)->momentum.phi() << " m="
                //           << h_bb->daughters.at(b_index)->momentum.mass()  << " index= " <<h_bb->daughters.at(b_index)->index << '\n';

                size_t max_i = 0;
                for (size_t i = 0; i < visible_daughters.size(); i++) {
                    auto dR_2 = ROOT::Math::VectorUtil::DeltaR(visible_daughters.at(i)->momentum, visibleMomentum);
                    anaData.dR_2().Fill(dR_2);
                    anaData.pt().Fill(visible_daughters.at(i)->momentum.pt());
                    anaData.h_2D().Fill(visible_daughters.at(i)->momentum.pt(), dR_2);
                    if(visible_daughters.at(i)->momentum.pt() > visible_daughters.at(max_i)->momentum.pt())
                        max_i = i;
                }
                auto dR_3 = ROOT::Math::VectorUtil::DeltaR(visible_daughters.at(max_i)->momentum, h_bb->daughters.at(b_index)->momentum);
                anaData.dR_max_pt_vs_original_momentum().Fill(dR_3);
                anaData.pt_max().Fill(visible_daughters.at(max_i)->momentum.pt());

                LorentzVectorXYZ visibleMomentumV2;
                for (size_t i = 0; i < visible_daughters.size(); i++) {
                    if(ROOT::Math::VectorUtil::DeltaR(visible_daughters.at(max_i)->momentum, visible_daughters.at(i)->momentum) < 0.4)
                        visibleMomentumV2 +=  visible_daughters.at(i)->momentum;

                    // std::cout << "visible_daughters_" << i << ": pt= " << visible_daughters.at(i)->momentum.pt() << " eta= "
                    //           << visible_daughters.at(i)->momentum.eta() << " phi= " << visible_daughters.at(i)->momentum.phi() << " m="
                    //           << visible_daughters.at(i)->momentum.mass() << " pdg= " << visible_daughters.at(i)->pdg<< '\n';
                }
                // std::cout << "visibleMomentumV2: " << ": pt= " << visibleMomentumV2.pt() << " eta= "
                //           << visibleMomentumV2.eta() << " phi= " << visibleMomentumV2.phi() << " m="
                //           << visibleMomentumV2.mass() << '\n';

                // std::cout << "visibleMomentumV1: " << ": pt= " << visibleMomentum.pt() << " eta= "
                //           << visibleMomentum.eta() << " phi= " << visibleMomentum.phi() << " m="
                //           << visibleMomentum.mass() << '\n';

                auto dR_VM2 = ROOT::Math::VectorUtil::DeltaR(h_bb->daughters.at(b_index)->momentum,visibleMomentumV2);
                anaData.dR_orignial_b_vs_VM().Fill(dR_VM2);

                for (size_t reco_b_index = 0; reco_b_index < event.jets_p4.size(); reco_b_index++) {
                    if(ROOT::Math::VectorUtil::DeltaR(visibleMomentum, event.jets_p4.at(reco_b_index)) < deltaR_value)
                        ++n_b_matches_v1;
                    if(ROOT::Math::VectorUtil::DeltaR(visibleMomentumV2, event.jets_p4.at(reco_b_index)) < deltaR_value)
                        ++n_b_matches_v2;

                    // std::cout << "reco_b_" << reco_b_index << ": pt= " << event.jets_p4.at(reco_b_index).pt() << " eta= "
                    //           << event.jets_p4.at(reco_b_index).eta() << " phi= " << event.jets_p4.at(reco_b_index).phi()
                    //           << " m= "<< event.jets_p4.at(reco_b_index).mass() << '\n';
                }

            }
            // if(B_star_zero.size() == 0)
            //     std::cout << event.evt << ", "  << event.lumi << ", " << event.run  << '\n';
            anaData.number_of_b().Fill(h_bb->daughters.size());
            anaData.h_tautau_matches().Fill(n_tau_matches);
            anaData.h_bb_matchesV1().Fill(n_b_matches_v1);
            anaData.h_bb_matchesV2().Fill(n_b_matches_v2);
            anaData.B_star_zero().Fill(B_star_zero.size());
            if(B_star_zero.size() != 2)
                anaData.num("Lambda_b0").Fill(Lambda_b0.size());
        }
    }

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    FileMerger anaData;
};

} // namespace analysis
PROGRAM_MAIN(analysis::GenStudy, Arguments)
