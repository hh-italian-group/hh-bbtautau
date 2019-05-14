#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "hh-bbtautau/Analysis/include/GenParticle.h"
#include "hh-bbtautau/Analysis/include/Particle.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"

#include "AnalysisTools/Core/include/EventIdentifier.h"


struct Arguments {
    REQ_ARG(std::string, inputPath);
    REQ_ARG(std::string, outputFile);
    // REQ_ARG(std::string, treeName);
};

namespace analysis {


using Event = ntuple::Event;
using EventPtr = std::shared_ptr<Event>;

class FileMerger : public root_ext::AnalyzerData {
 public:
     using AnalyzerData::AnalyzerData;
     // using particles::ParticleCode::pdg;

     TH1D_ENTRY(Higgs_pT, 100, 10, 200)
 };

class GenStudy {
public:
    GenStudy(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.outputFile())), anaData(output)
        {}

    void Run()
    {
        auto file = root_ext::OpenRootFile(args.inputPath());
        auto tuple = ntuple::CreateEventTuple("tauTau", file.get(), true, ntuple::TreeState::Full);

        // int n_matches = 0;
        // int n_no_matches = 0;
        for(const auto& event : *tuple) {

            const EventIdentifier EventId(event.run, event.lumi, event.evt);
            const EventIdentifier EventIdTest(1, 384, 383685);
            if(!(EventId == EventIdTest)) continue;
            std::cout << "pt_1 " << event.p4_1.pt() << " , eta: " << event.p4_1.eta()<< '\n';


            GenEvent genEvent(event);
            const GenParticleSet higgsPair = genEvent.GetParticles(particles::ParticleCode::higgs);
            // std::cout << "eventId" <<  << '\n';
            // if(higgsPair.size() == 0)
            //     std::cout << "higgsPair.size(): " << higgsPair.size() << " event.run: "<< event.run << " event.lumi: "<<event.lumi << "eventId: " <<event.evt << '\n';



                // throw analysis::exception("Higgs pair size must be 2."); //Stampare eventId e albero

            // const GenParticle* h_tautau = nullptr;
            // const GenParticle* h_bb = nullptr;
            //
            for(const auto& higgs : higgsPair){
            //     if(higgs->daughters.at(0)->pdg == particles::ParticleCode::tau && higgs->daughters.at(1)->pdg == particles::ParticleCode::tau)
            //         h_tautau = higgs;
            //     else if(higgs->daughters.at(0)->pdg == particles::ParticleCode::b && higgs->daughters.at(1)->pdg == particles::ParticleCode::b)
            //          h_bb = higgs;
            //      // else
            //      //     throw analysis::exception("Higgs pair doesn't decay in a pair of taus or b jets.");
            //
                 // genEvent.PrintChain(higgs);
            // }
            // // if(h_tautau == nullptr || h_bb == nullptr)
            // //     throw analysis::exception("The Higgs pairs are not filled");
            // if(h_tautau->daughters.size() == 0 || h_bb->daughters.size() == 0) continue;
            //
            // for (size_t tau_index = 0; tau_index < h_tautau->daughters.size(); tau_index++) {
            //     std::vector<const GenParticle*> visible_daughters;
                // std::cout << "Nutella?" << '\n';
                // auto visibleMomentum = GetFinalStateMomentum(particle, visible_daughters, true, false);
                // for (size_t reco_tau_index = 0; reco_tau_index < h_tautau->daughters.size(); reco_tau_index++) {
                    // auto reco_tau_pt = reco_tau_index == 1 ? event.p4_1 : event.p4_2;
                    // std::cout << "reco_tau_index: " << reco_tau_index << '\n';
                    // if(HasMatchWithMCObject(*h_tautau->daughters.at(tau_index), visible_daughters, reco_tau_pt, 0.2, true))
                    //     ++n_matches;
                    // else
                    //     ++n_no_matches;
                // }
            }
        }
        // anaData.Higgs_pT("sample_name").Fill(higgs->momentum.Pt(), event.weight_total);
    }

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    FileMerger anaData;
};

} // namespace analysis
PROGRAM_MAIN(analysis::GenStudy, Arguments)
