/*! Definition of data and event categories used in HH->bbTauTau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/GenParticle.h"
#include "h-tautau/Analysis/include/HHGenEvent.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Analysis/src/EventInfo.cpp"
#include "AnalysisTools/Print/include/PdfPrinter.h"
#include <Math/VectorUtil.h>

#include "TEfficiency.h"
#include "TStyle.h"
#include <TLegend.h>
#include <TCanvas.h>
#include <functional>
#include <iostream>
#include <algorithm>


struct Arguments {
    REQ_ARG(std::string, inputPath);
    REQ_ARG(analysis::Channel, channel);
    REQ_ARG(std::string, outputFile);
    REQ_ARG(std::string, particleNameTypeFile);
    OPT_ARG(bool, debug, false);
    OPT_ARG(std::string, eventIdBranches, "run:lumi:evt");

};
namespace analysis {

using Event = ntuple::Event;
using EventPtr = std::shared_ptr<Event>;

class GenStudyHist : public root_ext::AnalyzerData {
 public:
     using AnalyzerData::AnalyzerData;

    TH1D_ENTRY_EX(h_tautau_matches, 10, -0.5, 9.5, "number of matches", "events", true, 1, false, true)
    TH2D_ENTRY_EX(n_baryons_vs_n_mesons, 10, -0.5, 9.5, 10, -0.5, 9.5, "number of baryons", "number of mesons", false, 1, true)
    TH2D_ENTRY_EX(bm_vs_jets, 10, -0.5, 9.5, 10, -0.5, 9.5, "number of jets", "number of bm", false, 1, true)

    TH1D_ENTRY(double_match, 10, -0.5, 9.5)
    TH1D_ENTRY(double_match_gen, 10, -0.5, 9.5)

    TH2D_ENTRY_EX(leg_matches,3 , -0.5, 2.5, 3, -0.5, 2.5, "matches #tau_{1}", "matches #tau_{2}", false, 1, true)

    TH1D_ENTRY(jet_matches, 10, -0.5, 9.5)
    TH1D_ENTRY(n, 8, -0.5, 7.5)

    TH1D_ENTRY_EX(Higgs_Pt, 50, 0, 1600, "Pt [GeV]", "events", true, 1, false, true)
    TH1D_ENTRY_EX(jets_energy, 80, 0, 20, "Energy [GeV]", "events", true, 1, false, true)
    TH1D_ENTRY_EX(jets_momentum, 80, 0, 25, "Pt [GeV]", "events", true, 1, false, true)

    TH1D_ENTRY_EX(pt_charged_vis_tau, 50, 0, 500, "Pt [GeV]", "events", true, 1, false, true)

    static std::vector<double> x_bins()
    {
        std::vector<double> x_bins = {25,30,35,40,45,50,60,70,80,90,100,120,160,200,300,400,600,800,1000};
        return x_bins;
    }

    TH1D_ENTRY_CUSTOM_EX(pt_vis_gen_tau, x_bins(), "Pt vis [GeV]", "Efficiency", true, 1, false, true)

    TH1D_ENTRY_EX(inv_mass_jets, 50, 0, 150, "M_{jets} [GeV]", "events", true, 1, false, true)

    TH1D_ENTRY_EX(pt_tau, 50, 0, 500, "Pt [GeV]", "events", true, 1, false, true)
    TH1D_ENTRY_EX(eta_tau, 60, -3, 3,"#eta ", "events", true, 1, false, true)

    TH1D_ENTRY_EX(phi_tau, 60, -3, 3,"#phi ", "events", true, 1, false, true)
    TH1D_ENTRY_EX(deta_tau, 120, -0.05, 0.05,"#Delta#eta (#tau_{reco}, #tau_{gen}) ", "events", true, 1, false, true)
    TH1D_ENTRY_EX(deta_tau_win, 120, -0.05, 0.05,"#Delta#eta (#tau_{reco}, #tau_{gen}) ", "events", true, 1, false, true)

    TH1D_ENTRY_EX(dphi_tau, 120, -7, 7,"#Delta#phi (#tau_{reco}, #tau_{gen})  ", "events", true, 1, false, true)
    TH1D_ENTRY_EX(dphi_tau_win, 120, -0.25, 0.25,"#Delta#phi (#tau_{reco}, #tau_{gen})  ", "events", true, 1, false, true)

    TH2D_ENTRY_EX(delta_eta_vs_delta_phi, 25, -5, 5, 25, -5, 5, "#Delta#eta (#tau_{reco}, #tau_{gen}) ", "#Delta#phi (#tau_{reco}, #tau_{gen}) ", false, 1, true)
 };

using GenJet = std::vector<const GenParticle*>;

class GenStudy {
public:
    GenStudy(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.outputFile())), anaData(output), canvas("","", 600, 600)
        {
            GenEvent::InitializeParticleDataTable(args.particleNameTypeFile());

            gStyle->SetOptStat(0);
            canvas.Print((args.outputFile() + ".pdf[").c_str());
            canvas.Draw();
        }

    void CreateEfficiency(const TH1& passed, const TH1& total, const std::string& prefix, const std::string& channel,
                              const std::string& hist_name, bool print_info = false, const std::string& info_name = "")
    {
        constexpr static double oneSigma = 0.682689492137;

        if(!TEfficiency::CheckConsistency(passed, total))
            throw exception("passed TEfficiency objects do not have consistent bin contents.");
        TEfficiency eff(passed, total);

        eff.SetConfidenceLevel(oneSigma);
        eff.SetStatisticOption(TEfficiency::kFCP);

        if(print_info)
            std::cout << info_name <<"=" << eff.GetEfficiency(2) << " error_up=" << eff.GetEfficiencyErrorUp(2)
                      << " error_low=" << eff.GetEfficiencyErrorLow(2)  << std::endl;
        else{
            std::ostringstream ss_name;
            ss_name  << channel << "_" << hist_name;
            std::string name = ss_name.str();

            eff.SetTitle(name.c_str());
            eff.SetDirectory(output.get());
            eff.Write();
            root_ext::WriteObject(eff, output.get(), name);

            canvas.Print((args.outputFile()+".pdf").c_str(), ("Title:"+name).c_str());
            canvas.Clear();
        }
    }

    static void CreateEfficiencies(std::vector<TEfficiency> eff_plots){
        eff_plots.at(0).Draw();
        for(size_t plot_index = 1; plot_index < eff_plots.size(); ++plot_index){
            eff_plots.at(plot_index).Draw("SAME");
        }
    }

    static Channel genChannel(const HHGenEvent& hh_event)
    {
        if((hh_event.tau_decay.at(0) == TauGenDecayMode::Electron &&  hh_event.tau_decay.at(1)  == TauGenDecayMode::Hadrons) ||
                (hh_event.tau_decay.at(1) == TauGenDecayMode::Electron &&  hh_event.tau_decay.at(0)  == TauGenDecayMode::Hadrons))
            return Channel::ETau;

        else if((hh_event.tau_decay.at(0) == TauGenDecayMode::Muon &&  hh_event.tau_decay.at(1)  == TauGenDecayMode::Hadrons) ||
                    (hh_event.tau_decay.at(1) == TauGenDecayMode::Muon &&  hh_event.tau_decay.at(0)  == TauGenDecayMode::Hadrons))
            return Channel::MuTau;

        else if((hh_event.tau_decay.at(0) == TauGenDecayMode::Hadrons &&  hh_event.tau_decay.at(1)  == TauGenDecayMode::Hadrons))
            return Channel::TauTau;
//        else
//            throw exception ("Unsupported channel.");
    }

    static bool isTauInsideAcceptance(const HHGenEvent& hh_event, Channel channel)
    {
        static const std::map< Channel, std::vector<std::vector<int>> > channel_pt_cuts = {
            { Channel::ETau, { {34,20}, {26,35} } },
            { Channel::MuTau, { {25,20}, {21,32} } },
            { Channel::TauTau,{ {40,40} } } };

        for(size_t n = 0; n < channel_pt_cuts.at(channel).size(); ++n){
            auto pt_cuts= channel_pt_cuts.at(channel).at(n);
            if(hh_event.vis_tau[0].Pt() < pt_cuts.at(0) || hh_event.vis_tau[1].Pt() < pt_cuts.at(1))
                return false;
        }

        for(size_t tau_index = 0; tau_index < hh_event.h_tautau->daughters.size(); ++tau_index){
            double eta_cut = (hh_event.tau_decay.at(tau_index) == TauGenDecayMode::Electron ||
                              hh_event.tau_decay.at(tau_index) == TauGenDecayMode::Muon) ? 2.1 : 2.3;

            if(std::abs(hh_event.vis_tau[tau_index].Eta()) > eta_cut)
                return false;
        }
        if(hasDeltaRMatch(hh_event.vis_tau[0], hh_event.vis_tau[1], 0.1))
            return false;

        return true;
    }

    static bool isBInsideAcceptance(const HHGenEvent& hh_event)
    {
        double pt_cut = 20;
        double eta_cut = 2.4;

        for(size_t jet_index = 0; jet_index < hh_event.b_jets.size(); ++jet_index){
            for(size_t tau_index = 0; tau_index < hh_event.h_tautau->daughters.size(); ++tau_index){
                if(hasDeltaRMatch(hh_event.vis_tau[tau_index], hh_event.b_jets[jet_index], 0.5))
                    return false;
            }
            if(std::abs(hh_event.b_jets[jet_index].Eta()) > eta_cut || hh_event.b_jets[jet_index].Pt() < pt_cut)
                 return false;
        }
        return true;
    }

    void Run()
    {
        auto file = root_ext::OpenRootFile(args.inputPath());
        auto tuple = ntuple::CreateEventTuple(ToString(args.channel()), file.get(), true, ntuple::TreeState::Full);

        //output ntuple with only matched events
        new_output_file = root_ext::CreateRootFile("new_output_file.root");
        auto new_tuple = ntuple::CreateEventTuple("all_events", new_output_file.get(), false, ntuple::TreeState::Full);

        const double deltaR_value = 0.2;
        const double deltaR_jet_value = 0.4;

        for(const auto& event : *tuple) {
           if(static_cast<EventEnergyScale>(event.eventEnergyScale) != EventEnergyScale::Central) continue;

           if(args.debug()) {
               const EventIdentifier EventId(event.run, event.lumi, event.evt);
               const EventIdentifier EventIdTest(args.eventIdBranches());
               if(!(EventId == EventIdTest)) continue;
           }

            GenEvent genEvent(event);
//            genEvent.Print();
            const GenParticleSet higgsPair = genEvent.GetParticles(particles::ParticleCode::higgs, true);
            if(higgsPair.size() != 2)
                throw analysis::exception("Higgs pair size must be 2, insteed has a size of '%1%'.") %higgsPair.size();

            HHGenEvent HH_Gen_Event;

            HH_Gen_Event.h_tautau = nullptr;
            HH_Gen_Event.h_bb = nullptr;

            for(const auto& higgs : higgsPair){
                if(higgs->daughters.size() != 2)
                    throw analysis::exception("Higgs decay products must be 2, insteed has a size of '%1%'.") %higgs->daughters.size();

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

            if(HH_Gen_Event.h_bb->momentum.M() != 125)
                throw analysis::exception("The Higgs (h->bb) mass has a value of '%1%'.") %HH_Gen_Event.h_bb->momentum.M();
            if(HH_Gen_Event.h_tautau->momentum.M() != 125)
                throw analysis::exception("The Higgs (h->tautau) mass has a value of '%1%'.") %HH_Gen_Event.h_tautau->momentum.M();

            //H->tautau

            for (size_t tau_index = 0; tau_index < HH_Gen_Event.h_tautau->daughters.size(); tau_index++) {

                auto tau = HH_Gen_Event.h_tautau->daughters.at(tau_index);

                std::set<const GenParticle*> daughters_set;
                genEvent.FindFinalStateDaughters(*tau, daughters_set, particles::neutrinos());
                std::vector<const GenParticle*> daughters_vector (daughters_set.begin(), daughters_set.end());

                std::set<int> pdg_daughters;
                for(size_t daughter_index = 0 ; daughter_index < daughters_vector.size(); ++daughter_index )
                   pdg_daughters.insert(std::abs(daughters_vector.at(daughter_index)->pdg));

                if(pdg_daughters.count(particles::ParticleCode::e))
                    HH_Gen_Event.tau_decay[tau_index] = TauGenDecayMode::Electron;
                else if(pdg_daughters.count(particles::ParticleCode::mu))
                   HH_Gen_Event.tau_decay[tau_index] = TauGenDecayMode::Muon;
                else
                    HH_Gen_Event.tau_decay[tau_index] = TauGenDecayMode::Hadrons;
            }

            //Channel control
            if(genChannel(HH_Gen_Event) == args.channel())
                anaData.n("correct_channel").Fill(1);
            else {
               anaData.n("correct_channel").Fill(0);
               continue;
            }

            //order tau legs according to Pt and leptonic or hadronic tau
            size_t first_daughter_index = 0;
            if((args.channel() == Channel::TauTau && HH_Gen_Event.h_tautau->daughters.at(0)->momentum.Pt() < HH_Gen_Event.h_tautau->daughters.at(1)->momentum.Pt()) || ((args.channel() == Channel::ETau || args.channel() == Channel::MuTau) && HH_Gen_Event.tau_decay.at(0)  == TauGenDecayMode::Hadrons))
                ++first_daughter_index;
            HH_Gen_Event.tau[0] = HH_Gen_Event.h_tautau->daughters.at(first_daughter_index);
            HH_Gen_Event.tau[1] = HH_Gen_Event.h_tautau->daughters.at((first_daughter_index + 1) % 2);

            //Calculate sum visible momentum of particles
            for (size_t tau_index = 0; tau_index < HH_Gen_Event.h_tautau->daughters.size(); tau_index++) {
                std::vector<const GenParticle*> visible_daughters;
                HH_Gen_Event.vis_tau[tau_index] = genEvent.GetFinalStateMomentum(*HH_Gen_Event.tau[tau_index], visible_daughters, true, false);
                for(size_t vis_daughter_index = 0 ; vis_daughter_index <visible_daughters.size(); ++vis_daughter_index ){
                    if(GenEvent::GetParticleCharge(visible_daughters.at(vis_daughter_index)->pdg) != 0)
                         HH_Gen_Event.vis_charged_tau[tau_index] += visible_daughters.at(vis_daughter_index)->momentum;
                }
            }

            //H->bb

            GenParticleSet baryons_plus_mesons_set;
            genEvent.GetChosenParticlesTypes({particles::ParticleType::baryon, particles::ParticleType::meson},
                                         HH_Gen_Event.h_bb, baryons_plus_mesons_set);

            std::vector<const GenParticle*> baryons_plus_mesons(baryons_plus_mesons_set.begin(), baryons_plus_mesons_set.end());

            std::sort(baryons_plus_mesons.begin(), baryons_plus_mesons.end(), [](const GenParticle* a, const GenParticle* b) {
                return a->momentum.E() > b->momentum.E();
            });

            if(baryons_plus_mesons.size() < 2)
                throw analysis::exception("There're less than two barions or mesons.");

            const auto gen_jets = CreateGenJets(baryons_plus_mesons, 0.4);

            //control for one jet cases
            anaData.n("jets").Fill(gen_jets.size());
            anaData.n("one_jets").Fill(gen_jets.size() == 1);

            if(gen_jets.size() < 2) continue;

            // Sum momentum of all the gen b jets created
            for(size_t jet_index = 0; jet_index < gen_jets.size(); ++jet_index){
                for(size_t particle_index = 2; particle_index < gen_jets.at(jet_index).size(); ++particle_index)
                    HH_Gen_Event.h_bb_vis_all += gen_jets.at(jet_index).at(particle_index)->momentum;
            }
            // Sum momentum of the gen b jets used for the match
            for(size_t jet_index = 0; jet_index < 2; ++jet_index){
                for(size_t particle_index = 0; particle_index < gen_jets.at(jet_index).size(); ++particle_index)
                    HH_Gen_Event.b_jets[jet_index] += gen_jets.at(jet_index).at(particle_index)->momentum;      
                HH_Gen_Event.h_bb_vis += HH_Gen_Event.b_jets[0];
            }
            // Sum momentum of the gen b jets discarded
            for(size_t jet_index = 2; jet_index < gen_jets.size(); ++jet_index){
                for(size_t particle_index = 0; particle_index < gen_jets.at(jet_index).size(); ++particle_index)
                    HH_Gen_Event.b_jets_others[jet_index] += gen_jets.at(jet_index).at(particle_index)->momentum;
                HH_Gen_Event.h_bb_others_vis += HH_Gen_Event.b_jets_others[jet_index];
            }

            if(isBInsideAcceptance(HH_Gen_Event))
                 anaData.n("pass_b_acceptance").Fill(isBInsideAcceptance(HH_Gen_Event));

            if(!isTauInsideAcceptance(HH_Gen_Event, args.channel())) continue;
            anaData.n("pass_taus_acceptance").Fill(isTauInsideAcceptance(HH_Gen_Event, args.channel()));

            if(!isBInsideAcceptance(HH_Gen_Event)) continue;

            anaData.n("pass_b_plus_taus_acceptance").Fill(isBInsideAcceptance(HH_Gen_Event) &&isTauInsideAcceptance(HH_Gen_Event, args.channel()));

            //Matching Taus
            std::set<size_t> tau_matches_total;
            std::map<size_t, std::set<size_t>> tau_matches;

            for (size_t tau_index = 0; tau_index < HH_Gen_Event.tau.size(); tau_index++) {
                const auto& visibleMomentum = HH_Gen_Event.vis_tau[tau_index];
                HH_Gen_Event.h_tautau_vis += HH_Gen_Event.vis_tau[tau_index];
                for (size_t reco_tau_index = 0; reco_tau_index < event.lep_p4.size(); reco_tau_index++) {
                    if(!isCompatible(static_cast<LegType>(event.lep_type.at(reco_tau_index)), HH_Gen_Event.tau_decay[tau_index] )) continue;
//                    if(args.channel() != Channel::TauTau && reco_tau_index != 0 &&
//                        hasDeltaRMatch(event.lep_p4.at(0), event.lep_p4.at(reco_tau_index), 0.2)) continue;
//                        if(HH_Gen_Event.tau_decay[tau_index] != TauGenDecayMode::Hadrons) continue;

                    if(hasDeltaRMatch(visibleMomentum, event.lep_p4.at(reco_tau_index), deltaR_value)){
                        tau_matches_total.insert(reco_tau_index);
                        tau_matches[tau_index].insert(reco_tau_index);
                    }
                    auto dEta = visibleMomentum.Eta() - event.lep_p4.at(reco_tau_index).Eta();
                    auto dPhi = ROOT::Math::VectorUtil::DeltaPhi(visibleMomentum, event.lep_p4.at(reco_tau_index));
                    anaData.deta_tau(HH_Gen_Event.tau_decay[tau_index]).Fill(dEta);
                    anaData.dphi_tau(HH_Gen_Event.tau_decay[tau_index]).Fill(dPhi);
                    anaData.delta_eta_vs_delta_phi(HH_Gen_Event.tau_decay[tau_index]).Fill(dEta, dPhi);
                }
            }

            std::set<size_t> tau_jet_matches_total;
            std::map<size_t, std::set<size_t>> tau_jet_matches;
            for (size_t tau_index = 0; tau_index < HH_Gen_Event.h_tautau->daughters.size(); tau_index++) {
                if(HH_Gen_Event.tau_decay[tau_index] != TauGenDecayMode::Hadrons) continue;
                const auto& visibleMomentum = HH_Gen_Event.vis_tau[tau_index];
                for (size_t reco_tau_index = 0; reco_tau_index < event.jets_p4.size(); reco_tau_index++) {
                    if(hasDeltaRMatch(visibleMomentum, event.jets_p4.at(reco_tau_index), deltaR_jet_value)){
                        tau_jet_matches_total.insert(reco_tau_index);
                        tau_jet_matches[tau_index].insert(reco_tau_index);
                    }
                    auto dEta = visibleMomentum.Eta() - event.jets_p4.at(reco_tau_index).Eta();
                    auto dPhi = ROOT::Math::VectorUtil::DeltaPhi(visibleMomentum, event.jets_p4.at(reco_tau_index));

                    std::ostringstream ss_tau_leg;
                    ss_tau_leg << "jet_"<< HH_Gen_Event.tau_decay[tau_index];
                    std::string tau_leg = ss_tau_leg.str();

                    anaData.deta_tau(tau_leg).Fill(dEta);
                    anaData.dphi_tau(tau_leg).Fill(dPhi);
                    anaData.delta_eta_vs_delta_phi(tau_leg).Fill(dEta, dPhi);
                }
            }

            for (size_t tau_index = 0; tau_index < HH_Gen_Event.h_tautau->daughters.size(); tau_index++) {
                if(tau_matches[tau_index].size() > 1)
                    anaData.double_match_gen().Fill(tau_index);
            }

            if(tau_matches[0].size() == 1 && tau_matches[1].size() == 1
                    && *tau_matches[0].begin() == *tau_matches[1].begin())
                anaData.double_match().Fill(1);

            //Matching b's

            std::set<size_t> b_jet_matches_total;
            std::map<size_t, std::set<size_t>> b_matches;
            for (size_t jet_index = 0; jet_index < HH_Gen_Event.b_jets.size(); jet_index++) {
                for (size_t reco_b_index = 0; reco_b_index < event.jets_p4.size(); reco_b_index++) {
                    if(hasDeltaRMatch(HH_Gen_Event.b_jets[jet_index], event.jets_p4.at(reco_b_index), deltaR_value)){
                        b_jet_matches_total.insert(reco_b_index);
                        b_matches[jet_index].insert(reco_b_index);
                        anaData.n("flavour_reco_jets_matched").Fill(event.jets_hadronFlavour.at(reco_b_index));
                    }
                }
            }

            for (size_t jet_index = 0; jet_index < HH_Gen_Event.b_jets.size(); jet_index++) {
                if(b_matches[jet_index].size() > 1)
                    anaData.double_match_gen("b_jets").Fill(jet_index);
            }

            if(b_matches[0].size() == 1 && b_matches[1].size() == 1
                    && *b_matches[0].begin() == *b_matches[1].begin())
                anaData.double_match("b_jets").Fill(1);


            size_t hadron_flavour_b = 0;
            for (size_t reco_b_index = 0; reco_b_index < event.jets_p4.size(); reco_b_index++) {
                //number of b flavour reco jets
                if(event.jets_hadronFlavour.at(reco_b_index) == 5)
                    ++hadron_flavour_b;
            }
            anaData.n("b_flavour_reco_jets").Fill(hadron_flavour_b);

            anaData.leg_matches().Fill(tau_matches[0].size(), tau_matches[1].size());
            anaData.leg_matches("jets").Fill(tau_jet_matches[0].size(), tau_jet_matches[1].size());

//            if(tau_matches[tau_index].size() < 1){ // double-counting control
            anaData.h_tautau_matches().Fill(tau_matches_total.size());
            anaData.h_tautau_matches("jets").Fill(tau_jet_matches_total.size());

            for(size_t tau_index = 0; tau_index < 2; ++tau_index) {
                if(HH_Gen_Event.tau_decay.at(tau_index) != TauGenDecayMode::Hadrons) continue;
                if(tau_matches[tau_index].size() >= 1)
                    anaData.pt_vis_gen_tau("matched").Fill(HH_Gen_Event.vis_tau[tau_index].Pt());
                if(tau_jet_matches[tau_index].size() >= 1)
                    anaData.pt_vis_gen_tau("matched_jets").Fill(HH_Gen_Event.vis_tau[tau_index].Pt());

                anaData.pt_vis_gen_tau("all").Fill(HH_Gen_Event.vis_tau[tau_index].Pt());

            }
            if(gen_jets.size() ==1)
                anaData.Higgs_Pt("merged").Fill(HH_Gen_Event.h_bb->momentum.Pt());
            else
                anaData.Higgs_Pt("normal").Fill(HH_Gen_Event.h_bb->momentum.Pt());

            if(gen_jets.size() > 2) {
                anaData.jets_energy().Fill(HH_Gen_Event.h_bb_vis_all.E());
                anaData.jets_momentum().Fill(HH_Gen_Event.h_bb_vis_all.Pt());
                anaData.jets_energy("relative").Fill(HH_Gen_Event.h_bb_others_vis.E()/HH_Gen_Event.h_bb_vis_all.E());
                anaData.jets_momentum("relative").Fill(HH_Gen_Event.h_bb_others_vis.Pt()/HH_Gen_Event.h_bb_vis_all.Pt());
            }

            anaData.n("baryons_mesons").Fill(baryons_plus_mesons.size());
            anaData.bm_vs_jets().Fill(gen_jets.size(), baryons_plus_mesons.size());

            anaData.n("jets_matches").Fill(b_jet_matches_total.size());

            anaData.inv_mass_jets().Fill(HH_Gen_Event.h_bb_vis.M());
            anaData.inv_mass_jets("all").Fill(HH_Gen_Event.h_bb_vis_all.M());

            anaData.pt_tau("1").Fill(HH_Gen_Event.vis_tau[0].Pt());
            anaData.eta_tau("1").Fill(HH_Gen_Event.vis_tau[0].Eta());
            anaData.pt_tau("2").Fill(HH_Gen_Event.vis_tau[1].Pt());
            anaData.eta_tau("2").Fill(HH_Gen_Event.vis_tau[1].Eta());

            if(tau_matches_total.size() == 1) {
                for(size_t tau_index = 0; tau_index < 2; ++tau_index) {
                    if(tau_matches[tau_index].empty()){
                        anaData.pt_charged_vis_tau("missing").Fill(HH_Gen_Event.vis_charged_tau[tau_index].Pt());
                    }
                }
            }
            if(b_jet_matches_total.size() == 2 && tau_matches_total.size() ==2) {
                anaData.n("2_matches_for_all").Fill(1);
                for (size_t reco_tau_index = 0; reco_tau_index < event.lep_p4.size(); reco_tau_index++) {
                    int matched_gen_tau = -1;
                    for (int gen_tau_index = 0; gen_tau_index < HH_Gen_Event.tau.size(); gen_tau_index++) {
                        if(tau_matches.at(gen_tau_index).count(reco_tau_index)) {
                            matched_gen_tau = gen_tau_index;
                            break;
                        }
                    }
                    (*new_tuple)().lep_genTauIndex.push_back(matched_gen_tau);
                    new_tuple->Fill();
                }
                for (size_t reco_b_index = 0; reco_b_index < event.jets_p4.size(); reco_b_index++) {
                    int matched_gen_b = -1;
                    for (int gen_b_index = 0; gen_b_index < HH_Gen_Event.b_jets.size(); gen_b_index++) {
                        if(b_matches.at(gen_b_index).count(reco_b_index)) {
                            matched_gen_b = gen_b_index;
                            break;
                        }
                    }
                    (*new_tuple)().jets_genJetIndex.push_back(matched_gen_b);
                    new_tuple->Fill();
                }
                (*new_tuple)() = event;
                new_tuple->Fill();
            }
            else
                anaData.n("2_matches_for_all").Fill(0);
        }
        new_tuple->Write();

        CreateEfficiency(anaData.n("pass_taus_acceptance"),anaData.n("correct_channel"), " ", ToString(args.channel()),
                         "", true, "pass tau acceptance" );
        CreateEfficiency(anaData.n("pass_b_acceptance"),anaData.n("correct_channel"), " ", ToString(args.channel()),
                         "", true, "pass b acceptance" );
        CreateEfficiency(anaData.n("pass_b_plus_taus_acceptance"),anaData.n("correct_channel"), " ", ToString(args.channel()),
                         "", true, "pass tau plus b acceptance" );
        CreateEfficiency(anaData.n("pass_b_plus_taus_acceptance"), anaData.n("pass_taus_acceptance"), " ", ToString(args.channel()),
                         "", true, "pass tau plus b acceptance/pass tau acceptance" );

        anaData.h_tautau_matches("jets").SetLineColor(kBlue+2);
        anaData.h_tautau_matches("jets").SetLineWidth(2);
        anaData.h_tautau_matches("jets").Draw();

        anaData.h_tautau_matches().SetLineWidth(2);
        anaData.h_tautau_matches().SetLineColor(kRed+2);
        anaData.h_tautau_matches().Draw("same");

        std::shared_ptr<TLegend> legend(new TLegend (0.62, 0.65, 0.87, 0.80));
        legend->AddEntry(&anaData.h_tautau_matches(), "matches with taus");
        legend->AddEntry(&anaData.h_tautau_matches("jets"), "matches with jets");
        legend->SetTextSize(0.025f);
        legend->Draw();
        canvas.Print((ToString(args.channel()) +"_h_tautau" + "_matches.pdf]").c_str());

        CreateEfficiency(anaData.pt_vis_gen_tau("matched"), anaData.pt_vis_gen_tau("all"), "leptons", "tauTau", "eff match with taus");
        CreateEfficiency(anaData.pt_vis_gen_tau("matched_jets"), anaData.pt_vis_gen_tau("all"), "jets", "tauTau", "eff match with jets");

        canvas.Print((args.outputFile() + ".pdf]").c_str());
    }


private:
    std::vector<GenJet> CreateGenJets(const std::vector<const GenParticle*>& particles, double deltaR)
    {
        std::vector<GenJet> jets;
        std::set<const GenParticle*> taken_particles;
        while(taken_particles.size() < particles.size()) {
            GenJet jet;
            for(const GenParticle* particle : particles){
                if(taken_particles.count(particle)) continue;
                if(jet.empty() || hasDeltaRMatch(particle->momentum, jet.at(0)->momentum, deltaR)) {
                    jet.push_back(particle);
                    taken_particles.insert(particle);
                }
            }
            jets.push_back(jet);
        }
        return jets;
    }

    static bool isCompatible(const LegType& leg_type, const TauGenDecayMode& tau_decay)
    {
        static const std::map<LegType,TauGenDecayMode> reco_gen_decays = { {LegType::e, TauGenDecayMode::Electron},
                                                              {LegType::mu, TauGenDecayMode::Muon},
                                                              {LegType::tau, TauGenDecayMode::Hadrons},
                                                              {LegType::jet, TauGenDecayMode::Hadrons} };
        if(reco_gen_decays.count(leg_type)){
            if(reco_gen_decays.at(leg_type) == tau_decay) return true;
            else
                return false;
        }
        else
            throw exception ("reco leg type '%1%' now allowed") %leg_type ;
    }

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    std::shared_ptr<TFile> new_output_file;
    GenStudyHist anaData;
    TCanvas canvas;
};

} // namespace analysis
PROGRAM_MAIN(analysis::GenStudy, Arguments)
