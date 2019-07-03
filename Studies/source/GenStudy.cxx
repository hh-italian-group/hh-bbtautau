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

};

namespace analysis {

using Event = ntuple::Event;
using EventPtr = std::shared_ptr<Event>;

class FileMerger : public root_ext::AnalyzerData {
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
    TH1D_ENTRY_EX(jets_energy, 80, 0, 1.5, "Energy [GeV]", "events", true, 1, false, true)
    TH1D_ENTRY_EX(jets_momentum, 80, 0, 1.5, "Pt [GeV]", "events", true, 1, false, true)

    TH1D_ENTRY_EX(pt_charged_vis_tau, 50, 0, 500, "Pt [GeV]", "events", true, 1, false, true)

    const std::vector<double> x_bins = {25,30,35,40,45,50,60,70,80,90,100,120,160,200,300, 400, 600, 800,1000};
    TH1D_ENTRY_CUSTOM_EX(pt_vis_gen_tau, x_bins, "Pt [GeV]", "events", true, 1, false, true)

    TH1D_ENTRY_EX(inv_mass_jets, 50, 0, 150, "M_{jets} [GeV]", "events", true, 1, false, true)

    TH1D_ENTRY_EX(pt_tau, 50, 0, 500, "Pt [GeV]", "events", true, 1, false, true)
    TH1D_ENTRY_EX(eta_tau, 60, -3, 3,"#eta ", "events", true, 1, false, true)
    TH1D_ENTRY_EX(phi_tau, 60, -3, 3,"#phi ", "events", true, 1, false, true)
    TH1D_ENTRY_EX(deta_tau, 25, -5, 5,"#Delta#eta (#tau_{reco}, #tau_{gen}) ", "events", true, 1, false, true)
    TH1D_ENTRY_EX(dphi_tau, 20, -4, 4,"#Delta#phi (#tau_{reco}, #tau_{gen})  ", "events", true, 1, false, true)
    TH2D_ENTRY_EX(delta_eta_vs_delta_phi, 25, -5, 5, 25, -5, 5, "#Delta#eta (#tau_{reco}, #tau_{gen}) ", "#Delta#phi (#tau_{reco}, #tau_{gen}) ", false, 1, true)
 };

using GenJet = std::vector<const GenParticle*>;

class GenStudy {
public:
    GenStudy(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.outputFile())), anaData(output), canvas("","", 600, 600)
        {
            GenEvent::intializeNames(args.particleNameTypeFile());

            gStyle->SetOptStat(0);
            canvas.Print((args.outputFile() + ".pdf[").c_str());
            canvas.Draw();
        }

    void CreateEfficiency(const TH1& passed, const TH1& total, const std::string& prefix, const std::string& channel,
                              const std::string& hist_name, bool print_info = false, const std::string& info_name = "")
    {
            constexpr static double range_sf = 1.1;
            constexpr static double oneSigma = 0.682689492137;

            TFile* pFile = new TFile(("eff_"+prefix+"_"+ToString(args.channel())+".root").c_str(),"recreate");

            TH1D empty_hist("base", "", passed.GetNbinsX(),30, passed.GetBinLowEdge(passed.GetNbinsX()+1));

            if(!TEfficiency::CheckConsistency(passed, total))
                throw exception("passed TEfficiency objects do not have consistent bin contents");
            TEfficiency eff(passed, total);

            eff.SetConfidenceLevel(oneSigma);
            eff.SetStatisticOption(TEfficiency::kFCP);

            if(print_info)
                std::cout << info_name <<"=" << eff.GetEfficiency(2) << " error_up=" << eff.GetEfficiencyErrorUp(2)
                          << " error_low=" << eff.GetEfficiencyErrorLow(2) << std::endl;
            else{

                std::ostringstream ss_name;
                ss_name  << channel << " " << hist_name;
                std::string name = ss_name.str();

                eff.SetTitle(name.c_str());
                canvas.Clear();
                const std::vector<int> x_bins = {25,30,35,40,45,50,60,70,80,90,100,120,160,200,300,400,600,800,1000};
                double min = 1, max = 0;
                for(int n = 1; n <= passed.GetNbinsX(); ++n) {
                    min = std::min(eff.GetEfficiency(n), min);
                    max = std::max(eff.GetEfficiency(n), max);
                    empty_hist.GetXaxis()->SetBinLabel(n, ToString(x_bins.at(n)).c_str());
                    empty_hist.GetXaxis()->SetTitle(("Pt vis [GeV]"));
                    empty_hist.GetYaxis()->SetTitle(("Efficiency"));
                    empty_hist.GetYaxis()->SetTitleOffset(1.5f);
                    empty_hist.SetLabelSize(0.04f);
                }

                min = min > 0.2 ? min / range_sf : 0;
                empty_hist.GetYaxis()->SetRangeUser(min, max * range_sf);
                empty_hist.Draw();
                eff.Draw("SAME P");
                empty_hist.Write();
                canvas.SetLogx();

                empty_hist.SetTitle(name.c_str());
                eff.SetDirectory(pFile);
                empty_hist.SetDirectory(pFile);
                eff.Write();
                pFile->Write();

                canvas.Print((args.outputFile()+".pdf").c_str(), ("Title:"+name).c_str());
                canvas.Clear();
            }
    }

    void CreateEfficiencies(std::vector<TEfficiency> eff_plots){
        eff_plots.at(0).Draw();
        for(size_t plot_index = 0; plot_index < eff_plots.size(); ++plot_index){
            eff_plots.at(plot_index).Draw("SAME");
        }
    }

    bool isCorrectChannel(const HHGenEvent& hh_event, Channel channel)
    {
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
        return false;
    }


    bool isTauInsideAcceptance(const HHGenEvent& hh_event, Channel channel)
    {
        //vector[0]=tau_lep_cut, vector[1]=tau_had_cut,
        static const std::map<Channel, std::vector<int>> channel_pt_cuts = {
            { Channel::ETau, {24,20} },
            { Channel::MuTau, {20,20} },
            { Channel::TauTau, {35,35} } };

        for(size_t tau_index = 0; tau_index < hh_event.h_tautau->daughters.size(); ++tau_index){
            double eta_cut = (hh_event.tau_decay.at(tau_index) == GenDecayMode::Electron ||
                              hh_event.tau_decay.at(tau_index) == GenDecayMode::Muon) ? 2.1 : 2.3;

            auto pt_cuts= channel_pt_cuts.at(channel);
            if(std::abs(hh_event.vis_tau[tau_index].Eta()) > eta_cut || hh_event.vis_tau[0].Pt() < pt_cuts.at(0) ||
                    hh_event.vis_tau[1].Pt() < pt_cuts.at(1))
                return false;
        }
        if(HasMatchWithMCObject(hh_event.vis_tau[0], hh_event.vis_tau[1], 0.1))
            return false;

        return true;
    }

    bool isBInsideAcceptance(const HHGenEvent& hh_event)
    {
        double pt_cut = 20;
        double eta_cut = 2.4;

        for(size_t jet_index = 0; jet_index < hh_event.b_jets.size(); ++jet_index){
            for(size_t tau_index = 0; tau_index < hh_event.h_tautau->daughters.size(); ++tau_index){
                if(HasMatchWithMCObject(hh_event.vis_tau[tau_index], hh_event.b_jets[jet_index], 0.5))
                    return false;
            }
        }
        if(std::abs(hh_event.b_jets[0].Eta()) > eta_cut || hh_event.b_jets[0].Pt() < pt_cut ||
                std::abs(hh_event.b_jets[1].Eta()) > eta_cut || hh_event.b_jets[1].Pt() < pt_cut)
             return false;

        return true;
    }

    void Run()
    {
        auto file = root_ext::OpenRootFile(args.inputPath());
        auto tuple = ntuple::CreateEventTuple(ToString(args.channel()), file.get(), true, ntuple::TreeState::Full);
        const double deltaR_value = 0.2;
        const double deltaR_jet_value = 0.4;

        for(const auto& event : *tuple) {
           if(event.eventEnergyScale != 0) continue;

//           const EventIdentifier EventId(event.run, event.lumi, event.evt);
//           const EventIdentifier EventIdTest(1,101,100757);
//           if(!(EventId == EventIdTest)) continue;

            GenEvent genEvent(event);
//            genEvent.Print();
            const GenParticleSet higgsPair = genEvent.GetParticles(particles::ParticleCode::higgs, true);
            if(higgsPair.size() != 2)
                throw analysis::exception("Higgs pair size must be 2, insteed has a size of '%1%'") %higgsPair.size();

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

            if(HH_Gen_Event.h_bb->momentum.M() != 125)
                std::cout << HH_Gen_Event.h_bb->momentum.M() << std::endl;

            //H->tautau

            for (size_t tau_index = 0; tau_index < HH_Gen_Event.h_tautau->daughters.size(); tau_index++) {

               auto tau = HH_Gen_Event.h_tautau->daughters.at(tau_index);
               std::set<int> pdg_daughters;
               for(size_t daughter_index = 0 ; daughter_index < tau->daughters.size(); ++daughter_index )
                   pdg_daughters.insert(std::abs(tau->daughters.at(daughter_index)->pdg));

               if(pdg_daughters.count(particles::ParticleCode::e))
                   HH_Gen_Event.tau_decay[tau_index] = GenDecayMode::Electron;
               else if(pdg_daughters.count(particles::ParticleCode::mu))
                   HH_Gen_Event.tau_decay[tau_index] = GenDecayMode::Muon;
               else
                  HH_Gen_Event.tau_decay[tau_index] = GenDecayMode::Hadrons;
            }

            //Channel control
            std::vector<size_t> gen_reco_channel_check (3);
            if(args.channel() == Channel::ETau){
                if((HH_Gen_Event.tau_decay.at(0) == GenDecayMode::Electron &&  HH_Gen_Event.tau_decay.at(1)  == GenDecayMode::Hadrons) ||
                        (HH_Gen_Event.tau_decay.at(1) == GenDecayMode::Electron &&  HH_Gen_Event.tau_decay.at(0)  == GenDecayMode::Hadrons))
                    ++gen_reco_channel_check.at(0);
                 anaData.n("eTau").Fill(gen_reco_channel_check.at(0));
                 if((HH_Gen_Event.tau_decay.at(0) == GenDecayMode::Muon &&  HH_Gen_Event.tau_decay.at(1)  == GenDecayMode::Hadrons) ||
                        (HH_Gen_Event.tau_decay.at(1) == GenDecayMode::Muon &&  HH_Gen_Event.tau_decay.at(0)  == GenDecayMode::Hadrons))
                     ++gen_reco_channel_check.at(1);
                  anaData.n("eTau_miss_muTau").Fill(gen_reco_channel_check.at(1));
                if((HH_Gen_Event.tau_decay.at(0) == GenDecayMode::Hadrons &&  HH_Gen_Event.tau_decay.at(1)  == GenDecayMode::Hadrons) ||
                        (HH_Gen_Event.tau_decay.at(1) == GenDecayMode::Hadrons &&  HH_Gen_Event.tau_decay.at(0)  == GenDecayMode::Hadrons))
                    ++gen_reco_channel_check.at(2);
                 anaData.n("eTau_miss_tauTau").Fill(gen_reco_channel_check.at(2));
            }
            else if(args.channel() == Channel::MuTau){
                if((HH_Gen_Event.tau_decay.at(0) == GenDecayMode::Muon &&  HH_Gen_Event.tau_decay.at(1)  == GenDecayMode::Hadrons) ||
                        (HH_Gen_Event.tau_decay.at(1) == GenDecayMode::Muon &&  HH_Gen_Event.tau_decay.at(0)  == GenDecayMode::Hadrons))
                    ++gen_reco_channel_check.at(0);
                 anaData.n("muTau").Fill(gen_reco_channel_check.at(0));
                if((HH_Gen_Event.tau_decay.at(0) == GenDecayMode::Electron &&  HH_Gen_Event.tau_decay.at(1)  == GenDecayMode::Hadrons) ||
                        (HH_Gen_Event.tau_decay.at(1) == GenDecayMode::Electron &&  HH_Gen_Event.tau_decay.at(0)  == GenDecayMode::Hadrons))
                    ++gen_reco_channel_check.at(0);
                 anaData.n("muTau_miss_eTau").Fill(gen_reco_channel_check.at(0));
                if((HH_Gen_Event.tau_decay.at(0) == GenDecayMode::Hadrons &&  HH_Gen_Event.tau_decay.at(1)  == GenDecayMode::Hadrons) ||
                        (HH_Gen_Event.tau_decay.at(1) == GenDecayMode::Hadrons &&  HH_Gen_Event.tau_decay.at(0)  == GenDecayMode::Hadrons))
                    ++gen_reco_channel_check.at(0);
                 anaData.n("muTau_miss_tauTau").Fill(gen_reco_channel_check.at(0));
            }
            else if(args.channel() == Channel::TauTau){
                if((HH_Gen_Event.tau_decay.at(0) == GenDecayMode::Muon &&  HH_Gen_Event.tau_decay.at(1)  == GenDecayMode::Hadrons) ||
                        (HH_Gen_Event.tau_decay.at(1) == GenDecayMode::Muon &&  HH_Gen_Event.tau_decay.at(0)  == GenDecayMode::Hadrons))
                   ++gen_reco_channel_check.at(0);
                anaData.n("tauTau_miss_muTau").Fill(gen_reco_channel_check.at(0));
                if((HH_Gen_Event.tau_decay.at(0) == GenDecayMode::Electron &&  HH_Gen_Event.tau_decay.at(1)  == GenDecayMode::Hadrons) ||
                        (HH_Gen_Event.tau_decay.at(1) == GenDecayMode::Electron &&  HH_Gen_Event.tau_decay.at(0)  == GenDecayMode::Hadrons))
                   ++gen_reco_channel_check.at(1);
                anaData.n("tauTau_miss_eTau").Fill(gen_reco_channel_check.at(1));
                if((HH_Gen_Event.tau_decay.at(0) == GenDecayMode::Hadrons &&  HH_Gen_Event.tau_decay.at(1)  == GenDecayMode::Hadrons) ||
                        (HH_Gen_Event.tau_decay.at(1) == GenDecayMode::Hadrons &&  HH_Gen_Event.tau_decay.at(0)  == GenDecayMode::Hadrons))
                   ++gen_reco_channel_check.at(2);
                anaData.n("tauTau").Fill(gen_reco_channel_check.at(2));

            }

            if(!isCorrectChannel(HH_Gen_Event, args.channel())) continue;
            anaData.n("correct_channel").Fill(isCorrectChannel(HH_Gen_Event, args.channel()));

            //order tau legs according to Pt
            if(args.channel() == Channel::TauTau){
                if(HH_Gen_Event.h_tautau->daughters.at(0)->momentum.Pt() > HH_Gen_Event.h_tautau->daughters.at(1)->momentum.Pt()){
                    HH_Gen_Event.tau[0] = HH_Gen_Event.h_tautau->daughters.at(0);
                    HH_Gen_Event.tau[1] = HH_Gen_Event.h_tautau->daughters.at(1);
                }
                else{
                    HH_Gen_Event.tau[0] = HH_Gen_Event.h_tautau->daughters.at(1);
                    HH_Gen_Event.tau[1] = HH_Gen_Event.h_tautau->daughters.at(0);
                }
            }
            //order tau legs according to leptonic or hadronic tau
            else if(args.channel() == Channel::ETau || args.channel() == Channel::MuTau){
                if((HH_Gen_Event.tau_decay.at(0) == GenDecayMode::Electron &&  HH_Gen_Event.tau_decay.at(1)  == GenDecayMode::Hadrons) ||
                        (HH_Gen_Event.tau_decay.at(0) == GenDecayMode::Muon &&  HH_Gen_Event.tau_decay.at(1)  == GenDecayMode::Hadrons)){
                    HH_Gen_Event.tau[0] = HH_Gen_Event.h_tautau->daughters.at(0);
                    HH_Gen_Event.tau[1] = HH_Gen_Event.h_tautau->daughters.at(1);
                }
                else if((HH_Gen_Event.tau_decay.at(1) == GenDecayMode::Electron &&  HH_Gen_Event.tau_decay.at(0)  == GenDecayMode::Hadrons) ||
                        (HH_Gen_Event.tau_decay.at(1) == GenDecayMode::Muon &&  HH_Gen_Event.tau_decay.at(0)  == GenDecayMode::Hadrons)){
                    HH_Gen_Event.tau[0] = HH_Gen_Event.h_tautau->daughters.at(1);
                    HH_Gen_Event.tau[1] = HH_Gen_Event.h_tautau->daughters.at(0);
                }
             }

            for (size_t tau_index = 0; tau_index < HH_Gen_Event.h_tautau->daughters.size(); tau_index++) {
                std::vector<const GenParticle*> visible_daughters;
                HH_Gen_Event.vis_tau[tau_index] = genEvent.GetFinalStateMomentum(*HH_Gen_Event.tau[tau_index], visible_daughters, true, false);
                auto particle_charge = GenEvent::particleCharge();
                for(size_t vis_daughter_index = 0 ; vis_daughter_index <visible_daughters.size(); ++vis_daughter_index ){
                    if(particle_charge[visible_daughters.at(vis_daughter_index)->pdg] != 0)
                        HH_Gen_Event.vis_charged_tau[tau_index] += visible_daughters.at(vis_daughter_index)->momentum;
                }
            }

            //H->bb

            GenParticleSet baryons_plus_mesons_set;
            genEvent.GetTypesParticles({particles::ParticleCode::ParticleType::baryon, particles::ParticleCode::ParticleType::meson},
                                         HH_Gen_Event.h_bb, baryons_plus_mesons_set);

            std::vector<const GenParticle*> baryons_plus_mesons(baryons_plus_mesons_set.begin(), baryons_plus_mesons_set.end());

            std::sort(baryons_plus_mesons.begin(), baryons_plus_mesons.end(), [](const GenParticle* a, const GenParticle* b) {
                return a->momentum.E() > b->momentum.E();
            });

            if(baryons_plus_mesons.size() < 2)
                throw analysis::exception("There're less than two barions or mesons");

            const auto gen_jets = CreateGenJets(baryons_plus_mesons, 0.4);

            //control for one jet cases
            std::set<size_t> one_jet;
            if(gen_jets.size() == 1)
                one_jet.insert(1);

            anaData.n("jets").Fill(gen_jets.size());

            if(gen_jets.size() < 2) continue;

            for(size_t jet_index = 0; jet_index < 2; ++jet_index){
                for(size_t particle_index = 0; particle_index < gen_jets.at(jet_index).size(); ++particle_index){
                    HH_Gen_Event.b_jets[jet_index] += gen_jets.at(jet_index).at(particle_index)->momentum;
                }
            }

            for(size_t jet_index = 0; jet_index < gen_jets.size(); ++jet_index){
                for(size_t particle_index = 0; particle_index < gen_jets.at(jet_index).size(); ++particle_index){
                    HH_Gen_Event.b_jets_all[jet_index] += gen_jets.at(jet_index).at(particle_index)->momentum;
                }
            }
            if(isTauInsideAcceptance(HH_Gen_Event, args.channel()))
                anaData.n("pass_taus_acceptance").Fill(isTauInsideAcceptance(HH_Gen_Event, args.channel()));
            if(isBInsideAcceptance(HH_Gen_Event))
            anaData.n("pass_b_acceptance").Fill(isBInsideAcceptance(HH_Gen_Event));


            if(!isTauInsideAcceptance(HH_Gen_Event, args.channel())) continue;
            if(!isBInsideAcceptance(HH_Gen_Event)) continue;

            anaData.n("pass_b_plus_taus_acceptance").Fill(isBInsideAcceptance(HH_Gen_Event));

            HH_Gen_Event.h_bb_vis = (HH_Gen_Event.b_jets[0] + HH_Gen_Event.b_jets[1]);

            for(size_t jet_index = 2; jet_index < gen_jets.size(); ++jet_index){
                HH_Gen_Event.h_bb_vis_all += HH_Gen_Event.b_jets_all[jet_index];
            }

            //Matching Taus
            std::set<size_t> tau_matches_total;
            std::map<size_t, std::set<size_t>> tau_matches;

            for (size_t tau_index = 0; tau_index < HH_Gen_Event.tau.size(); tau_index++) {
                for (size_t reco_tau_index = 0; reco_tau_index < event.lep_p4.size(); reco_tau_index++) {
                    if(args.channel() != Channel::TauTau && reco_tau_index != 0 &&
                        HasMatchWithMCObject(event.lep_p4.at(0), event.lep_p4.at(reco_tau_index), 0.2)) continue;
                    const auto& visibleMomentum = HH_Gen_Event.vis_tau[tau_index];
                    HH_Gen_Event.h_tautau_vis += HH_Gen_Event.vis_tau[tau_index];
                    auto dEta = HH_Gen_Event.tau[tau_index]->momentum.Eta() - event.lep_p4.at(reco_tau_index).Eta();
                    auto dPhi = HH_Gen_Event.tau[tau_index]->momentum.Phi() - event.lep_p4.at(reco_tau_index).Phi();
                    anaData.deta_tau().Fill(dEta);
                    anaData.dphi_tau().Fill(dPhi);
                    anaData.delta_eta_vs_delta_phi().Fill(dEta, dPhi);
                    if(HasMatchWithMCObject(visibleMomentum, event.lep_p4.at(reco_tau_index), deltaR_value)){
                        tau_matches_total.insert(reco_tau_index);
                        tau_matches[tau_index].insert(reco_tau_index);
                    }
                }
            }

            std::set<size_t> tau_jet_matches_total;
            std::map<size_t, std::set<size_t>> tau_jet_matches;
            for (size_t tau_index = 0; tau_index < HH_Gen_Event.h_tautau->daughters.size(); tau_index++) {
                for (size_t reco_tau_index = 0; reco_tau_index < event.jets_p4.size(); reco_tau_index++) {
//                    if(args.channel() != Channel::TauTau && reco_tau_index != 0 &&
//                            ROOT::Math::VectorUtil::DeltaR(event.lep_p4.at(0), event.jets_p4.at(reco_tau_index)) < 0.2) continue;
                    const auto& visibleMomentum = HH_Gen_Event.vis_tau[tau_index];
                    if(HasMatchWithMCObject(visibleMomentum, event.jets_p4.at(reco_tau_index), deltaR_jet_value)){
                        tau_jet_matches_total.insert(reco_tau_index);
                        tau_jet_matches[tau_index].insert(reco_tau_index);
                    }
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
            for (size_t jet_index = 0; jet_index < HH_Gen_Event.b_jets.size(); jet_index++) {
                for (size_t reco_b_index = 0; reco_b_index < event.jets_p4.size(); reco_b_index++) {
                    if(HasMatchWithMCObject(HH_Gen_Event.b_jets[jet_index], event.jets_p4.at(reco_b_index), deltaR_value)){
                        b_jet_matches_total.insert(reco_b_index);
                        anaData.n("flavour_reco_jets_matched").Fill(event.jets_hadronFlavour.at(reco_b_index));
                    }
                }
            }

            size_t hadron_flavour_b = 0;
            for (size_t reco_b_index = 0; reco_b_index < event.jets_p4.size(); reco_b_index++) {
                //number of b flavour reco jets
                if(event.jets_hadronFlavour.at(reco_b_index) == 5)
                    ++hadron_flavour_b;
            }
            anaData.n("b_flavour_reco_jets").Fill(hadron_flavour_b);
            anaData.n("one_jet").Fill(one_jet.size());


            anaData.leg_matches().Fill(tau_matches[0].size(), tau_matches[1].size());
            anaData.leg_matches("jets").Fill(tau_jet_matches[0].size(), tau_jet_matches[1].size());

            anaData.h_tautau_matches().Fill(tau_matches_total.size());
            anaData.h_tautau_matches("jets").Fill(tau_jet_matches_total.size());

            for(size_t tau_index = 0; tau_index < 2; ++tau_index) {
                if(HH_Gen_Event.tau_decay.at(tau_index) != GenDecayMode::Hadrons) continue;
                if(tau_matches[tau_index].size() >= 1)
                    anaData.pt_vis_gen_tau("matched").Fill(HH_Gen_Event.vis_tau[tau_index].Pt());
                if(tau_jet_matches[tau_index].size() >= 1)
                    anaData.pt_vis_gen_tau("matched_jets").Fill(HH_Gen_Event.vis_tau[tau_index].Pt());

                anaData.pt_vis_gen_tau("all").Fill(HH_Gen_Event.vis_tau[tau_index].Pt());

            }

            anaData.Higgs_Pt("normal").Fill(HH_Gen_Event.h_bb->momentum.Pt());
            if(gen_jets.size() ==1)
                anaData.Higgs_Pt("merged").Fill(HH_Gen_Event.h_bb->momentum.Pt());

            if(gen_jets.size() > 2) {
                anaData.jets_energy().Fill(HH_Gen_Event.h_bb_vis_all.E());
                anaData.jets_momentum().Fill(HH_Gen_Event.h_bb_vis_all.Pt());
                anaData.jets_energy("relative").Fill(HH_Gen_Event.h_bb_vis_all.E()/HH_Gen_Event.h_bb->momentum.E());
                anaData.jets_momentum("relative").Fill(HH_Gen_Event.h_bb_vis_all.Pt()/HH_Gen_Event.h_bb->momentum.Pt());
            }

            anaData.n("baryons_mesons").Fill(baryons_plus_mesons.size());
            anaData.bm_vs_jets().Fill(gen_jets.size(), baryons_plus_mesons.size());

            anaData.n("jets_matches").Fill(b_jet_matches_total.size());

            if(b_jet_matches_total.size() == 2 && tau_matches_total.size() ==2)
                anaData.n("2_matches_for_all").Fill(1);
            else
                anaData.n("2_matches_for_all").Fill(0);

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
        }
        CreateEfficiency(anaData.n("pass_taus_acceptance"),anaData.n("correct_channel"), " ", ToString(args.channel()),
                         "", true, "pass tau acceptance" );
        CreateEfficiency(anaData.n("pass_b_acceptance"),anaData.n("correct_channel"), " ", ToString(args.channel()),
                         "", true, "pass b acceptance" );
        CreateEfficiency(anaData.n("pass_b_plus_taus_acceptance"),anaData.n("correct_channel"), " ", ToString(args.channel()),
                         "", true, "pass tau plus b acceptance" );
//        CreateEfficiency(anaData.n("pass_taus_acceptance"), anaData.n("pass_b_acceptance"), " ", ToString(args.channel()),
//                         "", true, "pass tau acceptance/pass b acceptance" );

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
                if(jet.empty() || HasMatchWithMCObject(particle->momentum, jet.at(0)->momentum, deltaR)) {
                    jet.push_back(particle);
                    taken_particles.insert(particle);
                }
            }
            jets.push_back(jet);
        }
        return jets;
    }

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    FileMerger anaData;
    TCanvas canvas;
};

} // namespace analysis
PROGRAM_MAIN(analysis::GenStudy, Arguments)
