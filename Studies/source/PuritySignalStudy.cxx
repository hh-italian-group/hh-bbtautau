/*! Definition of data and event categories used in HH->bbTauTau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Print//include/PlotPrimitives.h"
#include "h-tautau/Analysis/include/GenParticle.h"
#include "h-tautau/Analysis/include/HHGenEvent.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Analysis/src/EventInfo.cpp"
#include "AnalysisTools/Print/include/PdfPrinter.h"
#include <Math/VectorUtil.h>
#include <TTree.h>
#include "h-tautau/JetTools/include/BTagger.h"
#include "h-tautau/Analysis/include/SignalObjectSelector.h"
#include "h-tautau/Core/include/TauIdResults.h"
#include <TColor.h>

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
    OPT_ARG(bool, debug, false);
    OPT_ARG(std::string, eventIdBranches, "run:lumi:evt");
};
namespace analysis {

using Event = ntuple::Event;
using EventPtr = std::shared_ptr<Event>;
using JetCollection = std::vector<JetCandidate>;
using TauIdDiscriminator = analysis::TauIdDiscriminator;

class SignalPurityStudyHist : public root_ext::AnalyzerData {
 public:
     using AnalyzerData::AnalyzerData;

    TH1D_ENTRY_EX(jet_matches, 10, -0.5, 9.5, "number of matches", "events", true, 1, false, true)
    TH1D_ENTRY_EX(lep_matches, 10, -0.5, 9.5, "number of matches", "events", true, 1, false, true)
    TH1D_ENTRY(n, 8, -0.5, 7.5)
 };

class SignalPurityStudy {
public:
    SignalPurityStudy(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.outputFile())), anaData(output), canvas("","", 600, 600)
        {
            gStyle->SetOptStat(0);
            gStyle->SetPadLeftMargin(0.3f);
            gROOT->ForceStyle();
            canvas.Print(("match_plots_" + ToString(args.channel()) + ".pdf[").c_str());
            canvas.Draw();
        }

    void Run()
    {
        auto file = root_ext::OpenRootFile(args.inputPath());
        auto tuple = ntuple::CreateEventTuple("all_events", file.get(), true, ntuple::TreeState::Full);

        int count = 0;
        for(const auto& event : *tuple) {
            ++count;
//            if(count > 5) break;
           if(static_cast<EventEnergyScale>(event.eventEnergyScale) != EventEnergyScale::Central) continue;

           if(args.debug()) {
               const EventIdentifier EventId(event.run, event.lumi, event.evt);
               const EventIdentifier EventIdTest(args.eventIdBranches());
               if(!(EventId == EventIdTest)) continue;
           }

            size_t matches_deep_csv =  calculateMatchesJets(event, JetOrdering::DeepCSV);
            size_t matches_deep_flavour = calculateMatchesJets(event, JetOrdering::DeepFlavour);
            calculateMatchesJets(event, JetOrdering::Pt);
            calculateMatchesJets(event, JetOrdering::CSV);

            size_t matches_deep_tau = calculateMatchesTaus(event, TauIdDiscriminator::byDeepTau2017v2VSjet, args.channel());
            size_t matches_iso_mva_run2017v2 = calculateMatchesTaus(event, TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017, args.channel());
            calculateMatchesTaus(event, TauIdDiscriminator::byCombinedIsolationDeltaBetaCorr3Hits, args.channel(),false);
            calculateMatchesTaus(event, TauIdDiscriminator::byIsolationMVArun2v1DBoldDMwLT2016, args.channel());

            int n_tautau_matches = args.channel() == Channel::TauTau ? 2 : 1;
            if(matches_deep_flavour ==2 && matches_deep_tau == n_tautau_matches)
                anaData.n("both_two_matches").Fill(1);
            else
                anaData.n("both_two_matches").Fill(0);

        }
        printStats(anaData.n("both_two_matches"), 2, "both matches");
        printStats(anaData.jet_matches(ToString(JetOrdering::DeepFlavour)), 3, "DeepFlavour two matches");
        printStats(anaData.lep_matches(ToString(TauIdDiscriminator::byDeepTau2017v2VSjet)), 2, "deepTau two matches");

        drawHistoTogether({anaData.jet_matches(ToString(JetOrdering::Pt)),anaData.jet_matches(ToString(JetOrdering::CSV)),
                           anaData.jet_matches(ToString(JetOrdering::DeepCSV)), anaData.jet_matches(ToString(JetOrdering::DeepFlavour)) },
                           {ToString(JetOrdering::Pt), ToString(JetOrdering::CSV), ToString(JetOrdering::DeepCSV), ToString(JetOrdering::DeepFlavour)} );


        drawHistoTogether({anaData.lep_matches(ToString(TauIdDiscriminator::byDeepTau2017v2VSjet)),
                           anaData.lep_matches(ToString(TauIdDiscriminator::byCombinedIsolationDeltaBetaCorr3Hits)),
                           anaData.lep_matches(ToString(TauIdDiscriminator::byIsolationMVArun2v1DBoldDMwLT2016)),
                           anaData.lep_matches(ToString(TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017)) },
                           {"byDeepTau", "CombinedIsolation", "IsolationMVArun2v1", "IsolationMVArun2017v2"});

        canvas.Print(("match_plots_" + ToString(args.channel()) + ".pdf]").c_str(), "Title:jet_ordering");
    }

private:

    void drawHistoTogether(std::vector<TH1D> histos/*, std::vector<Color_t> colors*/, std::vector<std::string> histo_names)
    {
        auto legend = new TLegend(0.60, 0.65, 0.87, 0.80);
        std::vector<Color_t> colors = {kRed+1, kAzure+9, kCyan+3, kOrange+8};

        histos.at(0).SetLineColor(colors.at(0));
        histos.at(0).SetLineWidth(2);
        histos.at(0).GetYaxis()->SetTitleOffset(1.8f);
        histos.at(0).Draw();

        legend->AddEntry(&histos.at(0), ToString(histo_names.at(0)).c_str());
        legend->SetLineWidth(0);
        legend->Draw("same");

        for(size_t histo_index = 1; histo_index < histos.size(); ++histo_index){
            histos.at(histo_index).SetLineColor(colors.at(histo_index));
            histos.at(histo_index).SetLineWidth(2);
            histos.at(histo_index).Draw("same");

            legend->AddEntry(&histos.at(histo_index), ToString(histo_names.at(histo_index)).c_str());
            legend->SetTextSize(0.022f);
            legend->Draw("same");
        }
        gPad->SetLeftMargin(0.15f);
        gPad->Update();

        canvas.Print(("match_plots_" + ToString(args.channel()) + ".pdf").c_str(), "Title:jet_ordering");
    }

    size_t calculateMatchesJets(const Event& event, JetOrdering jet_ordering)
    {
        auto best_two = getBestTwoJetsIndexes(orderJets(event, jet_ordering ));

        // find indexes that refer to the elements with match
        std::vector<size_t> index_matched_jet;
        int counter = 0;
        for(auto i : event.jets_genJetIndex){
            if(i != 10) // 10 is the defult value of this branch on the ntuple
                index_matched_jet.push_back(counter);
            ++counter;
        }

        // find common elements between the indexes that refer to the two best indexes after ordering
        // and the indexes that refer to the elements with match
        auto common_indexes = intersection(index_matched_jet, best_two);
        anaData.jet_matches(ToString(jet_ordering)).Fill(common_indexes.size());

        return common_indexes.size();
    }

    size_t calculateMatchesTaus(const Event& event, TauIdDiscriminator tau_id_discriminator, Channel channel, bool is_descending_order = true)
    {
        std::vector<size_t> ordered_indexes = orderTaus(event, tau_id_discriminator, is_descending_order);
        std::vector<size_t> best_two = {ordered_indexes.at(0)};

        if(args.channel() == Channel::TauTau)
            best_two.push_back(ordered_indexes.at(1));

        // find indexes that refer to the elements with match
        std::vector<size_t> index_matched_tau;
        int counter = 0;
        for(auto i : event.lep_genTauIndex){
            if(i != 10) // 10 is the defult value of this branch on the ntuple
                index_matched_tau.push_back(counter);
            ++counter;
        }

        // find common elements between the indexes that refer to the two best indexes after ordering
        // and the indexes that refer to the elements with match
        auto common_indexes = intersection(index_matched_tau, best_two);

        anaData.lep_matches(ToString(tau_id_discriminator)).Fill(common_indexes.size());

        return common_indexes.size();
    }

    std::vector<size_t> getBestTwoJetsIndexes(std::vector<analysis::jet_ordering::JetInfo<LorentzVector>> input)
    {
        std::vector<size_t> output;
        int ctr = 0;

        for (auto input_i : input) {
            output.push_back(input_i.index);
            ++ctr;
            if (ctr==2)
                break;
        }
       return output;
    }

private:

    std::vector<analysis::jet_ordering::JetInfo<LorentzVector>> orderJets(const Event& event, JetOrdering jet_ordering/*, Period period*/)
    {
        std::vector<analysis::jet_ordering::JetInfo<LorentzVector>> jet_info_vector;
        BTagger bTagger(Period::Run2017,jet_ordering);

        for(size_t jet_index = 0; jet_index < event.jets_p4.size(); ++jet_index)
            jet_info_vector.emplace_back(static_cast<LorentzVectorE>(event.jets_p4.at(jet_index)), jet_index, bTagger.BTag(event,jet_index));


        return jet_ordering::OrderJets(jet_info_vector, false, 20, 2.1);
    }

    std::vector<size_t> orderTaus(const Event& event, TauIdDiscriminator tau_id_discriminator, bool is_descending_order = true)
    {
       std::vector<float> raw_values;
       for(size_t lep_index = 0; lep_index < event.lep_p4.size(); ++lep_index){
           if (static_cast<LegType>(event.lep_type.at(lep_index)) == analysis::LegType::tau ) {
               auto lepton = std::make_shared<ntuple::TupleLepton>(event, lep_index);
               auto raw = static_cast<float>(lepton->GetRawValue(tau_id_discriminator));
               raw_values.push_back(raw);
           }
           else
                raw_values.push_back(0);
       }


       std::vector<size_t> sorted_had_taus;
       std::vector<size_t> sorted = sort_indexes(raw_values, is_descending_order);

       for (const auto& index : sorted){
               if(isHadronicTau(static_cast<LegType>(event.lep_type.at(index))))
                   sorted_had_taus.push_back((index));
       }
       return sorted_had_taus;
    }

    std::vector<size_t> sort_indexes(const std::vector<float> &v, bool is_descending_order = true) {
      // initialize original index locations
      std::vector<size_t> idx(v.size());
      iota(idx.begin(), idx.end(), 0);

      // sort indexes based on comparing values in v
      sort(idx.begin(), idx.end(),
           [&v,is_descending_order](size_t i1, size_t i2) {
               if(!is_descending_order)
                  return v[i1] < v[i2];
               else
                  return v[i1] > v[i2];
      });

      return idx;
    }

    std::vector<size_t> intersection(std::vector<size_t> &v1, std::vector<size_t> &v2)
    {
        std::vector<size_t> v3;

        std::sort(v1.begin(), v1.end());
        std::sort(v2.begin(), v2.end());

        std::set_intersection(v1.begin(),v1.end(),
                              v2.begin(),v2.end(),
                              back_inserter(v3));
        return v3;
    }

    void printStats(const TH1D& histo, int binIndex, std::string tag)
    {
        double histoContent = histo.GetEntries();
        double binContent = histo.GetBinContent(binIndex);
        double ratio = binContent / histoContent;
        std::cout << "using bin #" << binIndex << ", rate of " << tag << " : " << binContent << "/" <<  histoContent
                  << " = " << ratio << std::endl;
    }

    static bool isHadronicTau(const LegType& leg_type)
    {
        std::map<LegType,bool> taus_legs = { {LegType::e, false},
                                                              {LegType::mu, false},
                                                              {LegType::tau, true}};
        if(taus_legs.count(leg_type))
            return taus_legs[leg_type];
        else
            throw exception ("reco leg type '%1%' now allowed") %leg_type ;
    }

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    SignalPurityStudyHist anaData;
    TCanvas canvas;

};

} // namespace analysis
PROGRAM_MAIN(analysis::SignalPurityStudy, Arguments)

