/*! Definition of data and event categories used in HH->bbTauTau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/GenParticle.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "h-tautau/JetTools/include/BTagger.h"

#include <iostream>
#include <algorithm>
#include "TEfficiency.h"
#include "TStyle.h"
#include <TLegend.h>
#include <TCanvas.h>

struct Arguments {
    REQ_ARG(std::string, inputPath);
    REQ_ARG(analysis::Channel, channel);
    REQ_ARG(std::string, outputFile);
    REQ_ARG(std::string, period);
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

        for(const auto& event : *tuple) {
           // if(static_cast<EventEnergyScale>(event.eventEnergyScale) != EventEnergyScale::Central) continue;

           if(args.debug()) {
               const EventIdentifier EventId(event.run, event.lumi, event.evt);
               const EventIdentifier EventIdTest(args.eventIdBranches());
               if(!(EventId == EventIdTest)) continue;
           }

            size_t matches_deep_flavour = calculateMatchesJets(event, JetOrdering::DeepFlavour);
            calculateMatchesJets(event, JetOrdering::DeepCSV);
            calculateMatchesJets(event, JetOrdering::Pt);
            calculateMatchesJets(event, JetOrdering::CSV);

            size_t matches_deep_tau = calculateMatchesTaus(event, TauIdDiscriminator::byDeepTau2017v2p1VSjet);
            calculateMatchesTaus(event, TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017);
            calculateMatchesTaus(event, TauIdDiscriminator::byCombinedIsolationDeltaBetaCorr3Hits,false);
            calculateMatchesTaus(event, TauIdDiscriminator::byIsolationMVArun2v1DBoldDMwLT2016);

            size_t n_tautau_matches = args.channel() == Channel::TauTau ? 2 : 1;
            if(matches_deep_flavour == 2 && matches_deep_tau == n_tautau_matches)
                anaData.n("both_two_matches").Fill(1);
            else
                anaData.n("both_two_matches").Fill(0);

        }
        printStats(anaData.n("both_two_matches"), 2, "both matches");
        printStats(anaData.jet_matches(ToString(JetOrdering::DeepFlavour)), 3, "DeepFlavour two matches");
        printStats(anaData.lep_matches(ToString(TauIdDiscriminator::byDeepTau2017v2p1VSjet)), 2, "deepTau two matches");

        std::vector<TH1D> jet_histos = {anaData.jet_matches(ToString(JetOrdering::Pt)),anaData.jet_matches(ToString(JetOrdering::CSV)),
                                         anaData.jet_matches(ToString(JetOrdering::DeepCSV)), anaData.jet_matches(ToString(JetOrdering::DeepFlavour))};
        drawHistoTogether(jet_histos, {ToString(JetOrdering::Pt), ToString(JetOrdering::CSV), ToString(JetOrdering::DeepCSV), ToString(JetOrdering::DeepFlavour)} );

        std::vector<TH1D> lep_histos = {anaData.lep_matches(ToString(TauIdDiscriminator::byDeepTau2017v2p1VSjet)),
                                        anaData.lep_matches(ToString(TauIdDiscriminator::byCombinedIsolationDeltaBetaCorr3Hits)),
                                        anaData.lep_matches(ToString(TauIdDiscriminator::byIsolationMVArun2v1DBoldDMwLT2016)),
                                        anaData.lep_matches(ToString(TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017))};
        drawHistoTogether(lep_histos, {"byDeepTau", "CombinedIsolation", "IsolationMVArun2v1", "IsolationMVArun2017v2"});

        canvas.Print(("match_plots_" + ToString(args.channel()) + ".pdf]").c_str(), "Title:jet_ordering");
    }

private:

    void drawHistoTogether(std::vector<TH1D>& histos, const std::vector<std::string>& histo_names)
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
        auto ordered_jets = orderJets(event, jet_ordering);
        const size_t n_max = std::min<size_t>(ordered_jets.size(), 2);
        size_t count = 0;
        for(size_t n = 0; n < n_max; ++n) {
            const size_t index = ordered_jets.at(n).index;
            if(event.jets_genJetIndex.at(index) != -1) ++count;
        }
        anaData.jet_matches(jet_ordering).Fill(count);
        return count;
    }

    size_t calculateMatchesTaus(const Event& event, TauIdDiscriminator tau_id_discriminator, bool is_descending_order = true)
    {
        std::vector<size_t> ordered_indexes = orderTaus(event, tau_id_discriminator, is_descending_order);
        //Cases for eTau & muTau with zero taus remaning after applied cut

        const size_t n_expected = args.channel() == Channel::TauTau ? 2 : 1;
        const size_t n_max = std::min(ordered_indexes.size(), n_expected);
        size_t count = 0;
        for(size_t n = 0; n < n_max; ++n) {
            const size_t index = ordered_indexes.at(n);
            if(event.lep_genTauIndex.at(index) != -1) ++count;
        }
        anaData.lep_matches(tau_id_discriminator).Fill(count);
        return count;
    }

private:

    std::vector<analysis::jet_ordering::JetInfo<LorentzVector>> orderJets(const Event& event, JetOrdering jet_ordering)
    {
        std::vector<analysis::jet_ordering::JetInfo<LorentzVector>> jet_info_vector;
        BTagger bTagger(Parse<analysis::Period>(args.period()), jet_ordering);

        UncertaintySource unc_source = analysis::UncertaintySource::None;
        UncertaintyScale unc_scale = analysis::UncertaintyScale::Central;
        for(size_t jet_index = 0; jet_index < event.jets_p4.size(); ++jet_index)
            jet_info_vector.emplace_back(LorentzVectorE(event.jets_p4.at(jet_index)), jet_index, bTagger.BTag(event,
                jet_index, unc_source, unc_scale, false));

        return jet_ordering::OrderJets(jet_info_vector, true, bTagger.PtCut(), bTagger.EtaCut());
    }

    std::vector<size_t> orderTaus(const Event& event, TauIdDiscriminator tau_id_discriminator, bool is_descending_order = true)
    {
        std::vector<float> raw_values;
        for(size_t lep_index = 0; lep_index < event.lep_p4.size(); ++lep_index){
           if (static_cast<LegType>(event.lep_type.at(lep_index)) == analysis::LegType::tau ) {
               ntuple::TupleLepton lepton(event, lep_index);

               const bool is_good_lepton = lepton.Passed(TauIdDiscriminator::byDeepTau2017v2p1VSe,DiscriminatorWP::VVVLoose) &&
                       lepton.Passed(TauIdDiscriminator::byDeepTau2017v2p1VSmu,DiscriminatorWP::VLoose);

               const float raw = is_good_lepton ? lepton.GetRawValue(tau_id_discriminator) : -1;
               raw_values.push_back(raw);
            }
           else
                raw_values.push_back(-1);
        }

        std::vector<size_t> index_sorted_had_taus;
        std::vector<size_t> sorted = sort_indexes(raw_values, is_descending_order);

        for (const auto& index : sorted){
            if(static_cast<LegType>(event.lep_type.at(index)) == LegType::tau && raw_values.at(index) >= 0)
                index_sorted_had_taus.push_back(index);
        }
        return index_sorted_had_taus;
    }

    std::vector<size_t> sort_indexes(const std::vector<float> &v, bool is_descending_order = true) {
      // initialize original index locations
      std::vector<size_t> idx(v.size());
      std::iota(idx.begin(), idx.end(), 0);

      // sort indexes based on comparing values in v
      std::sort(idx.begin(), idx.end(),
           [&v,is_descending_order](size_t i1, size_t i2) {
               return is_descending_order ? v[i1] > v[i2] : v[i1] < v[i2];
      });

      return idx;
    }

    void printStats(const TH1D& histo, int binIndex, std::string tag)
    {
        double histoContent = histo.GetEntries();
        double binContent = histo.GetBinContent(binIndex);
        double ratio = binContent / histoContent;
        std::cout << "using bin #" << binIndex << ", rate of " << tag << " : " << binContent << "/" <<  histoContent
                  << " = " << ratio << std::endl;
    }

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    SignalPurityStudyHist anaData;
    TCanvas canvas;
};
} // namespace analysis
PROGRAM_MAIN(analysis::SignalPurityStudy, Arguments)
