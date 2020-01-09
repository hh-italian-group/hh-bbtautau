/*! Study of the efficiency of the triggers made to the events at baseline selection.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include "AnalysisTools/Print/include/PdfPrinter.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/PropertyConfigReader.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "h-tautau/Analysis/include/SignalObjectSelector.h"

struct Arguments {
    run::Argument<std::string> output_file{"output_file", "Output pdf file"};
    run::Argument<analysis::Channel> channel{"channel", "channel on which we work"};
    run::Argument<analysis::Period> period{"period", "", analysis::Period::Run2017};
    run::Argument<std::string> sample_type{"sample_type", "", ""};
    run::Argument<std::string> file_cfg_name{"file_cfg_name", "file cfg"};
    run::Argument<std::string> input_path{"input_path", ""};
    run::Argument<analysis::SignalMode> mode{"mode", "signal mode"};
};

namespace analysis {

class EfficiencyStudy{
public:

    struct SampleDesc {
        std::vector<double> points_list;
        std::string file_name;
        std::string x_title;
        analysis::EllipseParameters massWindowParams;
    };

    struct TriggerDesc {

        std::vector<std::string> trigger;
        std::string title;
        root_ext::Color color;
        TauIdDiscriminator tauId;
        DiscriminatorWP tauIdWP;
    };

    EfficiencyStudy(const Arguments& _args) :
        args(_args), canvas("","", 600, 600), signalObjectSelector(args.mode())
    {
        gStyle->SetOptStat(0);
        canvas.SetGrid();
        canvas.Print((args.output_file() + ".pdf[").c_str());
        canvas.Draw();

        sample_desc = LoadConfig(trigger_descs);
    }

    static TEfficiency CreateEfficiency(const TH1D& passed, const TH1& total)
    {
        constexpr static double oneSigma = 0.682689492137;

        if(!TEfficiency::CheckConsistency(passed, total))
            throw exception("passed TEfficiency objects do not have consistent bin contents");

        TEfficiency eff (passed, total);
        eff.SetConfidenceLevel(oneSigma);
        eff.SetStatisticOption(TEfficiency::kFCP);

        return eff;
    }

    static std::shared_ptr<TGraphAsymmErrors> CreateTGraph(const TEfficiency& eff,
                                                           const std::vector<double>& mass_list)
    {
        std::vector<double> eff_value(mass_list.size());
        std::vector<double> eff_errorUp(mass_list.size());
        std::vector<double> eff_errorLow(mass_list.size());
        std::vector<double> eff_error_x(mass_list.size(), 0.);

        for(size_t n = 0; n < mass_list.size(); n++){
             eff_value.at(n) = eff.GetEfficiency(static_cast<int>(n+1));
             eff_errorUp.at(n) = eff.GetEfficiencyErrorUp(static_cast<int>(n+1));
             eff_errorLow.at(n) = eff.GetEfficiencyErrorLow(static_cast<int>(n+1));
         }

        return std::make_shared<TGraphAsymmErrors>(static_cast<int>(mass_list.size()), mass_list.data(),
               eff_value.data(), eff_error_x.data(), eff_error_x.data(),eff_errorLow.data(),
               eff_errorUp.data());
    }

    static std::shared_ptr<TGraph> CreateRatio(const TGraphAsymmErrors& nume, const TGraphAsymmErrors& deno)
    {
        const size_t n_points = static_cast<size_t>(nume.GetN());
        std::vector<double> x, y;
        for(size_t n = 0; n < n_points; n++){
            if(deno.GetY()[n] != 0) {
                x.push_back(deno.GetX()[n]);
                y.push_back(nume.GetY()[n] / deno.GetY()[n]);
            }
        }
        return std::make_shared<TGraph>(x.size(), x.data(), y.data());
    }

    void PrintTGraph(const std::vector<std::shared_ptr<TGraphAsymmErrors>>& gr,
                     const std::vector<std::shared_ptr<TGraph>>& ratioGraphs)
    {
        canvas.cd();

        TPad pad1("pad1", "", 0, 0.18, 1, 1);
        pad1.Draw();
        pad1.cd();
        pad1.SetGridy(1);
        pad1.SetGridx(1);
        pad1.SetLeftMargin(0.15f);

        auto mg = std::make_shared<TMultiGraph>();
        for (size_t i = 0; i < gr.size(); i++) {
            const auto& single_graph = gr.at(i);
            single_graph->SetTitle((trigger_descs.at(i).title).c_str());
            single_graph->SetName("gr");
            single_graph->SetMarkerColor(trigger_descs.at(i).color.GetColor_t());
            single_graph->SetMarkerStyle(7);
            single_graph->SetLineColor(trigger_descs.at(i).color.GetColor_t());
            mg->Add(single_graph.get());
        }

        mg->Draw("ALP");
        mg->SetTitle("");
        mg->GetYaxis()->SetTitle("Efficiency");
        mg->GetYaxis()->SetTitleOffset(1.8f);

        if (trigger_descs.size() > 1){
            mg->GetXaxis()->SetLabelSize(0.25f);
            mg->GetYaxis()->SetLabelSize(0.03f);
        }

        auto legend = pad1.BuildLegend(0.55, 0.7, 0.9, 0.9, "", "lep");
        legend->SetTextSize(0.03f);

        std::shared_ptr<TPad> pad2;
        std::shared_ptr<TMultiGraph> mgRatio;
        if (trigger_descs.size() > 1){

            canvas.cd();
            pad2 = std::make_shared<TPad>("ratio", "ratio", 0, 0.0, 1, 0.252);
            pad2->Draw();
            pad2->cd();
            pad2->SetGridy(1);
            pad2->SetLeftMargin(static_cast<float>(0.15));
            pad2->SetBottomMargin(0.25);

            mgRatio = std::make_shared<TMultiGraph>();

            for (size_t i = 0; i < ratioGraphs.size(); ++i) {
                auto single_RatioGraph = ratioGraphs.at(i);
                single_RatioGraph->SetMarkerStyle(7);
                single_RatioGraph->SetLineColor(trigger_descs.at(i+1).color.GetColor_t());
                single_RatioGraph->SetMarkerColor(trigger_descs.at(i+1).color.GetColor_t());
                mgRatio->Add(single_RatioGraph.get());
            }

            mgRatio->Draw("ALP");

            mgRatio->GetXaxis()->SetTitle(sample_desc.x_title.c_str());
            mgRatio->GetXaxis()->SetTitleOffset(0.75f);
            mgRatio->GetXaxis()->SetTitleSize(0.13f);
            mgRatio->GetXaxis()->SetLabelSize(0.09f);
            mgRatio->GetXaxis()->SetLabelOffset(0.015f);
            mgRatio->GetYaxis()->SetTitle("Ratio");
            mgRatio->GetYaxis()->SetTitleOffset(0.55f);
            mgRatio->GetYaxis()->SetTitleSize(0.1f);
            mgRatio->GetYaxis()->SetLabelSize(0.09f);
            mgRatio->GetYaxis()->SetLabelOffset(0.015f);
        }

        canvas.Print((args.output_file()+".pdf").c_str());
        canvas.Clear();

    }

    void Run()
    {
        const std::vector<double>& mass_list = sample_desc.points_list;
        std::vector<std::string> input_file;

        for (size_t i = 0; i < mass_list.size(); i++) {
            std::ostringstream ss_name;
            ss_name << std::setprecision(4) << mass_list.at(i);
            const std::string name = ss_name.str();

            const std::string file_string_name = boost::replace_all_copy(sample_desc.file_name, "{x}", name);

           std::ostringstream ss_full_input_file;
           ss_full_input_file << args.input_path() << "/" << file_string_name;
           input_file.emplace_back(ss_full_input_file.str());
        }

        TH1D deno("", "",static_cast<int>(mass_list.size()), 0, mass_list.size());
        std::vector<TH1D> numerators(trigger_descs.size());
        for(size_t desc_id = 0; desc_id < trigger_descs.size(); ++desc_id)
            numerators.at(desc_id) = TH1D(deno);

        for(size_t fileID = 0; fileID < input_file.size(); fileID++){
            auto file = root_ext::OpenRootFile(input_file.at(fileID));
            auto summaryTuple = ntuple::CreateSummaryTuple("summary", file.get(), true,
                                                           ntuple::TreeState::Full);
            auto summary = MergeSummaryTuple(*summaryTuple);
            std::shared_ptr<SummaryInfo> summaryInfo(new SummaryInfo(summary,args.channel()));
            auto originalTuple = ntuple::CreateEventTuple(ToString(args.channel()), file.get(), true,
                                                          ntuple::TreeState::Full);
            const Channel channel = args.channel();

            deno.SetBinContent(static_cast<int>(fileID + 1), summary.numberOfProcessedEvents);
            deno.SetBinError(static_cast<int>(fileID + 1), std::sqrt(summary.numberOfProcessedEvents));

            std::vector<size_t> passed(trigger_descs.size(), 0);
            for(const auto& event : *originalTuple) {
                const Period run_period = args.period();
                const JetOrdering jet_ordering = run_period == Period::Run2017
                                               ? JetOrdering::DeepCSV : JetOrdering::CSV;

                boost::optional<EventInfoBase> eventInfo = CreateEventInfo(event,signalObjectSelector,summaryInfo.get(), run_period, jet_ordering);
                if(!eventInfo.is_initialized()) continue;
                if((*eventInfo)->extraelec_veto || (*eventInfo)->extramuon_veto) continue;
                // if(eventInfo->GetEnergyScale() != EventEnergyScale::Central) continue;
                if(!eventInfo->HasBjetPair()) continue;
                if(!sample_desc.massWindowParams.IsInside(eventInfo->GetHiggsTTMomentum(true).M(),
                    eventInfo->GetHiggsBB().GetMomentum().M())) continue;
                if(eventInfo->GetLeg(1)->charge() == eventInfo->GetLeg(2)->charge()) continue;

                std::vector<TriggerDescriptorCollection::BitsContainer> reco_jet_matches;
                if(eventInfo->HasVBFjetPair()) {
                    reco_jet_matches.push_back(eventInfo->GetVBFJet(1)->triggerFilterMatch());
                    reco_jet_matches.push_back(eventInfo->GetVBFJet(2)->triggerFilterMatch());

                }

                for(size_t desc_id = 0; desc_id < trigger_descs.size(); ++desc_id) {
                    if(PassIsolation(channel, *eventInfo, trigger_descs.at(desc_id).tauId,
                                     trigger_descs.at(desc_id).tauIdWP)
                    && eventInfo->GetTriggerResults().AnyAcceptAndMatchEx(trigger_descs.at(desc_id).trigger,
                                                                          eventInfo->GetFirstLeg().GetMomentum().pt(),
                                                                          eventInfo->GetSecondLeg().GetMomentum().pt(),
                                                                          reco_jet_matches))
                        ++passed.at(desc_id);
                }
            }

            for(size_t desc_id = 0; desc_id < trigger_descs.size(); ++desc_id) {
                numerators.at(desc_id).SetBinContent(static_cast<int>(fileID+1), passed.at(desc_id));
                numerators.at(desc_id).SetBinError(static_cast<int>(fileID+1), std::sqrt(passed.at(desc_id)));
            }
        }

        std::vector<std::shared_ptr<TGraphAsymmErrors>> graphs;
        for(size_t desc_id = 0; desc_id < trigger_descs.size(); ++desc_id){
            TEfficiency eff = CreateEfficiency(numerators.at(desc_id), deno);
            graphs.emplace_back(CreateTGraph(eff, mass_list));
        }

        std::vector<std::shared_ptr<TGraph>> graphsRatio;
        for(size_t desc_id = 1; desc_id < trigger_descs.size(); ++desc_id)
            graphsRatio.emplace_back(CreateRatio(*graphs.at(desc_id), *graphs.at(0)));

        PrintTGraph(graphs, graphsRatio);
        canvas.Print((args.output_file() + ".pdf]").c_str());
    }

private:
    SampleDesc LoadConfig(std::vector<TriggerDesc>& trigger_descs) const
    {
        PropertyList properties;
        PropertyConfigReader reader;
        reader.Parse(args.file_cfg_name());
        const auto& items = reader.GetItems();

        SampleDesc sample_desc;
        const auto& sample_params = items.at(args.sample_type());
        sample_desc.file_name = sample_params.Get<>("file_name");
        sample_desc.points_list = sample_params.properties.GetList<double>("points", false);
        sample_desc.x_title = sample_params.Get<>("x_title");
        sample_desc.massWindowParams = sample_params.Get<EllipseParameters>("massWindowParams");

        for(const auto& item : items) {
            const auto& item_name = item.first;
            const auto& params = item.second;

            if(item_name.substr(0,4) != "trig") continue;

            Channel channel = params.Get<Channel>("channel");
            if(channel != args.channel()) continue;

            auto trigger_path = params.properties.GetList<>("trigger_paths", false);
            auto title_name = params.Get<>("title");
            auto color = params.Get<root_ext::Color>("color");
            TauIdDiscriminator tauId = params.Get<TauIdDiscriminator>("tauId");
            DiscriminatorWP tauIdWP = params.Get<DiscriminatorWP>("tauIdWP");

            trigger_descs.emplace_back(TriggerDesc{trigger_path, title_name, color, tauId, tauIdWP });
        }
        return sample_desc;
    }

    bool PassIsolation(Channel channel, EventInfoBase& event_base, TauIdDiscriminator discr,
                       DiscriminatorWP wp) const
    {
        if(channel == Channel::ETau) {
            const auto& tau = event_base.GetSecondLeg();
            return tau->Passed(discr, wp);
        }
        else if(channel == Channel::MuTau) {
            const auto& tau = event_base.GetSecondLeg();
            return tau->Passed(discr, wp);
        }
        else if(channel == Channel::TauTau) {
            const auto& tau1 = event_base.GetFirstLeg();
            const auto& tau2 = event_base.GetSecondLeg();
            return tau1->Passed(discr, wp) && tau2->Passed(discr, wp);
        }
        else
            throw exception("channel not supported");
    }

 private:
    Arguments args;
    TCanvas canvas;
    std::vector<TriggerDesc> trigger_descs;
    SampleDesc sample_desc;
    SignalObjectSelector signalObjectSelector;
};
}
PROGRAM_MAIN(analysis::EfficiencyStudy, Arguments)
