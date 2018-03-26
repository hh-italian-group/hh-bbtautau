#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "Instruments/include/SampleDescriptor.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/PropertyConfigReader.h"
#include "TEfficiency.h"
#include <TCanvas.h>
#include <exception.h>
#include "AnalysisTools/Print/include/PdfPrinter.h"

struct Arguments {
        run::Argument<std::string> channel_name{"channel_name", "channel on which we work"};
        run::Argument<std::string> input_file{"input_file", "Input file of the samples"};
        run::Argument<std::string> output_file{"output_file", "Output pdf file"};
   };


 namespace analysis {


 class EffiencyStudy {
 public:

     struct ChannelDesc {
         std::string channel;
         std::vector<std::string> histograms;
         double branching_ratio;

         ChannelDesc() {}
         ChannelDesc(const std::string& _channel, const std::vector<std::string>& _histograms, double br)
            : channel(_channel), histograms(_histograms), branching_ratio(br) {}
     };

     EffiencyStudy(const Arguments& _args) :
         args(_args), canvas("","", 600, 600)
     {
         gStyle->SetOptStat(0);
         canvas.SetGrid();
         canvas.Print((args.output_file() + ".pdf[").c_str());
         canvas.Draw();
     }

    void CreateEfficiency(const TH1& passed, const TH1& total, const std::string& prefix, const std::string& channel,
                          const std::string& hist_name)
    {
        if(!TEfficiency::CheckConsistency(passed, total))
            throw analysis::exception("passed TEfficiency objects do not have consistent bin contents");

        TEfficiency eff(passed, total);

        eff.SetConfidenceLevel(0.682689492137);
        eff.SetStatisticOption(TEfficiency::kFCP);

        std::ostringstream ss_name;
        ss_name << prefix << "_" << channel << "_" << hist_name;
        std::string name = ss_name.str();

        eff.SetTitle(name.c_str());
        canvas.Clear();
        TH1D empty_hist("", "", passed.GetNbinsX(), passed.GetBinLowEdge(1),
                        passed.GetBinLowEdge(passed.GetNbinsX()+1));
        double min = 1, max = 0;
        for(int n = 1; n <= passed.GetNbinsX(); ++n) {
            min = std::min(eff.GetEfficiency(n), min);
            max = std::max(eff.GetEfficiency(n), max);
            empty_hist.GetXaxis()->SetBinLabel(n, passed.GetXaxis()->GetBinLabel(n));
            empty_hist.SetLabelSize(0.04f);
        }

        static const double range_sf = 1.1;
        if(min <= 0.2){
            empty_hist.GetYaxis()->SetRangeUser(0, max * range_sf);
        }
        else{
            empty_hist.GetYaxis()->SetRangeUser(min / range_sf, max * range_sf);
        }
        empty_hist.Draw();
        eff.Draw("SAME P");
        LaTeXDraw(eff, empty_hist);
        empty_hist.SetTitle(name.c_str());


        canvas.Print((args.output_file()+".pdf").c_str(), ("Title:"+name).c_str());
        canvas.Clear();
    }

    void CreateEfficiencies(const TH1& selection_hist, const std::string& channel,
                            const std::string& hist_name, int shift, bool create_relative, double sf)
    {
        int N = selection_hist.GetNbinsX()-shift;
        TH1D deno("", "", N, 0, N);
        TH1D nume(deno), total(deno), empty(deno);

        for(int n = 1; n <= N; ++n){
            double totalWithWeight = std::ceil(selection_hist.GetBinContent(1)*sf);

            deno.SetBinContent(n, selection_hist.GetBinContent(n));
            deno.SetBinError(n, selection_hist.GetBinError(n));
            nume.SetBinContent(n, selection_hist.GetBinContent(n+shift));
            nume.SetBinError(n, selection_hist.GetBinError(n+shift));
            total.SetBinContent(n, totalWithWeight);
            total.SetBinError(n, std::sqrt(totalWithWeight));
            nume.GetXaxis()->SetBinLabel(n, selection_hist.GetXaxis()->GetBinLabel(n+shift));
            empty.GetXaxis()->SetBinLabel(n, selection_hist.GetXaxis()->GetBinLabel(n+shift));

        }
        CreateEfficiency(nume, total, "eff_abs", channel, hist_name);
        if(create_relative)
        CreateEfficiency(nume, deno, "eff_rel", channel, hist_name);
       }

     void Run()
     {
         auto file = root_ext::OpenRootFile(args.input_file());

         std::vector<ChannelDesc> Channels = {
             ChannelDesc("eTau", {"events","SignalElectrons_Central", "SignalTaus_Central", "jets_Central" }, 0.23),
             ChannelDesc("muTau", {"events","SignalMuons_Central", "SignalTaus_Central", "jets_Central"}, 0.23),
             ChannelDesc("tauTau", {"events","SignalTaus_Central", "jets_Central"}, 0.54),
         };

         //Creating the efficiency histogrmas after events

        for(const auto& entry : Channels) {
            const std::string& channel_name = entry.channel;
            const std::vector<std::string> hist_name = entry.histograms;
            const double divide_by_br = entry.branching_ratio;

            for(unsigned long i = 1; i < hist_name.size(); i++){

                const std::string full_hist_name = channel_name + "_stat/Selection_" + hist_name.at(i);
                auto selection_hist = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*file, full_hist_name));
                    
                CreateEfficiencies(*selection_hist, channel_name, hist_name.at(i), 1, true, 1);

            }

            //Creating the efficiency histogrmas for events

            const std::string Events = channel_name + "_stat/Selection_" + hist_name.at(0);
            auto selection_events = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*file, Events));

            CreateEfficiencies(*selection_events, channel_name, hist_name.at(0), 1, true, 1);
            CreateEfficiencies(*selection_events, channel_name, hist_name.at(0) + "BR" , 2, false, divide_by_br);
        }


        canvas.Print((args.output_file() + ".pdf]").c_str());
     }
    private:

     void LaTeXDraw(TEfficiency eff, TH1D empty_histo){
         TLatex latex;
         latex.SetTextSize(0.024f);
         latex.SetTextAlign(13);

         for(int n = 1; n <= empty_histo.GetNbinsX(); ++n){

             std::ostringstream ss_GetEfficiency;
             double Efficiency= ((eff.GetEfficiency(n))*100);

              double ErrorUp= ((eff.GetEfficiencyErrorUp(n))*100);
              double ErrorLow= ((eff.GetEfficiencyErrorLow(n))*100);

            ss_GetEfficiency << std::setprecision(3) << Efficiency << "^{+"
                             << std::fixed << std::setprecision(1) << ErrorUp
                              << "}_{-" << std::setprecision(1) << ErrorLow << "}"
                              << "%" ;
            std::string GetEfficiency = ss_GetEfficiency.str();

            double max = eff.GetEfficiency(1);
            max = std::max(eff.GetEfficiency(n), max);

            latex.DrawLatex(n - 0.8, max*1.09, GetEfficiency.c_str() );
         }

     }
 private:

     Arguments args;
     std::shared_ptr<TFile> output;
     TCanvas canvas;
     };


 } //namespace analysis

 PROGRAM_MAIN(analysis::EffiencyStudy, Arguments)

