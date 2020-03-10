/*! Study of the efficiency of the cuts made to the events at baseline selection.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include "AnalysisTools/Print/include/PdfPrinter.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"

struct Arguments {
        run::Argument<std::string> input_file{"input_file", "Input file of the samples"};
        run::Argument<std::string> output_file{"output_file", "Output pdf file"};
        run::Argument<std::string> channels{"channels", "channels on which we work", "eTau,muTau,tauTau"};
   };

namespace analysis {

class EffiencyStudy {
 public:

     struct ChannelDesc {
         std::string channel;
         std::vector<std::string> histograms;
         std::vector<bool> divideOrNot;
         double branching_ratio;


         ChannelDesc() {}
         ChannelDesc(const std::string& _channel, const std::vector<std::string>& _histograms, const std::vector<bool>& _divideOrNot,double br)
            : channel(_channel), histograms(_histograms), divideOrNot(_divideOrNot) ,branching_ratio(br){}
     };

     EffiencyStudy(const Arguments& _args) :
         args(_args), canvas("","", 600, 600)
     {
         gStyle->SetOptStat(0);
         canvas.SetGrid();
         canvas.Print((args.output_file() + ".pdf[").c_str());
         canvas.Draw();
         auto channel_list = SplitValueList(args.channels(), false, ", ");
         active_channels.insert(channel_list.begin(), channel_list.end());
     }

    void CreateEfficiency(const TH1& passed, const TH1& total, const std::string& prefix, const std::string& channel,
                          const std::string& hist_name)
    {
        constexpr static double range_sf = 1.1;
        constexpr static double oneSigma = 0.682689492137;

        TH1D empty_hist("", "", passed.GetNbinsX(), passed.GetBinLowEdge(1),
                        passed.GetBinLowEdge(passed.GetNbinsX()+1));

        if(!TEfficiency::CheckConsistency(passed, total))
            throw exception("passed TEfficiency objects do not have consistent bin contents");

        TEfficiency eff(passed, total);

        eff.SetConfidenceLevel(oneSigma);
        eff.SetStatisticOption(TEfficiency::kFCP);

        std::ostringstream ss_name;
        ss_name << prefix << "_" << channel << "_" << hist_name;
        std::string name = ss_name.str();

        eff.SetTitle(name.c_str());
        canvas.Clear();

        double min = 1, max = 0;
        for(int n = 1; n <= passed.GetNbinsX(); ++n) {
            min = std::min(eff.GetEfficiency(n), min);
            max = std::max(eff.GetEfficiency(n), max);
            empty_hist.GetXaxis()->SetBinLabel(n, passed.GetXaxis()->GetBinLabel(n));
            empty_hist.SetLabelSize(0.04f);
        }

        min = min > 0.2 ? min / range_sf : 0;
        empty_hist.GetYaxis()->SetRangeUser(min, max * range_sf);
        empty_hist.Draw();
        eff.Draw("SAME P");
        int nbins = passed.GetNbinsX();
        LaTeXDraw(eff, nbins, max);
        empty_hist.SetTitle(name.c_str());

        canvas.Print((args.output_file()+".pdf").c_str(), ("Title:"+name).c_str());
        canvas.Clear();
    }

    void CreateEfficiencies(const TH1& selection_hist, const std::string& channel,
                            const std::string& hist_name, int shift, bool create_relative, double sf)
    {
        int N = selection_hist.GetNbinsX()-shift;
        TH1D deno("", "", N, 0, N);
        TH1D nume(deno), total(deno)/*, empty(deno)*/;

        for(int n = 1; n <= N; ++n){
            double totalWithWeight = std::ceil(selection_hist.GetBinContent(1)*sf);

            deno.SetBinContent(n, selection_hist.GetBinContent(n + (shift -1)));
            deno.SetBinError(n, selection_hist.GetBinError(n + (shift -1)));
            nume.SetBinContent(n, selection_hist.GetBinContent(n+shift));
            nume.SetBinError(n, selection_hist.GetBinError(n+shift));
            total.SetBinContent(n, totalWithWeight);
            total.SetBinError(n, std::sqrt(totalWithWeight));
            nume.GetXaxis()->SetBinLabel(n, selection_hist.GetXaxis()->GetBinLabel(n+shift));

        }
        CreateEfficiency(nume, total, "eff_abs", channel, hist_name);
        if(create_relative)
            CreateEfficiency(nume, deno, "eff_rel", channel, hist_name);
    }

     void Run()
     {
         auto file = root_ext::OpenRootFile(args.input_file());

         std::vector<ChannelDesc> Channels = {
             ChannelDesc("eTau", {"events","SignalElectrons_Central", "SignalTaus_Central", "jets_Central"},
                         {true, false, false, false}, 0.231),
             ChannelDesc("muTau", {"events","SignalMuons_Central", "SignalTaus_Central", "jets_Central"},
                         {true, false, false, false}, 0.225),
             ChannelDesc("tauTau", {"events","SignalTaus_Central", "jets_Central"},
                         {true, false, false}, 0.420)
         };

         for(const auto& desc : Channels) {
             if(!active_channels.count(desc.channel)) continue;
            for(size_t i = 0; i < desc.histograms.size(); i++){
                const std::string full_hist_name = desc.channel + "_stat/Selection_" + desc.histograms.at(i);
                auto selection_hist = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*file, full_hist_name));
                CreateEfficiencies(*selection_hist, desc.channel, desc.histograms.at(i), 1, true, 1);
                if(desc.divideOrNot.at(i))
                    CreateEfficiencies(*selection_hist, desc.channel, desc.histograms.at(i) + "BR", 2, false,
                                       desc.branching_ratio);
            }
        }


        canvas.Print((args.output_file() + ".pdf]").c_str());
     }
    private:
        void LaTeXDraw(const TEfficiency& eff, int nBins, double maxSize) const
       {
           TLatex latex;
           latex.SetTextSize(0.020f);
           latex.SetTextAlign(23);

           for(int n = 1; n <= nBins; ++n){

               std::ostringstream ss_GetEfficiency;
               double Efficiency= ((eff.GetEfficiency(n))*100);
               double ErrorUp= ((eff.GetEfficiencyErrorUp(n))*100);
               double ErrorLow= ((eff.GetEfficiencyErrorLow(n))*100);

               const analysis::StVariable stValue(Efficiency, ErrorUp, ErrorLow);
               std::string GetEfficiency = stValue.ToLatexString();

               int decimals_to_print = stValue.decimals_to_print();


                if (decimals_to_print <= 3){
                   latex.SetTextAngle(0);
               }
               else if(decimals_to_print > 3){
                   latex.SetTextAngle(35);
               }
               latex.DrawLatex(n - 0.5, maxSize*1.09, GetEfficiency.c_str() );
           }
       }
 private:
     Arguments args;
     TCanvas canvas;
     std::set<std::string> active_channels;
};
}
 PROGRAM_MAIN(analysis::EffiencyStudy, Arguments)
