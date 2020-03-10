/*! Study of the efficiency of the cuts made to the events at baseline selection.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "TEfficiency.h"
#include <TCanvas.h>
#include "AnalysisTools/Print/include/PdfPrinter.h"
#include <boost/algorithm/string.hpp>
#include "AnalysisTools/Core/include/TextIO.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"

struct Arguments {
        run::Argument<std::vector<std::string>> input_file{"input_file", "Input file of the samples"};
        run::Argument<std::string> channels{"channels", "channels on which we work", "eTau,muTau,tauTau"};
        run::Argument<std::vector<std::string>> Cuts{"Cuts", "Vetos applied"};
     };

namespace analysis {

class EffiencyStudy {
 public:

     using TEffPtr = std::shared_ptr<TEfficiency>;
     using EffDesc = std::map<std::string,analysis::StVariable>;
     using FileDescCollection = std::map<std::string, EffDesc>;

     struct ChannelDesc {
         std::string channel;
         std::vector<std::string> histograms;
         std::vector<std::string> binNames;

         ChannelDesc() {}
         ChannelDesc(const std::string& _channel, const std::vector<std::string>& _histograms,
             const std::vector<std::string>& _binNames)
            : channel(_channel), histograms(_histograms), binNames(_binNames) {}
     };

     EffiencyStudy(const Arguments& _args) :
         args(_args), canvas("","", 600, 600)
     {
         auto channel_list = SplitValueList(args.channels(), false, ", ");
         active_channels.insert(channel_list.begin(), channel_list.end());
     }

    TEffPtr CreateEfficiency(const TH1& passed, const TH1& total, const std::string& prefix, const std::string& channel,
                          const std::string& hist_name)
    {
        constexpr static double oneSigma = 0.682689492137;

        TH1D empty_hist("", "", passed.GetNbinsX(), passed.GetBinLowEdge(1),
                        passed.GetBinLowEdge(passed.GetNbinsX()+1));

        if(!TEfficiency::CheckConsistency(passed, total))
            throw exception("passed TEfficiency objects do not have consistent bin contents");

        auto eff = std::make_shared<TEfficiency>(passed, total);

        eff->SetConfidenceLevel(oneSigma);
        eff->SetStatisticOption(TEfficiency::kFCP);

        std::ostringstream ss_name;
        ss_name << prefix << "_" << channel << "_" << hist_name;
        std::string name = ss_name.str();

        eff->SetTitle(name.c_str());

    return eff;
    }

    std::vector<TEffPtr> CreateEfficiencies(const TH1& selection_hist, const std::string& channel,
                            const std::string& hist_name, int shift, bool create_relative,
                            const std::string& bin_name, int& bin_id)
    {
        int N = selection_hist.GetNbinsX()-shift;
        TH1D deno("", "", N, 0, N);
        TH1D nume(deno), total(deno);

        for(int n = 1; n <= N; ++n){
            double totalWithWeight = std::ceil(selection_hist.GetBinContent(1)*1);

            deno.SetBinContent(n, selection_hist.GetBinContent(n + (shift -1)));
            deno.SetBinError(n, selection_hist.GetBinError(n + (shift -1)));
            nume.SetBinContent(n, selection_hist.GetBinContent(n+shift));
            nume.SetBinError(n, selection_hist.GetBinError(n+shift));
            total.SetBinContent(n, totalWithWeight);
            total.SetBinError(n, std::sqrt(totalWithWeight));
            const std::string bin_label = selection_hist.GetXaxis()->GetBinLabel(n+shift);
            if(bin_label == bin_name)
                bin_id = n;
        }

        std::vector<TEffPtr> effHistos;

        auto eff_abs = CreateEfficiency(nume, total, "eff_abs", channel, hist_name);
        effHistos.push_back(eff_abs);
        if(create_relative) {
            auto eff_rel = CreateEfficiency(nume, deno, "eff_rel", channel, hist_name);
            effHistos.push_back(eff_rel);
        }
        return effHistos;
   }

     void Run()
     {
         FileDescCollection fileInformation;

         std::vector<ChannelDesc> Channels = {
             ChannelDesc("eTau", {"events", "SignalTaus_Central"} , {"taus", "againstElectron"})
        };

         for(size_t cutId = 0; cutId < args.Cuts().size(); cutId++){

             std::string Cut = args.Cuts().at(cutId);
             auto file = root_ext::OpenRootFile(args.input_file().at(cutId));

             for(const auto& desc : Channels) {
                 if(!active_channels.count(desc.channel)) continue;

                 for(size_t i = 0; i < desc.histograms.size(); i++){
                    const std::string full_hist_name = desc.channel + "_stat/Selection_" + desc.histograms.at(i);
                    auto selection_hist = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*file, full_hist_name));

                    int bin_id = -1;
                    auto eff_vec = CreateEfficiencies(*selection_hist, desc.channel, desc.histograms.at(i), 1, true,
                                                      desc.binNames.at(i), bin_id);
                    if(bin_id == -1)
                        throw exception("Bin '%1%' not found.") % desc.binNames.at(i);
                    for(auto eff : eff_vec) {

                             double value = eff->GetEfficiency(bin_id);
                              double error_up = eff->GetEfficiencyErrorUp(bin_id);
                              double error_low = eff->GetEfficiencyErrorLow(bin_id);
                        fileInformation[eff->GetTitle()][Cut] = analysis::StVariable(value, error_up, error_low);

                    }
                }
            }
        }
        PrintValues(fileInformation);
    }

    private:

    void PrintValues(const FileDescCollection& fileInformation) const
    {
        std::cout << " \\begin{table}\n";
        std::cout << "      \\begin{center}\n";
        std::cout << "          \\begin{tabular}{|c|c|c|c|c|c|}   \\hline\n";

        for(size_t n = 0; n < args.Cuts().size(); n++){
          std::cout << "    &"<< args.Cuts().at(n);
        }
         std::cout << " \\\\ \\hline \\hline\n";
         std::cout << "    \\textbf{SM} & &  &  & &  \\\\ \n";
         for(const auto& item: fileInformation) {
             const std::string& hist_name = item.first;
             std::string result = boost::replace_all_copy(hist_name, "_", "\\_");
             std::cout << "    " << result;
             for(size_t cutId = 0; cutId < args.Cuts().size(); cutId++){
                 const std::string& Cut = args.Cuts().at(cutId);
                 const auto& eff = item.second.at(Cut);

                 const analysis::StVariable stValue(eff.value*100, eff.error_up*100, eff.error_low*100);
                 std::string GetEfficiency = stValue.ToLatexString();
                 std::cout << " &   $ " <<GetEfficiency << "\\% $" ;
            }
             std::cout << " \\\\"<< '\n';
         }
         std::cout << "    \\hline \\hline \n";
         std::cout << "          \\end{tabular}\n" ;
         std::cout << "          \\caption{ Channel: eTau} \n" ;
         std::cout << "      \\end{center}\n" ;
         std::cout << "\\end{table}\n" ;
     }

    private:
        Arguments args;
        TCanvas canvas;
        std::set<std::string> active_channels;
};
}
 PROGRAM_MAIN(analysis::EffiencyStudy, Arguments)
