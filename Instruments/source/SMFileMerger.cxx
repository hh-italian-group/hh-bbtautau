/*! Merge BSM files.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "Instruments/include/SMFileConfigEntryReader.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "Instruments/include/SampleDescriptor.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"


struct Arguments {
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> input_path{"input_path", "Input path of the samples"};
    run::Argument<std::string> file_cfg_name{"file_cfg_name", "SM file cfg"};
    run::Argument<std::string> output_file{"output_file", "Output root file"};
};

namespace analysis {

namespace sample_merging{

class SMFileMergerData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;


//    const std::vector<double> mhh_Bins{ 245, 270, 300, 330, 360, 390, 420, 450, 500, 550, 600, 700, 800, 1000, 10000 };
//    const std::vector<double> cosTheta_Bins{ 0.0, 0.4, 0.6, 0.8, 0.9, 1.0 };

    const std::vector<double> mhh_Bins{ 250,260,270,280,290,300,310,320,330,340,
                                        350,360,370,380,390,400,410,420,430,440,
                                        450,460,470,480,490,
                                        500,510,520,530,540,550,600,610,620,630,
                                        640,650,660,670,680,690,700,750,800,850,
                                        900,950,1000,1100,1200,1300,1400,1500.,1750,2000,50000};
    const std::vector<double> cosTheta_Bins{0.0, 0.4, 0.6, 0.8, 1.0} ;

    TH1D_ENTRY_CUSTOM(lhe_hh_m, mhh_Bins)
    TH1D_ENTRY_CUSTOM(lhe_hh_cosTheta, cosTheta_Bins)
    TH2D_ENTRY_CUSTOM(lhe_hh_cosTheta_vs_m, mhh_Bins, cosTheta_Bins)
    TH2D_ENTRY_CUSTOM(weight, mhh_Bins, cosTheta_Bins)
};


class SMFileMerger {
public:
    SMFileMerger(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.output_file())), anaData(output)
    {
        LoadInputs();
    }

    void Run() {}

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    SMFileMergerData anaData;

    void LoadInputs()
    {
        analysis::ConfigReader config_reader;

        SMBinDescriptorCollection file_descriptors;
        SMFileConfigEntryReader file_entry_reader(file_descriptors);
        config_reader.AddEntryReader("FILE", file_entry_reader, true);

        config_reader.ReadConfig(args.file_cfg_name());
        std::string name_sm;

        for (const auto& file_descriptor : file_descriptors) {
            const SMBinDescriptor file_descriptor_element = file_descriptor.second;
            const std::string& name = file_descriptor_element.name;

            for (auto single_file_path : file_descriptor_element.file_paths){ //loop on files
                std::cout << "File descriptor characteristics: " << file_descriptor.first << ", " <<
                             single_file_path << ", " << file_descriptor_element.fileType
                          << std::endl;
                const std::string filename = args.input_path()  + "/" + single_file_path;
                auto inputFile = root_ext::OpenRootFile(filename);


                std::shared_ptr<ntuple::ExpressTuple> AllEventTuple(new ntuple::ExpressTuple(args.tree_name(), inputFile.get(), true));
                const Long64_t n_entries = AllEventTuple->GetEntries();

                for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on entries

                    AllEventTuple->GetEntry(current_entry);
                    std::shared_ptr<ntuple::ExpressEvent> event(new ntuple::ExpressEvent(AllEventTuple->data()));

                    //std::cout << "Mhh: " << event->lhe_hh_m << ", cosTheta: " << event->lhe_hh_cosTheta << std::endl;
                    anaData.lhe_hh_m(name).Fill(event->lhe_hh_m);
                    anaData.lhe_hh_cosTheta(name).Fill(std::abs(event->lhe_hh_cosTheta));
                    anaData.lhe_hh_cosTheta_vs_m(name).Fill(event->lhe_hh_m,std::abs(event->lhe_hh_cosTheta));
                    anaData.lhe_hh_cosTheta_vs_m().Fill(event->lhe_hh_m,std::abs(event->lhe_hh_cosTheta));
                } //end loop on entries

            } //end loop on files
//            const double scale = 1 / anaData.lhe_hh_cosTheta_vs_m(name).Integral(0,anaData.lhe_hh_cosTheta_vs_m(name).GetNbinsX()+1,
//                                                                                 0, anaData.lhe_hh_cosTheta_vs_m(name).GetNbinsY()+1);
//            anaData.lhe_hh_cosTheta_vs_m(name).Scale(scale);
            RenormalizeHistogram(anaData.lhe_hh_cosTheta_vs_m(name), 1, true);
            if (file_descriptor_element.fileType == FileType::sm)
                name_sm = name;
        } //end loop n file_descriptors
        RenormalizeHistogram(anaData.lhe_hh_cosTheta_vs_m(), 1, true);


        const std::vector<std::string> names = { "", name_sm };
        for(const auto& name : names){
            const auto& hist = anaData.lhe_hh_cosTheta_vs_m(name);
            const Int_t N = hist.GetNbinsX() + 1;
            const Int_t H = hist.GetNbinsY() + 1;
            for (Int_t n = 1; n <= N; ++n){
                for (Int_t h = 0; h <= H; ++h){
                    const bool zero_content = anaData.lhe_hh_cosTheta_vs_m(name).GetBinContent(n,h) == 0;
                    const bool overflow_bin = n == 0 || n == N || h == 0 || h == H;
                    if((zero_content && !overflow_bin) || (!zero_content && overflow_bin)) {
                        std::ostringstream ss;
                        const std::string prefix = overflow_bin ? "Non empty" : "Empty";
                        ss << prefix << " bin in " << name << ": (" << n << ", " << h << ") with center at ("
                           << hist.GetXaxis()->GetBinCenter(n) << ", " << hist.GetYaxis()->GetBinCenter(h) << ").";
                        throw exception(ss.str());
                    }
                }
            }
            anaData.weight(name).CopyContent(anaData.lhe_hh_cosTheta_vs_m(name_sm));
            anaData.weight(name).Divide(&anaData.lhe_hh_cosTheta_vs_m(name));
        }

    }
};

} //namespace sample_merging

} //namespace analysis

PROGRAM_MAIN(analysis::sample_merging::SMFileMerger, Arguments)
