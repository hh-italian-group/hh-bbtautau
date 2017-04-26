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
    TH1D_ENTRY(lhe_hh_m, 100, 0, 1000)
    TH1D_ENTRY(lhe_hh_cosTheta, 6, -1.0, 1.0)
    //TH2D_ENTRY(lhe_hh_cosTheta_vs_m, 25, 200, 2000, 30, -1.1, 1.1)
    TH2D_ENTRY_CUSTOM(lhe_hh_cosTheta_vs_m, mhh_Bins(), cosTheta_Bins())
    TH2D_ENTRY_CUSTOM(weight, mhh_Bins(), cosTheta_Bins())

    virtual const std::vector<double>& cosTheta_Bins() const
    {
        static const std::vector<double> bins = { -1.0, -0.6, -0.2, 0.2, 0.6, 1.0 };
        return bins;
    }

    virtual const std::vector<double>& mhh_Bins() const
    {
        static const std::vector<double> bins = { 250.,270.,300.,330.,360.,390., 420.,450.,500.,550.,600.,700.,800.,1000. };
        return bins;
    }
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
                    anaData.lhe_hh_cosTheta(name).Fill(event->lhe_hh_cosTheta);
                    anaData.lhe_hh_cosTheta_vs_m(name).Fill(event->lhe_hh_m,event->lhe_hh_cosTheta);
                } //end loop on entries

            } //end loop on files
            const double scale = 1 / anaData.lhe_hh_cosTheta_vs_m(name).Integral(0,anaData.lhe_hh_cosTheta_vs_m(name).GetNbinsX()+1,
                                                                                 0, anaData.lhe_hh_cosTheta_vs_m(name).GetNbinsY()+1);
            anaData.lhe_hh_cosTheta_vs_m(name).Scale(scale);
            if (file_descriptor_element.fileType == FileType::sm)
                name_sm = name;
        } //end loop n file_descriptors

        for (const auto& file_descriptor : file_descriptors) {
            const std::string& name = file_descriptor.second.name;

            for (unsigned n = 0; n <= anaData.lhe_hh_cosTheta_vs_m(name).GetNbinsX()+1; ++n){
                for (unsigned h = 0; h <= anaData.lhe_hh_cosTheta_vs_m(name).GetNbinsY()+1; ++h){
                    if (anaData.lhe_hh_cosTheta_vs_m(name).GetBinContent(n,h) == 0)
                        std::cout << name <<  " - Empty Bin! (" << n << "," << h << ")" << std::endl ;
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
