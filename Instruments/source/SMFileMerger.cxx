/*! Merge DYJets files.
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
    TH1D_ENTRY(lhe_hh_cosTheta, 300, -1.5, 1.5)
    TH2D_ENTRY(lhe_hh_cosTheta_vs_m, 100, 0, 1000, 300, -1.5, 1.5)
    TH2D_ENTRY(weight, 100, 0, 1000, 300, -1.5, 1.5)
};


class SMFileMerger {
public:

    SMFileMerger(const Arguments& _args) : args(_args), output(root_ext::CreateRootFile(args.output_file())),
                                                               anaData(output)
    {
        LoadInputs();
    }

public:

    void Run()
    {
    }


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
        std::vector<root_ext::SmartHistogram<TH2D>> histograms_sm;
        std::vector<root_ext::SmartHistogram<TH2D>> histograms_bsm;

        for (auto file_descriptor : file_descriptors){ //loop on DYJets files
            const SMBinDescriptor file_descriptor_element = file_descriptor.second;

                std::cout << "File descriptor characteristics: " << file_descriptor.first << ", " <<
                             file_descriptor_element.file_path << ", " << file_descriptor_element.fileType
                          << std::endl;
                std::string filename = args.input_path()  + "/" + file_descriptor_element.file_path;
                auto inputFile = root_ext::OpenRootFile(filename);
                std::string name = file_descriptor_element.name;

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

                double scale = 1/anaData.lhe_hh_cosTheta_vs_m(name).Integral();
                anaData.lhe_hh_cosTheta_vs_m(name).Scale(scale);
                if (file_descriptor_element.fileType == FileType::sm)
                    histograms_sm.push_back(anaData.lhe_hh_cosTheta_vs_m(name));
                else
                    histograms_bsm.push_back(anaData.lhe_hh_cosTheta_vs_m(name));


        } //end loop n file_descriptors

        std::cout << "Hist sm size: " << histograms_sm.size() << std::endl;
        std::cout << "Hist bsm size: " << histograms_bsm.size() << std::endl;

        for (unsigned n = 0; n < histograms_bsm.size(); ++n){
            root_ext::SmartHistogram<TH2D> hist_sm = histograms_sm.at(0);
            root_ext::SmartHistogram<TH2D> hist_bsm = histograms_bsm.at(n);

            std::shared_ptr<root_ext::SmartHistogram<TH2D>> ratio_histogram =
                    std::shared_ptr<root_ext::SmartHistogram<TH2D>>(new root_ext::SmartHistogram<TH2D>(hist_sm));
            std::shared_ptr<root_ext::SmartHistogram<TH2D>> histogram_bsm =
                    std::shared_ptr<root_ext::SmartHistogram<TH2D>>(new root_ext::SmartHistogram<TH2D>(hist_bsm));
            ratio_histogram->Divide(histogram_bsm.get());

            anaData.weight(n).Clone("ratio_histogram");
        }
    }

};

} //namespace sample_merging

} //namespace analysis

PROGRAM_MAIN(analysis::sample_merging::SMFileMerger, Arguments)
