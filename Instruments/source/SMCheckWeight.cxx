/*! Check SM weight.
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
    run::Argument<std::string> weight_file{"weight_file", "file to reweight BSM samples"};
    run::Argument<std::string> SM_file{"SM_file", "SM file"};
    run::Argument<std::string> BSM_file{"BSM_file", "BSM file"};
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> output_file{"output_file", "Output root file"};
};

namespace analysis {

namespace sample_merging{

class SMCheckWeightData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    TH1D_ENTRY(lhe_hh_m, 100, 0, 1000)
    TH1D_ENTRY(lhe_hh_cosTheta, 300, -1.5, 1.5)
    //TH2D_ENTRY(lhe_hh_cosTheta_vs_m, 25, 200, 2000, 30, -1.1, 1.1)
    TH2D_ENTRY_CUSTOM(lhe_hh_cosTheta_vs_m, mhh_Bins(), cosTheta_Bins())
    TH2D_ENTRY_CUSTOM(weight, mhh_Bins(), cosTheta_Bins())
};


class SMCheckWeight {
public:
    SMCheckWeight(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.output_file())), anaData(output)
    {
        LoadInputs();
    }

    void Run() {}

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    SMCheckWeightData anaData;

    void LoadInputs()
    {
        auto inputFile_weight = root_ext::OpenRootFile(args.weight_file());
        auto inputFile_SM = root_ext::OpenRootFile(args.SM_file());
        auto inputFile_BSM = root_ext::OpenRootFile(args.BSM_file());

        TH2D *weight = (TH2D*)inputFile_weight->Get("weight_node_BSM");

        Int_t bin_x = weight->GetXaxis()->FindBin(350);
        Int_t bin_y = weight->GetYaxis()->FindBin(0.5);

        for (unsigned n = 0; n < weight->GetNbinsX(); ++n){
            for (unsigned h = 0; h < weight->GetNbinsY(); ++h){
                std::cout << " Bin Content: " << weight->GetBinContent(bin_x,bin_y) << std::endl ;
            }
        }

        analysis::ConfigReader config_reader;

        SMBinDescriptorCollection file_descriptors;
        SMFileConfigEntryReader file_entry_reader(file_descriptors);
        config_reader.AddEntryReader("FILE", file_entry_reader, true);

        config_reader.ReadConfig(args.file_cfg_name());
        std::string name_sm;

        for (auto file_descriptor : file_descriptors){ //loop on files
            const SMBinDescriptor file_descriptor_element = file_descriptor.second;
            const std::string& name = file_descriptor_element.name;

            const Channel channel = Parse<Channel>(args.tree_name());
            const Channel descriptor_channel = Parse<Channel>(file_descriptor_element.channel);
            if (descriptor_channel != channel) continue;

            //double count = 0;
            for (auto single_file_path : file_descriptor_element.file_paths){ //loop on files

                std::cout << "File descriptor characteristics: " << file_descriptor.first << ", " <<
                             single_file_path << ", " << file_descriptor_element.fileType
                          << std::endl;
                std::string filename = args.input_path()  + "/" + single_file_path;
                auto inputFile = root_ext::OpenRootFile(filename);
                ntuple::EventTuple eventTuple(args.tree_name(), inputFile.get(), true, {},
                                              GetEnabledBranches());
                const Long64_t n_entries = eventTuple.GetEntries();

                for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on entries
                    eventTuple.GetEntry(current_entry);
                    const ntuple::Event& event = eventTuple.data();
                    if (static_cast<EventEnergyScale>(event.eventEnergyScale) != analysis::EventEnergyScale::Central)
                        continue;
                    GenEventType genEventType = static_cast<GenEventType>(event.genEventType);

                        sample_desc.gen_counts[genEventType] += event.genEventWeight;
                        global_map.gen_counts[genEventType] += event.genEventWeight;
                        if (file_descriptor_element.fileType == FileType::inclusive)
                            inclusive.gen_counts[genEventType] += event.genEventWeight;


                } //end loop on entries
            } // end loop on files
            all_samples.push_back(sample_desc);
        } //end loop n file_descriptors

    }
};

} //namespace sample_merging

} //namespace analysis

PROGRAM_MAIN(analysis::sample_merging::SMCheckWeight, Arguments)
