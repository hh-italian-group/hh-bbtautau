/*! Merge DYJets files.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "Instruments/include/TTFileConfigEntryReader.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "Instruments/include/SampleDescriptor.h"


struct Arguments {
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> input_path{"input_path", "Input path of the samples"};
    run::Argument<std::string> cfg_name{"cfg_name", "cfg bin splitting"};
    run::Argument<std::string> file_cfg_name{"file_cfg_name", "TT file cfg"};
    run::Argument<std::string> output_file{"output_file", "Output file"};
};

namespace analysis {

namespace sample_merging{

class TTFileMerger {
public:
    using GenEventTypeMap = std::map<GenEventType, double>;
    using VectorSampleDescriptor = std::vector<SampleDescriptor<TTBinDescriptor, GenEventTypeMap>>;
    using VectorDYBinDescriptor = std::vector<TTBinDescriptor>;
    TTFileMerger(const Arguments& _args) : args(_args)
    {
        LoadInputs();
        output_bins = TTBinDescriptor::LoadConfig(args.cfg_name());
    }

public:

    void Run()
    {
        for(auto& output_bin : output_bins)
        {
            CalculateWeight(output_bin);
        }
        TTBinDescriptor::SaveCfg(args.output_file(), output_bins);
    }


private:
    SampleDescriptor<TTBinDescriptor, GenEventTypeMap> global_map;
    SampleDescriptor<TTBinDescriptor, GenEventTypeMap> inclusive;
    VectorSampleDescriptor all_samples;
    VectorDYBinDescriptor output_bins;
    Arguments args;


    static const std::set<std::string>& GetEnabledBranches()
    {
        static const std::set<std::string> EnabledBranches_read = {
            "eventEnergyScale", "genEventType", "genEventWeight"
        };
        return EnabledBranches_read;
    }

    
    void LoadInputs()
    {
        analysis::ConfigReader config_reader;

        TTBinDescriptorCollection file_descriptors;
        TTFileConfigEntryReader file_entry_reader(file_descriptors);
        config_reader.AddEntryReader("FILE", file_entry_reader, true);

        config_reader.ReadConfig(args.file_cfg_name());

        for (auto file_descriptor : file_descriptors){ //loop on DYJets files
            const TTBinDescriptor file_descriptor_element = file_descriptor.second;

            SampleDescriptor<TTBinDescriptor, GenEventTypeMap> sample_desc;
            sample_desc.bin = file_descriptor_element;
            //global_map.bin = file_descriptor_element;
            if (file_descriptor_element.fileType == FileType::inclusive)
                inclusive.bin = file_descriptor_element;

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


    void CalculateWeight(TTBinDescriptor& output_bin) const
    {
        double all_events = global_map.Integral(output_bin);
        double sample_contribution = inclusive.Integral(output_bin);
        // formula 2
        
        PhysicalValue nu_incl(sample_contribution,
                              sqrt(sample_contribution));
        PhysicalValue weight (nu_incl.GetValue()/all_events,
                              (all_events - sample_contribution)/std::pow(all_events,2)*
                              sqrt(sample_contribution));
        output_bin.nu = nu_incl;
        output_bin.weight = weight;
        output_bin.inclusive_integral = inclusive.Integral();

        if(output_bin.nu.GetStatisticalError() == std::numeric_limits<double>::infinity())
            throw exception("ref not found");
    }



};

} //namespace sample_merging

} //namespace analysis

PROGRAM_MAIN(analysis::sample_merging::TTFileMerger, Arguments)
