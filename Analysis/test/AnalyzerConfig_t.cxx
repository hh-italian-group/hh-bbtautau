/*! Check Config file.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "Analysis/include/AnalyzerConfigEntryReader.h"
#include "Analysis/include/SampleDescriptorConfigEntryReader.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"


struct Arguments {
    run::Argument<std::string> file_cfg_name{"file_cfg_name", "Source file cfg"};
};

namespace analysis {

class AnalyzerConfig_t {
public:

    AnalyzerConfig_t(const Arguments& _args) :
        args(_args)
    {}

    void Run()
    {
        analysis::ConfigReader config_reader;

        AnalyzerSetupCollection ana_setup;
        AnalyzerConfigEntryReader ana_entry_reader(ana_setup);
        config_reader.AddEntryReader("ANA_DESC", ana_entry_reader, false);

        SampleDescriptorCollection sample_descriptors;
        SampleDescriptorConfigEntryReader sample_entry_reader(sample_descriptors);
        config_reader.AddEntryReader("SAMPLE", sample_entry_reader, true);

        CombineSampleDescriptorCollection combine_descriptors;
        CombineSampleDescriptorConfigEntryReader combine_entry_reader(combine_descriptors,sample_descriptors);
        config_reader.AddEntryReader("SAMPLE_CMB", combine_entry_reader, false);

        config_reader.ReadConfig(args.file_cfg_name());

        for (auto ana_descriptor : ana_setup){ //loop ana_descriptors
            const AnalyzerSetup ana_descriptor_element = ana_descriptor.second;
            const std::string& name = ana_descriptor_element.name;

            std::cout << "Analyzer descriptor characteristics: " << name << ", " <<
                            ana_descriptor_element.int_lumi << ", " <<
                         ana_descriptor_element.apply_mass_cut << ", " <<
                         ana_descriptor_element.energy_scales.size() << ", " <<
                         ana_descriptor_element.final_variables.size() << ", " <<
                         ana_descriptor_element.int_lumi << ", " << std::endl;

        } //end loop ana_descriptors

        for (auto sample_descriptor : sample_descriptors){ //loop sample_descriptors
            const SampleDescriptor sample_descriptor_element = sample_descriptor.second;
            const std::string& name = sample_descriptor_element.name;

            std::cout << "Sample descriptor characteristics: " << name << ", " <<
                            sample_descriptor_element.categoryType << ", " <<
                         sample_descriptor_element.channels.size() << ", " <<
                         sample_descriptor_element.color << ", " <<
                         sample_descriptor_element.file_path_pattern << ", " <<
                         sample_descriptor_element.cross_section << ", " <<
                         sample_descriptor_element.signal_points_raw.size() << ", " <<
                         sample_descriptor_element.GetNSignalPoints() << ", " << std::endl;

            if(!sample_descriptor_element.signal_points_raw.size()) continue;
            else
                std::cout<< "FileName: " << sample_descriptor_element.GetFileName(sample_descriptor_element.GetNSignalPoints() - 1) << std::endl; //

        } //end loop sample_descriptors

        for (auto combine_descriptor : combine_descriptors){ //loop combine_descriptors
            const CombineSampleDescriptor combine_descriptor_element = combine_descriptor.second;
            const std::string& name = combine_descriptor_element.name;

            std::cout << "Combine descriptor characteristics: " << name << ", " <<
                            combine_descriptor_element.sample_descriptors.size() << std::endl;
        } //end loop combine_descriptors

    }

private:
    Arguments args;

};

} //namespace analysis

PROGRAM_MAIN(analysis::AnalyzerConfig_t, Arguments)
