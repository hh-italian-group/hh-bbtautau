/*! Merge DYJets files.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "Instruments/include/DYFileConfigEntryReader.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "Instruments/include/SampleDescriptor.h"


struct Arguments {
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> cfg_name{"cfg_name", "cfg bin splitting"};
    run::Argument<std::string> file_cfg_name{"file_cfg_name", "DY file cfg"};
    run::Argument<std::string> output_file{"output_file", "Output file"};
};

namespace analysis {

namespace sample_merging{

class DYFileMerger {
public:
    using GenMap = ntuple::GenEventCountMap;
    using VectorSampleDescriptor = std::vector<SampleDescriptor<DYBinDescriptor, ntuple::GenEventCountMap>>;
    using VectorDYBinDescriptor = std::vector<DYBinDescriptor>;
    DYFileMerger(const Arguments& _args) : args(_args)
    {
        LoadInputs();
        output_bins = DYBinDescriptor::LoadConfig(args.cfg_name());
    }

public:

    void Run()
    {
        for(auto& output_bin : output_bins)
        {
            CalculateWeight(output_bin);
        }
        DYBinDescriptor::SaveCfg(args.output_file(), output_bins);
    }


private:
    SampleDescriptor<DYBinDescriptor, ntuple::GenEventCountMap> global_map;
    SampleDescriptor<DYBinDescriptor, ntuple::GenEventCountMap> inclusive;
    VectorSampleDescriptor all_samples;
    VectorDYBinDescriptor output_bins;
    Arguments args;


    void LoadInputs()
    {
        analysis::ConfigReader config_reader;

        DYBinDescriptorCollection file_descriptors;
        DYFileConfigEntryReader file_entry_reader(file_descriptors);
        config_reader.AddEntryReader("FILE", file_entry_reader, true);

        config_reader.ReadConfig(args.file_cfg_name());

        for (auto file_descriptor : file_descriptors){ //loop on DYJets files
            const DYBinDescriptor file_descriptor_element = file_descriptor.second;
            std::cout << "File descriptor characteristics: " << file_descriptor.first << ", " <<
                         file_descriptor_element.file_path << ", " << file_descriptor_element.fileType
                      << std::endl;
            SampleDescriptor<DYBinDescriptor, ntuple::GenEventCountMap> sample_desc;
            sample_desc.bin = file_descriptor_element;
            //global_map.bin = file_descriptor_element;
            if (file_descriptor_element.fileType == FileType::inclusive)
                inclusive.bin = file_descriptor_element;

            auto inputFile = root_ext::OpenRootFile(file_descriptor_element.file_path);
            ntuple::SummaryTuple summaryTuple(args.tree_name(), inputFile.get(), true);
            const Long64_t n_entries = summaryTuple.GetEntries();

            for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on entries
                summaryTuple.GetEntry(current_entry);
                for (size_t i = 0; i < summaryTuple.data().lhe_n_partons.size(); ++i){ // loop on gen event info
                    const ntuple::GenId genId(summaryTuple.data().lhe_n_partons.at(i),
                                              summaryTuple.data().lhe_n_b_partons.at(i),
                                              summaryTuple.data().lhe_ht10_bin.at(i));
                    size_t nevents = summaryTuple.data().lhe_n_events.at(i);
                    sample_desc.gen_counts[genId] += nevents;
                    global_map.gen_counts[genId] += nevents;
                    if (file_descriptor_element.fileType == FileType::inclusive)
                        inclusive.gen_counts[genId] += nevents;

                } // end loop on gen event info
            } //end loop on entries
            all_samples.push_back(sample_desc);
        } //end loop n file_descriptors
    }


    void CalculateWeight(DYBinDescriptor& output_bin) const
    {
        size_t all_events = global_map.Integral(output_bin);
        for(auto& sample : all_samples) {
            size_t contribution = sample.Integral(output_bin);
            if(!contribution) continue;
            //formula 2
            PhysicalValue nu ( contribution / double(sample.Integral()), sqrt(contribution)/double(sample.Integral()));
            PhysicalValue weight (nu.GetValue()/double(all_events), (double(all_events - contribution)/std::pow(all_events,2))*sqrt(contribution)/double(sample.Integral()));

            if(!(sample.bin.fileType == FileType::inclusive)) {
                size_t sample_contribution = inclusive.Integral(sample.bin);
                // formula 3
                PhysicalValue nu_incl(sample_contribution/double(inclusive.Integral()),
                                      sqrt(sample_contribution)/double(inclusive.Integral()));
                nu *= nu_incl;
                weight *= nu_incl;
            }

            if(output_bin.nu.GetStatisticalError() > nu.GetStatisticalError()) {
                output_bin.nu = nu;
                output_bin.ref_sample = sample.bin.name;
                output_bin.weight = weight;
                output_bin.fileType = sample.bin.fileType;
            }

        }

        if(output_bin.nu.GetStatisticalError() == std::numeric_limits<double>::infinity())
            throw exception("ref not found");
    }



};

} //namespace sample_merging

} //namespace analysis

PROGRAM_MAIN(analysis::sample_merging::DYFileMerger, Arguments)
