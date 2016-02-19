/*! Limit configuration producer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Core/include/RootExt.h"
#include "Analysis/include/UncertaintyCalculatorCollection.h"

class LimitConfigurationProducer {
public:
    LimitConfigurationProducer(const std::string& uncConfigName, const std::string& srcConfigName,
                               const std::string& anaDataFileName,
                               const std::string& _outputPath)
        : anaDataReader(anaDataFileName), outputPath(_outputPath),
          dataCategories(srcConfigName, "", analysis::Channel::TauTau),
          calculators(uncertainties, dataCategories, anaDataReader,
                      analysis::EventSubCategory::KinematicFitConvergedWithMassWindow,
                      analysis::FlatAnalyzerData::m_ttbb_kinfit_Name())
    {
        using namespace analysis;
        using namespace analysis::limits;

        ConfigReader configReader(uncConfigName);
        SampleCategoryCollectionReader sampleReader(samples);
        CategoryDescriptorReader categoryReader(samples, categories);
        UncertaintyDescriptorReader uncertaintyReader(uncertainties);
        configReader.AddEntryReader("SAMPLES", sampleReader, false);
        configReader.AddEntryReader("CATEGORY", categoryReader, false);
        configReader.AddEntryReader("UNC", uncertaintyReader, true);

        configReader.ReadConfig();
    }

    void Run()
    {
        using namespace analysis::limits;

        for(const auto& categoryIter : categories) {
            const CategoryDescriptor& categoryDescriptor = categoryIter.second;
            std::cout << "Processing '" << categoryDescriptor.name << "' category.\n";

            std::cout << "Creating cgs configuration..." << std::endl;
            ProduceCgsConfig(categoryDescriptor);

            std::cout << "Creating unc configuration..." << std::endl;
            ProduceUncConfig(categoryDescriptor);

            std::cout << "Creating unc values configuration..." << std::endl;
            ProduceUncValuesConfig(categoryDescriptor);
        }
    }

private:
    void ProduceCgsConfig(const analysis::limits::CategoryDescriptor& categoryDescriptor) const
    {
        using namespace analysis::limits;
        const std::string cfgFileName = outputPath + "/cgs-Hhh-8TeV-" + categoryDescriptor.index + ".conf";
        std::ofstream cfg(cfgFileName);
        if(cfg.fail())
            throw analysis::exception("Unable to create '") << cfgFileName << "'.";
        cfg.exceptions(std::ofstream::failbit | std::ofstream::badbit);

        cfg << "# cgs.config: Specification of groups, categories, and samples\n";

        WriteSampleGroup(cfg, SampleCategory::Signal, categoryDescriptor.samples);
        WriteSampleGroup(cfg, SampleCategory::Background, categoryDescriptor.samples);
        WriteSampleGroup(cfg, SampleCategory::StandardModel, categoryDescriptor.samples);

        cfg << "\ncategories: " << categoryDescriptor.name
            << "\nsignals: " << SampleCategory::Signal
            << "\nbackgrounds: " << SampleCategory::Background
            << "\ndata: " << *categoryDescriptor.samples.GetCategorisedSamples().at(SampleCategory::Data).begin();
    }

    void WriteSampleGroup(std::ostream& cfg, analysis::limits::SampleCategory category,
                          const analysis::limits::SampleCategoryCollection& sampleCollection) const
    {
        if(!sampleCollection.GetCategorisedSamples().count(category))
            throw analysis::exception("Sample category '") << category << "' not found.";
        const auto& sampleNames = sampleCollection.GetCategorisedSamples().at(category);
        if(!sampleNames.size())
            throw analysis::exception("Sample category '") << category << "' is empty.";
        auto nameIter = sampleNames.begin();
        cfg << "$ GROUP " << category << " " << *nameIter;
        ++nameIter;
        for(; nameIter != sampleNames.end(); ++nameIter)
            cfg << "," << *nameIter;
        cfg << "\n";
    }

    void ProduceUncConfig(const analysis::limits::CategoryDescriptor& categoryDescriptor) const
    {
        using namespace analysis::limits;

        static const std::vector<int> column_widths = { 52, 0 };

        const std::string cfgFileName = outputPath + "/unc-Hhh-8TeV-" + categoryDescriptor.index + ".conf";
        std::ofstream cfg(cfgFileName);
        if(cfg.fail())
            throw analysis::exception("Unable to create '") << cfgFileName << "'.";
        cfg.exceptions(std::ofstream::failbit | std::ofstream::badbit);

        cfg << "# unc.config: Specification of uncertainty parameter names and types\n";
        cfg << std::left;
        for(const UncertaintyDescriptor* uncertaintyDescriptor : uncertainties.GetOrderedCollection()) {
            const std::string full_name = uncertaintyDescriptor->FullName(categoryDescriptor.channel_name,
                                                                          categoryDescriptor.category_name);
            cfg << std::setw(column_widths.at(0)) << full_name
                << std::setw(column_widths.at(1)) << uncertaintyDescriptor->type << "\n";
        }
    }

    void ProduceUncValuesConfig(const analysis::limits::CategoryDescriptor& categoryDescriptor)
    {
        using namespace analysis::limits;
        using analysis::EventCategory;

        static const std::map<std::string, EventCategory> category_name_map = {
            { "2jet0tag", EventCategory::TwoJets_ZeroBtag },
            { "2jet1tag", EventCategory::TwoJets_OneBtag },
            { "2jet2tag", EventCategory::TwoJets_TwoBtag }
        };
        const EventCategory eventCategory = category_name_map.at(categoryDescriptor.category_name);

        const std::string cfgFileName = outputPath + "/unc-Hhh-8TeV-" + categoryDescriptor.index + ".vals";
        std::ofstream cfg(cfgFileName);
        if(cfg.fail())
            throw analysis::exception("Unable to create '") << cfgFileName << "'.";
        cfg.exceptions(std::ofstream::failbit | std::ofstream::badbit);

        cfg << "# unc.vals: specification of uncertainty values by category, sample, and uncertainty name\n";
        cfg << std::left << std::fixed << std::setprecision(3);

        for(const UncertaintyDescriptor* uncertaintyDescriptor : uncertainties.GetOrderedCollection()) {
            const std::string full_name = uncertaintyDescriptor->FullName(categoryDescriptor.channel_name,
                                                                          categoryDescriptor.category_name);
            SampleNameSet common_samples;

            for(const std::string& sample : uncertaintyDescriptor->samples) {
                if(!categoryDescriptor.samples.GetAllSamples().count(sample)) continue;
                if(uncertaintyDescriptor->sample_values.count(sample)) {
                    const double value = uncertaintyDescriptor->sample_values.at(sample);
                    WriteUncValue(cfg, categoryDescriptor.name, sample, full_name, value);
                } else if(uncertaintyDescriptor->calculate_value) {
                    std::vector<UncertaintyInterval> unc_vector;
                    const auto samples_to_process = categoryDescriptor.samples.GenerateSampleListToProcess(sample);
                    const bool single_sample = samples_to_process.size() == 1;
                    for(const std::string& sub_sample : samples_to_process) {
                        const auto unc = calculators.Calculate(uncertaintyDescriptor->name, eventCategory, sub_sample);
                        if(!single_sample)
                            std::cout << "    " << sub_sample << ": " << uncertaintyDescriptor->name
                                      << " = " << unc << ".\n";
                        unc_vector.push_back(unc);
                    }
                    const UncertaintyInterval average_unc = UncertaintyInterval::WeightedAverage(unc_vector);
                    const auto unc_value = UncertaintyCalculatorCollection::CombineUpDownUncertainties(average_unc);
                    std::cout << sample << ": " << uncertaintyDescriptor->name << " = " << average_unc
                              << " -> " << unc_value << ".\n";
                    if(std::abs(1. - unc_value.GetValue()) >= uncertaintyDescriptor->threshold)
                        WriteUncValue(cfg, categoryDescriptor.name, sample, full_name, unc_value.GetValue());
                } else {
                    common_samples.insert(sample);
                }
            }
            if(common_samples.size()) {
                const std::string sample_list = uncertaintyDescriptor->SampleList(common_samples);
                WriteUncValue(cfg, categoryDescriptor.name, sample_list, full_name, uncertaintyDescriptor->value);
            }
        }
    }

    static void WriteUncValue(std::ostream& cfg, const std::string& category_name, const std::string& sample_list,
                              const std::string& full_unc_name, double value)
    {
        static const std::vector<int> column_widths = { 20, 36, 48, 0 };

        cfg << std::setw(column_widths.at(0)) << category_name
            << std::setw(column_widths.at(1)) << sample_list
            << std::setw(column_widths.at(2)) << full_unc_name
            << std::setw(column_widths.at(3)) << value << "\n";
    }

private:
    analysis::FlatAnalyzerDataCollectionReader anaDataReader;
    std::string outputPath;
    analysis::limits::SampleCategoryCollectionMap samples;
    analysis::limits::CategoryDescriptorMap categories;
    analysis::limits::UncertaintyDescriptorCollection uncertainties;
    analysis::DataCategoryCollection dataCategories;
    analysis::limits::UncertaintyCalculatorCollection calculators;
};
