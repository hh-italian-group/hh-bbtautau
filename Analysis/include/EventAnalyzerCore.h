/*! Definition of core functionality for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "SampleDescriptorConfigEntryReader.h"

namespace analysis {

struct CoreAnalyzerArguments {
    REQ_ARG(std::string, sources);
    REQ_ARG(std::string, setup);
    OPT_ARG(std::string, mva_setup, "");
    OPT_ARG(std::string, working_path, "");
    OPT_ARG(unsigned, n_threads, 1);

    CoreAnalyzerArguments() {}
    CoreAnalyzerArguments(const CoreAnalyzerArguments&) = default;
    virtual ~CoreAnalyzerArguments() {}
};

class EventAnalyzerCore {
public:
    EventAnalyzerCore(const CoreAnalyzerArguments& args, Channel _channel) :
        channelId(_channel), working_path(args.working_path())
    {
        ROOT::EnableThreadSafety();
        if(args.n_threads() > 1)
            ROOT::EnableImplicitMT(args.n_threads());

        ConfigReader config_reader;

        AnalyzerSetupCollection ana_setup_collection;
        AnalyzerConfigEntryReader ana_entry_reader(ana_setup_collection);
        config_reader.AddEntryReader("ANA_DESC", ana_entry_reader, false);

        MvaReaderSetupCollection mva_setup_collection;
        MvaReaderSetupEntryReader mva_entry_reader(mva_setup_collection);
        config_reader.AddEntryReader("MVA", mva_entry_reader, false);

        SampleDescriptorConfigEntryReader sample_entry_reader(sample_descriptors);
        config_reader.AddEntryReader("SAMPLE", sample_entry_reader, true);

        CombinedSampleDescriptorConfigEntryReader combined_entry_reader(cmb_sample_descriptors, sample_descriptors);
        config_reader.AddEntryReader("SAMPLE_CMB", combined_entry_reader, false);

        config_reader.ReadConfig(args.sources());

        if(!ana_setup_collection.count(args.setup()))
            throw exception("Setup '%1%' not found in the configuration file '%2%'.") % args.setup() % args.sources();
        ana_setup = ana_setup_collection.at(args.setup());

        if(args.mva_setup().size()) {
            const auto mva_setup_names = SplitValueList(args.mva_setup(), false, ", \t", true);
            std::vector<MvaReaderSetup> mva_setups;
            for(const auto& name : mva_setup_names) {
                if(!mva_setup_collection.count(name))
                    throw exception("MVA setup '%1%' not found.") % name;
                mva_setups.push_back(mva_setup_collection.at(name));
            }
            mva_setup = mva_setups.size() == 1 ? mva_setups.front() : MvaReaderSetup::Join(mva_setups);
            CreateMvaSelectionAliases();
        }
        RemoveUnusedSamples();

        CreateEventSubCategoriesToProcess();
    }

    EventAnalyzerCore(const EventAnalyzerCore&) = delete;
    EventAnalyzerCore& operator=(const EventAnalyzerCore&) = delete;
    virtual ~EventAnalyzerCore() {}

    const std::string& ChannelNameLatex() const { return __Channel_names_latex.EnumToString(channelId); }

    std::string FullPath(const std::string& path) const
    {
        return working_path.empty() ? path : working_path + "/" + path;
    }

    static bool FixNegativeContributions(TH1D& histogram, std::string& debug_info, std::string& negative_bins_info)
    {
        static const double correction_factor = 0.0000001;

        std::ostringstream ss_debug;

        const PhysicalValue original_Integral = Integral(histogram, true);
        ss_debug << "\nSubtracted hist for '" << histogram.GetName() << ".\n";
        ss_debug << "Integral after bkg subtraction: " << original_Integral << ".\n";
        debug_info = ss_debug.str();
        if (original_Integral.GetValue() < 0) {
            std::cout << debug_info << std::endl;
            std::cout << "Integral after bkg subtraction is negative for histogram '"
                      << histogram.GetName() << std::endl;
            return false;
        }

        std::ostringstream ss_negative;

        for (Int_t n = 1; n <= histogram.GetNbinsX(); ++n) {
            if (histogram.GetBinContent(n) >= 0) continue;
            const std::string prefix = histogram.GetBinContent(n) + histogram.GetBinError(n) >= 0 ? "WARNING" : "ERROR";

            ss_negative << prefix << ": " << histogram.GetName() << " Bin " << n << ", content = "
                        << histogram.GetBinContent(n) << ", error = " << histogram.GetBinError(n)
                        << ", bin limits=[" << histogram.GetBinLowEdge(n) << "," << histogram.GetBinLowEdge(n+1)
                        << "].\n";
            const double error = correction_factor - histogram.GetBinContent(n);
            const double new_error = std::sqrt(std::pow(error,2) + std::pow(histogram.GetBinError(n),2));
            histogram.SetBinContent(n, correction_factor);
            histogram.SetBinError(n, new_error);
        }
        analysis::RenormalizeHistogram(histogram, original_Integral.GetValue(), true);
        negative_bins_info = ss_negative.str();
        return true;
    }

private:
    void CreateEventSubCategoriesToProcess()
    {
        sub_categories_to_process.insert(ana_setup.sub_categories.begin(), ana_setup.sub_categories.end());
        if(mva_setup.is_initialized()) {
            for(const auto& base : ana_setup.sub_categories) {
                SelectionCut predefined_cut;
                if(base.TryGetLastMvaCut(predefined_cut)) continue;
                for(const auto& mva_sel : mva_setup->selections) {
                    auto param = mva_sel.second;
                    std::cout<<"Selection cut: "<<mva_sel.first<<" - name: "<<param.name
                            <<" spin: "<<param.spin<<" mass: "<<param.mass<<" cut: "<<param.cut<<std::endl;
                    auto sub_category = EventSubCategory(base).SetCutResult(mva_sel.first, true);
                    sub_categories_to_process.insert(sub_category);
                }
            }
        }
    }

    void CreateMvaSelectionAliases()
    {
        mva_sel_aliases.clear();
        if(!mva_setup) return;
        for(const auto& entry : mva_setup->selections)
        {
            const SelectionCut sel = entry.first;
            const auto& params = entry.second;
            std::ostringstream ss_name;
            ss_name << params.name << "S" << params.spin << "M" << params.mass << "C" << params.cut;
            const std::string name = ss_name.str();
            for(const auto& alias : mva_sel_aliases) {
                if(alias.second == name)
                    throw exception("Duplicated mva selection alias = '%1%'.") % name;
            }
            mva_sel_aliases[sel] = name;
        }
    }

    template<typename SampleCollection>
    std::vector<std::string> FilterInactiveSamples(const SampleCollection& samples,
                                                   const std::vector<std::string>& sample_names) const
    {
        std::vector<std::string> filtered;
        for(const auto& name : sample_names) {
            if(!samples.count(name))
                throw exception("Sample '%1%' not found.") % name;
            const auto& sample = samples.at(name);
            if(sample.channels.size() && !sample.channels.count(channelId)) continue;
            filtered.push_back(name);
        }
        return filtered;
    }

    void RemoveUnusedSamples()
    {
        RemoveUnusedSamples(sample_descriptors, { &ana_setup.signals, &ana_setup.backgrounds, &ana_setup.data });
        RemoveUnusedSamples(cmb_sample_descriptors, { &ana_setup.cmb_samples });
    }

    template<typename SampleCollection>
    void RemoveUnusedSamples(SampleCollection& samples,
                             const std::vector< std::vector<std::string>*> selected_sample_names) const
    {
        std::set<std::string> used_samples;
        for(auto selected_collection : selected_sample_names) {
            *selected_collection = FilterInactiveSamples(samples, *selected_collection);
            used_samples.insert(selected_collection->begin(), selected_collection->end());
        }

        const std::set<std::string> all_sample_names = tools::collect_map_keys(samples);
        for(const auto& sample_name : all_sample_names) {
            if(!used_samples.count(sample_name))
                samples.erase(sample_name);
        }
    }

protected:
    AnalyzerSetup ana_setup;
    boost::optional<MvaReaderSetup> mva_setup;
    SampleDescriptorCollection sample_descriptors;
    CombinedSampleDescriptorCollection cmb_sample_descriptors;
    EventSubCategorySet sub_categories_to_process;
    Channel channelId;
    std::map<SelectionCut, std::string> mva_sel_aliases;
    std::string working_path;
};

} // namespace analysis
