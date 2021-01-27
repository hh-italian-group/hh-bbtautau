/*! Definition of core functionality for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/EventAnalyzerCore.h"

#include "hh-bbtautau/Analysis/include/SampleDescriptorConfigEntryReader.h"

namespace analysis {

EventAnalyzerCore::EventAnalyzerCore(const CoreAnalyzerArguments& args, Channel _channel, bool use_base_categories) :
    channelId(_channel), working_path(args.working_path()), signalObjectSelector(SignalMode::HH)
{

    ROOT::EnableThreadSafety();
    if(args.n_threads() > 1)
        ROOT::EnableImplicitMT(args.n_threads());

    ConfigReader config_reader;

    std::cout << "Loading configurations... " << std::flush;
    AnalyzerSetupCollection ana_setup_collection;
    AnalyzerConfigEntryReader ana_entry_reader(ana_setup_collection);
    config_reader.AddEntryReader("ANA_DESC", ana_entry_reader, false);

    SampleDescriptorConfigEntryReader sample_entry_reader(sample_descriptors);
    config_reader.AddEntryReader("SAMPLE", sample_entry_reader, true);

    CombinedSampleDescriptorConfigEntryReader combined_entry_reader(cmb_sample_descriptors, sample_descriptors);
    config_reader.AddEntryReader("SAMPLE_CMB", combined_entry_reader, false);

    config_reader.ReadConfig(args.sources());

    if(!ana_setup_collection.count(args.setup()))
        throw exception("Setup '%1%' not found in the configuration file '%2%'.") % args.setup() % args.sources();
    ana_setup = ana_setup_collection.at(args.setup());
    std::cout << "done.\nInitializing core variables... " << std::flush;

    signalObjectSelector = SignalObjectSelector(ana_setup.mode);
    EventRegion::Initialize(signalObjectSelector.GetTauVSjetDiscriminator().second,
                            signalObjectSelector.GetTauVSjetSidebandWPRange().first,
                            signalObjectSelector.GetTauVSjetDiscriminator().second);

    ana_setup.ConvertToEventRegion();

    RemoveUnusedSamples();

    bTagger = std::make_shared<BTagger>(ana_setup.period, ana_setup.jet_ordering);

    CreateEventSubCategoriesToProcess(use_base_categories);

    if(!ana_setup.xs_cfg.empty())
        crossSectionProvider = std::make_shared<tuple_skimmer::CrossSectionProvider>(ana_setup.xs_cfg);

    std::vector<std::string> unc_sources_group_string = SplitValueList(args.unc_sources_groups(),false, ",", true);
    for (auto& s : unc_sources_group_string){
        unc_sources_group.insert(ana_setup.unc_sources.at(s).begin(), ana_setup.unc_sources.at(s).end());
    }

}

const std::string& EventAnalyzerCore::ChannelNameLatex() const { return __Channel_names_latex.EnumToString(channelId); }

std::string EventAnalyzerCore::FullPath(const std::string& path) const
{
    return working_path.empty() ? path : working_path + "/" + path;
}

bool EventAnalyzerCore::FixNegativeContributions(TH1D& histogram, std::string& debug_info,
                                                 std::string& negative_bins_info)
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

void EventAnalyzerCore::ProcessCombinedSamples(AnaDataCollection& anaDataCollection,
                                               const EventSubCategory& subCategory,
                                               const std::vector<std::string>& sample_names)
{
    for(const std::string& sample_name : sample_names) {
        if(!cmb_sample_descriptors.count(sample_name))
            throw exception("Combined sample '%1%' not found.") % sample_name;
        CombinedSampleDescriptor& sample = cmb_sample_descriptors.at(sample_name);
        if(sample.channels.size() && !sample.channels.count(channelId)) continue;
        std::cout << "\t\t\t" << sample.name << '\n';
        for(const std::string& sub_sample_name : sample.sample_descriptors) {
            if(!sample_descriptors.count(sub_sample_name))
                throw exception("Unable to create '%1%': sub-sample '%2%' not found.")
                    % sample_name % sub_sample_name;
            SampleDescriptor& sub_sample =  sample_descriptors.at(sub_sample_name);
            AddSampleToCombined(anaDataCollection, subCategory, sample, sub_sample);
        }
    }
}

void EventAnalyzerCore::AddSampleToCombined(AnaDataCollection& anaDataCollection, const EventSubCategory& subCategory,
                                            CombinedSampleDescriptor& sample, SampleDescriptor& sub_sample)
{   for(const EventAnalyzerDataId& metaDataId : EventAnalyzerDataId::MetaLoop(ana_setup.categories,
            ana_setup.regions, unc_sources_group, GetAllUncertaintyScales())) {
        const auto anaDataId = metaDataId.Set(sample.name).Set(subCategory);
        auto& anaData = anaDataCollection.Get(anaDataId);
        for(const auto& sub_sample_wp : sub_sample.working_points) {
            const auto subDataId = metaDataId.Set(sub_sample_wp.full_name).Set(subCategory);
            auto& subAnaData = anaDataCollection.Get(subDataId);
            for(const auto& sub_entry : subAnaData.GetEntriesEx<TH1D>()) {
                auto& entry = anaData.GetEntryEx<TH1D>(sub_entry.first);
                for(const auto& hist : sub_entry.second->GetHistograms()) {
                    entry(hist.first).AddHistogram(*hist.second);
                }
            }
        }
    }

}

void EventAnalyzerCore::CreateEventSubCategoriesToProcess(bool use_base_categories)
{
    const auto& sub_categories = use_base_categories ? ana_setup.sub_categories_base : ana_setup.sub_categories;
    sub_categories_to_process.insert(sub_categories.begin(), sub_categories.end());
}

void EventAnalyzerCore::RemoveUnusedSamples()
{
    RemoveUnusedSamples(sample_descriptors, { &ana_setup.signals, &ana_setup.backgrounds, &ana_setup.data });
    RemoveUnusedSamples(cmb_sample_descriptors, { &ana_setup.cmb_samples });
}

} // namespace analysis
