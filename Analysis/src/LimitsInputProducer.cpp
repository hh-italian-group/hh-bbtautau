/*! Definition of the limits input producer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/LimitsInputProducer.h"

namespace analysis {

std::string LimitsInputProducer::FullDataCardName(const std::string& datacard_name, UncertaintySource unc_source,
                                                  UncertaintyScale unc_scale)
{
    if(unc_source == UncertaintySource::None)
        return datacard_name;

    std::ostringstream full_name;
    full_name << datacard_name << "_CMS_";
    if(unc_source == UncertaintySource::TauES)
        full_name << "scale_t";
    else if(unc_source == UncertaintySource::Total)
        full_name << "scale_j";
    else if(unc_source == UncertaintySource::TopPt)
        full_name << "topPt";
    else
        throw exception("Unsupported uncertainty source %1%.") % unc_source;
    full_name << "_13TeV" << unc_scale;
    return full_name.str();
}

std::string LimitsInputProducer::EventRegionSuffix(EventRegion region)
{
    static const std::map<EventRegion, std::string> regions {
        { EventRegion::OS_AntiIsolated(), "_OS_antiiso" },
        { EventRegion::SS_AntiIsolated(), "_SS_antiiso" },
        { EventRegion::SS_Isolated(), "_SS_iso" },
    };

    if(region == EventRegion::SignalRegion())
        return "";
    if(regions.count(region))
        return regions.at(region);
    return "_" + ToString(region);
}


void LimitsInputProducer::Produce(const std::string& outputFileNamePrefix, const std::string& setup_name,
                                  const std::map<EventCategory, std::string>& eventCategories,
                                  EventSubCategory eventSubCategory,
                                  const std::set<UncertaintySource>& uncertaintySources,
                                  const EventRegionSet& eventRegions, const std::map<SelectionCut,
                                  std::string>& sel_aliases)
{
    // static constexpr double tiny_value = 1e-9;
    // static constexpr double tiny_value_error = tiny_value;
    static const std::string dirNamePrefix = ToString(anaDataCollection->ChannelId()) + "_";

    std::ostringstream s_file_name;
    s_file_name << outputFileNamePrefix << "_" << setup_name;
    if(eventSubCategory != EventSubCategory::NoCuts())
        s_file_name << "_" << eventSubCategory.ToString(sel_aliases);
    const std::string file_name = s_file_name.str();
    auto outputFile = root_ext::CreateRootFile(file_name + ".root");
    std::set<EventAnalyzerDataId> empty_histograms;

    for(const EventAnalyzerDataId& metaId : EventAnalyzerDataId::MetaLoop(eventCategories, uncertaintySources,
                                                                          GetAllUncertaintyScales(),
                                                                          sampleWorkingPoints, eventRegions))
    {
        if(!GetActiveUncertaintyScales(metaId.Get<UncertaintySource>()).count(metaId.Get<UncertaintyScale>()))
            continue;
        const std::string directoryName = dirNamePrefix + ToString(metaId.Get<EventCategory>()) +
                                          EventRegionSuffix(metaId.Get<EventRegion>());
        TDirectory* directory = root_ext::GetDirectory(*outputFile, directoryName, true);
        const SampleWP& sampleWP = sampleWorkingPoints.at(metaId.Get<std::string>());
        const auto anaDataId = metaId.Set(eventSubCategory);
        const auto& anaData = anaDataCollection->Get(anaDataId);
        auto& hist_entry = anaData.GetEntryEx<TH1D>(eventCategories.at(anaDataId.Get<EventCategory>()));
        std::shared_ptr<TH1D> hist;
        if(hist_entry.GetHistograms().count(""))
            hist = std::make_shared<TH1D>(hist_entry());
        if(hist)
            hist->Scale(sampleWP.datacard_sf);
        if(!hist || hist->Integral() == 0.) continue;
        // {
        //     bool print_warning;
        //     if(CanHaveEmptyHistogram(anaDataId, print_warning)) continue;
        //     if(print_warning)
        //         empty_histograms.insert(anaDataId);
        //     if(!hist)
        //         hist = std::make_shared<TH1D>(hist_entry());
        //     const Int_t central_bin = hist->GetNbinsX() / 2;
        //     hist->SetBinContent(central_bin, tiny_value);
        //     hist->SetBinError(central_bin, tiny_value_error);
        // }
        const auto datacard_name = FullDataCardName(sampleWP.datacard_name, metaId.Get<UncertaintySource>(),
                                                    metaId.Get<UncertaintyScale>());
        root_ext::WriteObject(*hist, directory, datacard_name);
    }

    if(empty_histograms.size()) {
        const std::string of_name = file_name + "_emptyShapes.txt";
        std::ofstream of(of_name);
        of.exceptions(std::ios::failbit);
        for(const auto& id : empty_histograms)
            of << id << "\n";
        std::cout << "\t\t\tWarning: some datacard histograms are empty.\n"
                  << "\t\t\tThey are replaced with histograms with a tiny yield in the central bin.\n"
                  << "\t\t\tSee '" << of_name << "' for details." << std::endl;
    }
}

bool LimitsInputProducer::CanHaveEmptyHistogram(const EventAnalyzerDataId& id, bool& print_warning) const
{
    const auto& unc_source = id.Get<UncertaintySource>();
    const SampleWP& sampleWP = sampleWorkingPoints.at(id.Get<std::string>());
    print_warning = id.Get<EventRegion>() == EventRegion::SignalRegion();
    if(unc_source == UncertaintySource::None)
        return false;
    if(sampleWP.sampleType == SampleType::Data || sampleWP.sampleType == SampleType::QCD)
        return true;
    if(unc_source == UncertaintySource::TopPt)
        return sampleWP.sampleType != SampleType::TT;
    return false;
}

} // namespace analysis
