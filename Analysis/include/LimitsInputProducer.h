/*! Definition of BaseEventAnalyzer class, the base class for event analyzers.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

namespace analysis {

template<typename _FirstLeg, typename _SecondLeg>
class BaseEventAnalyzer {
public:
    using FirstLeg = _FirstLeg;
    using SecondLeg = _SecondLeg;
    using EventInfo = ::analysis::EventInfo<FirstLeg, SecondLeg>;
    using EventAnalyzerData = ::analysis::EventAnalyzerData<FirstLeg, SecondLeg>;
    using PhysicalValueMap = std::map<EventRegion, PhysicalValue>;

    std::string FullDataCardName(const std::string& datacard_name, EventEnergyScale eventEnergyScale) const
    {
        if(eventEnergyScale == EventEnergyScale::Central)
            return datacard_name;

        std::string channel_name = ChannelName();
        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);
        std::ostringstream full_name;
        full_name << datacard_name << "_CMS_scale_";
        if(eventEnergyScale == EventEnergyScale::TauUp || eventEnergyScale == EventEnergyScale::TauDown)
            full_name << "t_" << channel_name;
        else if(eventEnergyScale == EventEnergyScale::JetUp || eventEnergyScale == EventEnergyScale::JetDown)
            full_name << "j";
        else if(eventEnergyScale == EventEnergyScale::BtagEfficiencyUp
                || eventEnergyScale == EventEnergyScale::BtagEfficiencyDown)
            full_name << "btagEff";
        else if(eventEnergyScale == EventEnergyScale::BtagFakeUp || eventEnergyScale == EventEnergyScale::BtagFakeDown)
            full_name << "btagFake";
        else
            throw exception("Unsupported event energy scale %1%.") % eventEnergyScale;
        full_name << "_13TeV";
        if(eventEnergyScale == EventEnergyScale::TauUp || eventEnergyScale == EventEnergyScale::JetUp ||
                eventEnergyScale == EventEnergyScale::BtagEfficiencyUp ||
                eventEnergyScale == EventEnergyScale::BtagFakeUp )
            full_name << "Up";
        else
            full_name << "Down";
        return full_name.str();
    }

    void ProduceFileForLimitsCalculation(const std::string& hist_name, EventSubCategory eventSubCategory,
                                         typename EventAnalyzerData::HistogramAccessor histogramAccessor)
    {
        static const std::map<EventCategory, std::string> categoryToDirectoryNameSuffix = {
            { EventCategory::Inclusive, "inclusive" }, { EventCategory::TwoJets_ZeroBtag, "2jet0tag" },
            { EventCategory::TwoJets_OneBtag, "2jet1tag" }, { EventCategory::TwoJets_TwoBtag, "2jet2tag" }
        };

        static const std::map<std::string, std::string> channelNameForFolder = {
            { "eTau", "eleTau" }, { "muTau", "muTau" }, { "tauTau", "tauTau" }
        };

        static const double tiny_value = 1e-9;
        static const double tiny_value_error = tiny_value;

        std::ostringstream s_file_name;
        s_file_name << args.outputFileName() << "_" << hist_name;
        if(eventSubCategory != EventSubCategory::NoCuts)
            s_file_name << "_" << eventSubCategory;
        s_file_name << ".root";
        const std::string file_name = s_file_name.str();

        auto outputFile = root_ext::CreateRootFile(file_name);
        for(EventCategory eventCategory : EventCategoriesToProcess()) {
            if(!categoryToDirectoryNameSuffix.count(eventCategory)) continue;
            const std::string directoryName = channelNameForFolder.at(ChannelName()) + "_"
                    + categoryToDirectoryNameSuffix.at(eventCategory);
            outputFile->mkdir(directoryName.c_str());
            TDirectory* directory = outputFile->GetDirectory(directoryName.c_str());
            for(const DataCategory* dataCategory : dataCategoryCollection.GetCategories(DataCategoryType::Limits)) {
                if(!dataCategory->datacard.size())
                    throw exception("Empty datacard name for data category '%1%'.") % dataCategory->name;
                for(const EventEnergyScale& eventEnergyScale : EventEnergyScaleToProcess()) {
                    const EventAnalyzerDataMetaId_noRegion_noName meta_id(eventCategory, eventSubCategory,
                                                                         eventEnergyScale);
                    std::shared_ptr<TH1D> hist;
                    if(auto hist_orig = GetSignalHistogram(meta_id, dataCategory->name, hist_name))
                        hist = std::shared_ptr<TH1D>(new TH1D(*hist_orig));
                    else {
                        std::cout << "Warning - Datacard histogram '" << hist_name
                                  << "' not found for data category '" << dataCategory->name << "' in '"
                                  << eventCategory << "/" << eventSubCategory << "/" << eventEnergyScale
                                  << "'. Using histogram with a tiny yield in the central bin instead.\n";

                        EventAnalyzerData& anaData = GetAnaData(meta_id.MakeId(EventRegion::OS_Isolated,
                                                                              dataCategory->name));
                        auto& new_hist = (anaData.*histogramAccessor)();
                        hist = std::shared_ptr<TH1D>(new TH1D(new_hist));
                        const Int_t central_bin = hist->GetNbinsX() / 2;
                        hist->SetBinContent(central_bin, tiny_value);
                        hist->SetBinError(central_bin, tiny_value_error);
                    }
                    const std::string full_datacard_name = FullDataCardName(dataCategory->datacard, eventEnergyScale);
                    hist->Scale(dataCategory->limits_sf);
                    root_ext::WriteObject(*hist, directory, full_datacard_name);

                    if(eventEnergyScale == EventEnergyScale::Central && dataCategory->datacard == "ZL") {
                        std::string channel_name = ChannelName();
                        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);
                        const std::string name_syst_prefix = dataCategory->datacard + "_CMS_htt_"
                                + dataCategory->datacard + "Scale_" + channel_name + "_8TeV";
                        const std::string name_syst_up = name_syst_prefix + "Up";
                        const std::string name_syst_down = name_syst_prefix + "Down";
                        std::shared_ptr<TH1D> hist_syst_up(new TH1D(*hist));
                        hist_syst_up->Scale(1.02);
                        root_ext::WriteObject(*hist_syst_up, directory, name_syst_up);
                        std::shared_ptr<TH1D> hist_syst_down(new TH1D(*hist));
                        hist_syst_down->Scale(0.98);
                        root_ext::WriteObject(*hist_syst_down, directory, name_syst_down);
                    }
                    //added shape systematics for QCD
                    if(eventEnergyScale == EventEnergyScale::Central && dataCategory->datacard == "QCD_alternative") {
                        std::string channel_name = ChannelName();
                        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);
                        const std::string name_syst_prefix = "QCD_CMS_htt_" + dataCategory->datacard
                                + "Shape_" + channel_name + "_8TeV";
                        const std::string name_syst_up = name_syst_prefix + "Up";
                        const std::string name_syst_down = name_syst_prefix + "Down";
                        std::shared_ptr<TH1D> hist_syst_up(new TH1D(*hist));
                        root_ext::WriteObject(*hist_syst_up, directory, name_syst_up);
                        std::shared_ptr<TH1D> hist_syst_down(new TH1D(*hist));
                        root_ext::WriteObject(*hist_syst_down, directory, name_syst_down);
                    }
                }
            }
        }
    }
};

} // namespace analysis
