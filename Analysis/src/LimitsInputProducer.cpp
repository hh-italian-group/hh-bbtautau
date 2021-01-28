/*! Definition of the limits input producer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/LimitsInputProducer.h"

namespace analysis {

std::string LimitsInputProducer::FullDataCardName(std::string datacard_name, UncertaintySource unc_source,
                                                  UncertaintyScale unc_scale, Period period)
{
    boost::replace_all(datacard_name, ".", "p");
    if(unc_source == UncertaintySource::None)
        return datacard_name;

    std::ostringstream full_name;
    full_name << datacard_name << "_CMS_";
    full_name << UncSourceSuffix(unc_source, period);
    /*
    if(unc_source == UncertaintySource::TauES)
        full_name << "scale_t";
    else if(unc_source == UncertaintySource::JetReduced_Total)
        full_name << "scale_j";
    else if(unc_source == UncertaintySource::TopPt)
        full_name << "topPt";
    else if(unc_source == UncertaintySource::TopPt)
        full_name << UncSourceSuffix(unc_source);
    else
        throw exception("Unsupported uncertainty source %1%.") % unc_source;
         */
    // << static_cast<int>(period)
    full_name  << "_13TeV"  << unc_scale;

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

std::string LimitsInputProducer::UncSourceSuffix(UncertaintySource unc_source, Period period)
{
    std::map<UncertaintySource, std::string> unc_sources {
        { UncertaintySource::JetReduced_Absolute, "JES_Abs" },
        { UncertaintySource::JetReduced_BBEC1, "JES_BBEC1"},
        { UncertaintySource::JetReduced_EC2, "JES_EC2" },
        { UncertaintySource::JetReduced_FlavorQCD, "JES_FlavQCD" },
        { UncertaintySource::JetReduced_HF, "JES_HF" },
        { UncertaintySource::JetReduced_RelativeBal, "JES_RelBal" },
        { UncertaintySource::TauES, "scale_t"},
        { UncertaintySource::JetReduced_Total, "res_j"},
        { UncertaintySource::TauES_DM0, "scale_t_dm0"},
        { UncertaintySource::TauES_DM1, "scale_t_dm1"},
        { UncertaintySource::TauES_DM10, "scale_t_dm10"},
        { UncertaintySource::TauES_DM11, "scale_t_dm11"},
        { UncertaintySource::EleFakingTauES_DM0, "e_FakeTau_dm0"},
        { UncertaintySource::EleFakingTauES_DM1, "e_FakeTau_dm1"},
        { UncertaintySource::MuFakingTauES, "m_FakeTau_dm0"},
        { UncertaintySource::EleTriggerUnc, "eff_trigger_e"},
        { UncertaintySource::MuonTriggerUnc, "eff_trigger_m"},
        { UncertaintySource::TauTriggerUnc_DM0, "eff_trigger_t_dm0"},
        { UncertaintySource::TauTriggerUnc_DM1, "eff_trigger_t_dm1"},
        { UncertaintySource::TauTriggerUnc_DM10, "eff_trigger_t_dm10"},
        { UncertaintySource::TauTriggerUnc_DM11, "eff_trigger_t_dm11"},
        { UncertaintySource::TauVSjetSF_DM0, "eff_t_dm0"},
        { UncertaintySource::TauVSjetSF_DM1, "eff_t_dm1"},
        { UncertaintySource::TauVSjetSF_3prong, "eff_t_3prong"},
        { UncertaintySource::TauVSjetSF_pt20to25, "eff_t_pt20to25"},
        { UncertaintySource::TauVSjetSF_pt25to30, "eff_t_pt25to30"},
        { UncertaintySource::TauVSjetSF_pt30to35, "eff_t_pt30to35"},
        { UncertaintySource::TauVSjetSF_pt35to40, "eff_t_pt35to40"},
        { UncertaintySource::TauVSjetSF_ptgt40, "eff_t_ptgt40"},
        { UncertaintySource::TauVSeSF_barrel, "e_FakeTau_barrel"},
        { UncertaintySource::TauVSeSF_endcap, "e_FakeTau_endcap"},
        { UncertaintySource::TauVSmuSF_etaLt0p4, "m_FakeTau_etalt0p4"},
        { UncertaintySource::TauVSmuSF_eta0p4to0p8, "m_FakeTau_eta0p4to0p8"},
        { UncertaintySource::TauVSmuSF_eta0p8to1p2, "m_FakeTau_eta0p8to1p2"},
        { UncertaintySource::TauVSmuSF_eta1p2to1p7, "m_FakeTau_eta1p2to1p7"},
        { UncertaintySource::TauVSmuSF_etaGt1p7, "m_FakeTau_etagt1p7"},
        { UncertaintySource::EleIdIsoUnc,  "eff_e"},
        { UncertaintySource::MuonIdIsoUnc, "eff_mu"},
        { UncertaintySource::TopPt,  "shape_topPt"},
        { UncertaintySource::L1_prefiring, "prefiring"},
        { UncertaintySource::PileUp, "PU"},
        { UncertaintySource::PileUpJetId_eff, "PUJET_ID_eff"},
        { UncertaintySource::PileUpJetId_mistag, "PUJET_ID_mistag"},
        { UncertaintySource::TauCustomSF_DM0,  "hbbhtt_eff_t_dm0"},
        { UncertaintySource::TauCustomSF_DM1, "hbbhtt_eff_t_dm1"},
        { UncertaintySource::TauCustomSF_DM10, "hbbhtt_eff_t_dm10"},
        { UncertaintySource::TauCustomSF_DM11, "hbbhtt_eff_t_dm11"},
        { UncertaintySource::VBFTriggerUnc_jets, "hbbhtt_eff_VBFtrigger_j" },
        { UncertaintySource::VBFTauTriggerUnc_DM0, "hbbhtt_eff_VBFtrigger_t_DM0" },
        { UncertaintySource::VBFTauTriggerUnc_DM1, "hbbhtt_eff_VBFtrigger_t_DM1" },
        { UncertaintySource::VBFTauTriggerUnc_3prong, "hbbhtt_eff_VBFtrigger_t_3prong" },
        { UncertaintySource::btag_lf, "btag_LF" },
        { UncertaintySource::btag_hf, "btag_HF" },
        { UncertaintySource::btag_hfstats1, "btag_hfstats1" },
        { UncertaintySource::btag_hfstats2, "btag_hfstats2" },
        { UncertaintySource::btag_lfstats1, "btag_lfstats1" },
        { UncertaintySource::btag_lfstats2, "btag_lfstats2" },
        { UncertaintySource::btag_cferr1, "btag_cferr1" },
        { UncertaintySource::btag_cferr2, "btag_cferr2" },
    };
    if(period==Period::Run2016){
        unc_sources.emplace(UncertaintySource::JetReduced_Absolute_year , "JES_Abs_2016");
        unc_sources.emplace(UncertaintySource::JetReduced_BBEC1_year , "JES_BBEC1_2016");
        unc_sources.emplace(UncertaintySource::JetReduced_EC2_year , "JES_EC2_2016");
        unc_sources.emplace(UncertaintySource::JetReduced_HF_year , "JES_HF_2016");
        unc_sources.emplace(UncertaintySource::JetReduced_RelativeSample_year , "JES_RelSample_2016");
        unc_sources.emplace(UncertaintySource::JetReduced_Absolute_year , "JES_Abs_2016");
    }
    else if(period==Period::Run2017){
        unc_sources.emplace(UncertaintySource::JetReduced_Absolute_year , "JES_Abs_2017");
        unc_sources.emplace(UncertaintySource::JetReduced_BBEC1_year , "JES_BBEC1_2017");
        unc_sources.emplace(UncertaintySource::JetReduced_EC2_year , "JES_EC2_2017");
        unc_sources.emplace(UncertaintySource::JetReduced_HF_year , "JES_HF_2017");
        unc_sources.emplace(UncertaintySource::JetReduced_RelativeSample_year , "JES_RelSample_2017");
        unc_sources.emplace(UncertaintySource::JetReduced_Absolute_year , "JES_Abs_2017");
    }
    else if(period==Period::Run2018){
        unc_sources.emplace(UncertaintySource::JetReduced_Absolute_year , "JES_Abs_2018");
        unc_sources.emplace(UncertaintySource::JetReduced_BBEC1_year , "JES_BBEC1_2018");
        unc_sources.emplace(UncertaintySource::JetReduced_EC2_year , "JES_EC2_2018");
        unc_sources.emplace(UncertaintySource::JetReduced_HF_year , "JES_HF_2018");
        unc_sources.emplace(UncertaintySource::JetReduced_RelativeSample_year , "JES_RelSample_2018");
        unc_sources.emplace(UncertaintySource::JetReduced_Absolute_year , "JES_Abs_2018");
    }
    if(unc_sources.count(unc_source))
        return unc_sources.at(unc_source);
    else
        return analysis::ToString(unc_source);
}

std::pair<std::string,std::string> LimitsInputProducer::ProdCatSuffix(std::string category)
{
    static const std::map<std::string, std::pair<std::string,std::string>> categories {
        {"2j1bR_noVBF", {"ggHH", "res1b"}},
        {"2j2b+R_noVBF", {"ggHH", "res2b"}},
        {"2j2Lb+B_noVBF", {"ggHH", "boosted"}},
        {"2j1b+_VBF_qqHH", {"qqHH", "classGGF"}},
        {"2j1b+_VBF_ggHH", {"qqHH", "classVBF"}},
        {"2j1b+_VBF_TT", {"qqHH", "classTT"}},
        {"2j1b+_VBF_ttH", {"qqHH", "classttH"}},
        {"2j1b+_VBF_DY", {"qqHH", "classDY"}},
    };

    if(categories.count(category))
        return categories.at(category);
    else
        return std::make_pair("","");
}


void LimitsInputProducer::Produce(const std::string& outputFileNamePrefix, const std::string& setup_name,
                                  const std::map<EventCategory, std::string>& eventCategories,
                                  EventSubCategory eventSubCategory,
                                  const std::set<UncertaintySource>& uncertaintySources,
                                  const EventRegionSet& eventRegions, const std::map<SelectionCut,
                                  std::string>& sel_aliases, Period period)
{
    //static constexpr double tiny_value = 1e-9;
    //static constexpr double tiny_value_error = tiny_value;
    //productionmode_decaychannel_year_channel_category
    // e.g. ggHH_hbbhtt_2016_muTau_res1b
    static const std::string dirNamePrefix = boost::str(boost::format("_hbbhtt_%1%_%2%_")
            %  anaDataCollection->ChannelId() % static_cast<int>(period));

    std::ostringstream s_file_name;
    s_file_name << outputFileNamePrefix << "_" << setup_name;
    if(eventSubCategory != EventSubCategory::NoCuts())
        s_file_name << "_" << eventSubCategory.ToString(sel_aliases);
    const std::string file_name = s_file_name.str();
    auto outputFile = root_ext::CreateRootFile(file_name + ".root",  ROOT::kLZMA, 9);
    std::set<EventAnalyzerDataId> empty_histograms;

    for(const EventAnalyzerDataId& metaId : EventAnalyzerDataId::MetaLoop(eventCategories, uncertaintySources,
                                                                          GetAllUncertaintyScales(),
                                                                          sampleWorkingPoints, eventRegions))
    {
        if(!GetActiveUncertaintyScales(metaId.Get<UncertaintySource>()).count(metaId.Get<UncertaintyScale>()))
            continue;
        const SampleWP& sampleWP = sampleWorkingPoints.at(metaId.Get<std::string>());
        if(sampleWP.datacard_name.empty()) continue;
        const std::string directoryName = ProdCatSuffix(ToString(metaId.Get<EventCategory>())).first+ dirNamePrefix + ProdCatSuffix(ToString(metaId.Get<EventCategory>())).second+
                                          EventRegionSuffix(metaId.Get<EventRegion>());

        TDirectory* directory = root_ext::GetDirectory(*outputFile, directoryName, true);
        const auto anaDataId = metaId.Set(eventSubCategory);
        const auto& anaData = anaDataCollection->Get(anaDataId);
        auto& hist_entry = anaData.GetEntryEx<TH1D>(eventCategories.at(anaDataId.Get<EventCategory>()));
        std::shared_ptr<TH1D> hist;
        if(hist_entry.GetHistograms().count("")|| sampleWP.datacard_name == "data_obs")
            hist = std::make_shared<TH1D>(hist_entry());
        if(hist)
            hist->Scale(sampleWP.datacard_sf);
        if(!(hist && (hist->Integral() > 0. || sampleWP.datacard_name == "data_obs"))) continue;
        //{
        //     bool print_warning;
        //     if(CanHaveEmptyHistogram(anaDataId, print_warning)) continue;
        //     if(print_warning)
        //         empty_histograms.insert(anaDataId);
        //     if(!hist)
        //         hist = std::make_shared<TH1D>(hist_entry());
        //     const Int_t central_bin = hist->GetNbinsX() / 2;
        //     hist->SetBinContent(central_bin, tiny_value);
        //     hist->SetBinError(central_bin, tiny_value_error);
        //}
        const auto datacard_name = FullDataCardName(sampleWP.datacard_name, metaId.Get<UncertaintySource>(),
                                                    metaId.Get<UncertaintyScale>(), period);
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
