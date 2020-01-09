/*! Definition of the file configuration entry reader for analyzer.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/SampleDescriptorConfigEntryReader.h"

namespace analysis {

void AnalyzerConfigEntryReader::EndEntry()
{
    CheckReadParamCounts("int_lumi", 1, Condition::less_equal);
    CheckReadParamCounts("period", 1, Condition::less_equal);
    CheckReadParamCounts("mode", 1, Condition::less_equal);
    CheckReadParamCounts("qcd_method", 1, Condition::less_equal);
    CheckReadParamCounts("tauID_wp", 1, Condition::less_equal);
    CheckReadParamCounts("final_variables", 1, Condition::less_equal);
    CheckReadParamCounts("apply_mass_cut", 1, Condition::less_equal);
    CheckReadParamCounts("apply_os_cut", 1, Condition::less_equal);
    CheckReadParamCounts("apply_iso_cut", 1, Condition::less_equal);
    CheckReadParamCounts("use_kinFit", 1, Condition::less_equal);
    CheckReadParamCounts("use_svFit", 1, Condition::less_equal);
    CheckReadParamCounts("applyTauID", 1, Condition::less_equal);
    CheckReadParamCounts("unc_sources", 1, Condition::less_equal);
    CheckReadParamCounts("categories", 1, Condition::less_equal);
    CheckReadParamCounts("sub_categories", 1, Condition::less_equal);
    CheckReadParamCounts("regions", 1, Condition::less_equal);
    CheckReadParamCounts("data", 1, Condition::less_equal);
    CheckReadParamCounts("signals", 1, Condition::less_equal);
    CheckReadParamCounts("backgrounds", 1, Condition::less_equal);
    CheckReadParamCounts("cmb_samples", 1, Condition::less_equal);
    CheckReadParamCounts("draw_sequence", 1, Condition::less_equal);
    CheckReadParamCounts("mva_setup", 1, Condition::less_equal);
    CheckReadParamCounts("hist_cfg", 1, Condition::less_equal);
    CheckReadParamCounts("trigger_path", 1, Condition::less_equal);
    CheckReadParamCounts("syncDataIds", 1, Condition::less_equal);
    CheckReadParamCounts("plot_cfg", 1, Condition::less_equal);
    CheckReadParamCounts("plot_page_opt", 1, Condition::less_equal);
    CheckReadParamCounts("massWindowParams", 0, Condition::greater_equal);
    CheckReadParamCounts("trigger", 0, Condition::greater_equal);
    CheckReadParamCounts("unc_cfg", 1, Condition::less_equal);
    CheckReadParamCounts("jet_ordering", 1, Condition::less_equal);
    CheckReadParamCounts("limit_setup", 0, Condition::greater_equal);
    CheckReadParamCounts("qcd_ss_os_sf",1,Condition::less_equal);
    CheckReadParamCounts("qcd_ss_os_err",1, Condition::less_equal);
    current.CreateLimitSetups();
    ConfigEntryReaderT<AnalyzerSetup>::EndEntry();
}

void AnalyzerConfigEntryReader::ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                                              std::istringstream& /*ss*/)
{
    ParseEntry("int_lumi", current.int_lumi);
    ParseEntry("period", current.period);
    ParseEntry("mode", current.mode);
    ParseEntry("qcd_method", current.qcd_method);
    ParseEntry("tauID_wp", current.tauID_wp);
    ParseEntryList("final_variables", current.final_variables);
    ParseEntry("apply_mass_cut", current.apply_mass_cut);
    ParseEntry("apply_os_cut", current.apply_os_cut);
    ParseEntry("apply_iso_cut", current.apply_iso_cut);
    ParseEntry("use_kinFit", current.use_kinFit);
    ParseEntry("use_svFit", current.use_svFit);
    ParseEntry("applyTauID", current.applyTauID);
    ParseEnumList("unc_sources", current.unc_sources);
    ParseEnumList("categories", current.categories);
    ParseEnumList("sub_categories", current.sub_categories);
    ParseEnumList("regions", current.regions);
    ParseEntryList("data", current.data);
    ParseEntryList("signals", current.signals);
    ParseEntryList("backgrounds", current.backgrounds);
    ParseEntryList("cmb_samples", current.cmb_samples);
    ParseEntryList("draw_sequence", current.draw_sequence);
    ParseEntry("limit_category", current.limit_categories);
    ParseEntry("mva_setup", current.mva_setup);
    ParseEntry("hist_cfg", current.hist_cfg);
    ParseEntry("trigger_path", current.trigger_path);
    ParseEntryList("syncDataIds", current.syncDataIds);
    ParseEntry("plot_cfg", current.plot_cfg);
    ParseEntry("plot_page_opt", current.plot_page_opt);
    ParseEntry("massWindowParams", current.massWindowParams);
    ParseMappedEntryList("trigger", current.trigger,false);
    ParseEntry("unc_cfg", current.unc_cfg);
    ParseEntry("jet_ordering", current.jet_ordering);
    ParseMappedEntryList("limit_setup", current.limit_setup_raw,false);
    ParseEntry("qcd_ss_os_sf",current.qcd_ss_os_sf);
    ParseEntry("qcd_ss_os_err",current.qcd_ss_os_err);
}

void MvaReaderSetupEntryReader::EndEntry()
{
    CheckReadParamCounts("training", 0, Condition::greater_equal);
    CheckReadParamCounts("variables", 0, Condition::greater_equal);
    CheckReadParamCounts("masses", 0, Condition::greater_equal);
    CheckReadParamCounts("spins", 0, Condition::greater_equal);
    CheckReadParamCounts("cuts", 0, Condition::greater_equal);
    CheckReadParamCounts("legacy", 0, Condition::greater_equal);
    CheckReadParamCounts("training_range", 0, Condition::greater_equal);
    CheckReadParamCounts("samples", 0, Condition::greater_equal);

    current.CreateSelections();
    ConfigEntryReaderT<MvaReaderSetup>::EndEntry();
}

void MvaReaderSetupEntryReader::ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                           std::istringstream& /*ss*/)
{
    ParseEntry("training", current.trainings);
    ParseMappedEntryList("variables", current.variables, false);
    ParseMappedEntryList("masses", current.masses, true);
    ParseMappedEntryList("spins", current.spins, true);
    ParseMappedEntryList("cuts", current.cuts, false);
    ParseEntry("training_range", current.training_ranges);
    ParseMappedEntryList("samples", current.samples, false);
    ParseEntry("legacy", current.legacy);
}

void SampleDescriptorConfigEntryReader::EndEntry()
{
    CheckReadParamCounts("name_suffix", 1, Condition::less_equal);
    CheckReadParamCounts("file_path", 1, Condition::less_equal);
    CheckReadParamCounts("cross_section", 1, Condition::less_equal);
    CheckReadParamCounts("points", 0, Condition::greater_equal);
    CheckReadParamCounts("draw_ex", 0, Condition::greater_equal);
    CheckReadParamCounts("norm_sf", 1, Condition::less_equal);
    CheckReadParamCounts("reference_pu_sample", 1, Condition::less_equal);


    Base::EndEntry();
}

void SampleDescriptorConfigEntryReader::ReadParameter(const std::string& param_name, const std::string& param_value,
                                                      std::istringstream& ss)
{
    ParseEntry("name_suffix", current.name_suffix);
    ParseEntry("file_path", current.file_path);
    ParseEntry<double, NumericalExpression>("cross_section", current.cross_section,
                                            [](double xs){ return xs > 0; });
    ParseMappedEntryList("points", current.points, true);
    ParseEntry("draw_ex", current.draw_ex);
    ParseEntryList("norm_sf", current.norm_sf, true);
    ParseEntry("reference_pu_sample", current.reference_pu_sample);

    Base::ReadParameter(param_name,param_value,ss);
}


CombinedSampleDescriptorConfigEntryReader::CombinedSampleDescriptorConfigEntryReader(
        CombinedSampleDescriptorCollection& _descriptors,
        const SampleDescriptorCollection& _sampleDescriptors) :
   Base(_descriptors), sampleDescriptorCollection(&_sampleDescriptors)
{
}

void CombinedSampleDescriptorConfigEntryReader::EndEntry()
{
    CheckReadParamCounts("sample_descriptors", 1, Condition::equal_to);
    Base::EndEntry();
}

void CombinedSampleDescriptorConfigEntryReader::ReadParameter(const std::string& param_name,
                                                              const std::string& param_value,
                                                              std::istringstream& ss)
{
    ParseEntryList("sample_descriptors", current.sample_descriptors, false, " \t",
               [&](const std::string& name){return sampleDescriptorCollection->count(name);});

    Base::ReadParameter(param_name,param_value,ss);
}

ModellingUncertaintyEntryReader::ModellingUncertaintyEntryReader(ModellingUncertaintyCollection& _items) : items(&_items) {}

void ModellingUncertaintyEntryReader::StartEntry(const std::string& name, const std::string& reference_name)
{
    ConfigEntryReader::StartEntry(name, reference_name);
    current = ModellingUncertainty();
    current_name = name;
}

void ModellingUncertaintyEntryReader::EndEntry()
{
    CheckReadParamCounts("unc", 1, Condition::less_equal);
    CheckReadParamCounts("sf", 1, Condition::less_equal);
    CheckReadParamCounts("ref_category", 1, Condition::less_equal);

    current.CreateSampleUncMap();
    items->Add(current_name, current);
}

void ModellingUncertaintyEntryReader::ReadParameter(const std::string& /*param_name*/,
                                                    const std::string& /*param_value*/,
                                                    std::istringstream& /*ss*/)
{
    ParseEntry("unc", current.uncertainties);
    ParseEntry("sf", current.scale_factors);
    ParseEntry("ref_category", current.ref_category);
}

} // namespace analysis
