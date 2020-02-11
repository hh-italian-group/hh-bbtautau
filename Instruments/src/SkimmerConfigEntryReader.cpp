/*! Definition of classes to read TupleSkimmer configuration.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Instruments/include/SkimmerConfigEntryReader.h"

namespace analysis {
namespace tuple_skimmer {

void SetupEntryReader::EndEntry()
{
    CheckReadParamCounts("unc_sources", 1, Condition::less_equal);
    CheckReadParamCounts("channels", 1, Condition::less_equal);
    CheckReadParamCounts("period", 1, Condition::less_equal);
    CheckReadParamCounts("mode", 1, Condition::less_equal);
    CheckReadParamCounts("btag_wp", 1, Condition::less_equal);
    CheckReadParamCounts("use_cache", 1, Condition::less_equal);
    CheckReadParamCounts("common_weights", 1, Condition::less_equal);
    CheckReadParamCounts("n_splits", 1, Condition::less_equal);
    CheckReadParamCounts("split_seed", 1, Condition::less_equal);
    CheckReadParamCounts("jet_ordering", 1, Condition::less_equal);

    CheckReadParamCounts("apply_mass_cut", 1, Condition::less_equal);
    CheckReadParamCounts("apply_charge_cut", 1, Condition::less_equal);
    CheckReadParamCounts("apply_bb_cut", 1, Condition::less_equal);
    CheckReadParamCounts("apply_tau_iso", 1, Condition::less_equal);
    CheckReadParamCounts("keep_genJets", 1, Condition::less_equal);
    CheckReadParamCounts("keep_genParticles", 1, Condition::less_equal);
    CheckReadParamCounts("keep_MET_cov",1,Condition::less_equal);
    CheckReadParamCounts("tau_id_cuts", 1, Condition::less_equal);
    CheckReadParamCounts("massWindowParams", 0, Condition::greater_equal);
    CheckReadParamCounts("apply_kinfit", 1, Condition::less_equal);
    CheckReadParamCounts("applyTauId", 1, Condition::less_equal);

    ConfigEntryReaderT<Setup>::EndEntry();
}

void SetupEntryReader::ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                                     std::istringstream& /*ss*/)
{
    ParseEnumList("unc_sources", current.unc_sources);
    ParseEnumList("channels", current.channels);
    ParseEntry("period", current.period);
    ParseEntry("mode", current.mode);
    ParseEntry("btag_wp", current.btag_wp);
    ParseEntry("use_cache", current.use_cache);
    ParseEntryList("common_weights", current.common_weights);
    ParseEntry("n_splits", current.n_splits);
    ParseEntry("split_seed", current.split_seed);
    ParseEntry("jet_ordering", current.jet_ordering);

    ParseEntry("apply_mass_cut", current.apply_mass_cut);
    ParseEntry("apply_charge_cut", current.apply_charge_cut);
    ParseEntry("apply_bb_cut", current.apply_bb_cut);
    ParseEntry("apply_tau_iso", current.apply_tau_iso);
    ParseEntry("keep_genJets", current.keep_genJets);
    ParseEntry("keep_genParticles", current.keep_genParticles);
    ParseEntry("keep_MET_cov",current.keep_MET_cov);
    ParseEntryList("tau_id_cuts", current.tau_id_cuts);
    ParseEntry("massWindowParams", current.massWindowParams);
    ParseEntry("apply_kinfit",current.apply_kinfit);
    ParseEntry("applyTauId",current.applyTauId);
}

void SkimJobEntryReader::EndEntry()
{
    CheckReadParamCounts("merged_output", 1, Condition::less_equal);
    CheckReadParamCounts("apply_common_weights", 1, Condition::less_equal);
    CheckReadParamCounts("weights", 1, Condition::less_equal);
    CheckReadParamCounts("isData", 1, Condition::less_equal);

    const size_t n_files = GetReadParamCounts("file");
    const size_t n_files_ex = GetReadParamCounts("file_ex");
    const size_t n_files_xs = GetReadParamCounts("file_xs");

    if(n_files_xs && n_files + n_files_ex)
        throw exception("Cross-section should be specified for all files within a job or for non of them.");

    if(current.ProduceMergedOutput() && n_files_ex)
        throw exception("Entries merged_output and file_ex can't be both specified within single job.");

    ConfigEntryReaderT<SkimJob>::EndEntry();
}

void SkimJobEntryReader::ReadParameter(const std::string& param_name, const std::string& param_value,
                                       std::istringstream& /*ss*/)
{
    ParseEntry("merged_output", current.merged_output);
    ParseFileDescriptor(param_name, param_value);
    ParseEntry("apply_common_weights", current.apply_common_weights);
    ParseEntryList("weights", current.weights);
    ParseEntry("isData", current.isData);
}

void SkimJobEntryReader::ParseFileDescriptor(const std::string& param_name, const std::string& param_value)
{
    std::vector<std::string> inputs;
    if(param_name == "file") {
        inputs = SplitValueList(param_value, false);
        current.files.emplace_back(inputs);
    } else if(param_name == "file_ex") {
        auto columns = SplitValueList(param_value, true);
        if(columns.size() < 2)
            throw exception("Invalid extended file description.");
        inputs.insert(inputs.end(), columns.begin() + 1, columns.end());
        current.files.emplace_back(inputs, columns.at(0));
    } else if(param_name == "file_xs") {
        auto columns = SplitValueList(param_value, false);
        if(columns.size() < 2)
            throw exception("Invalid description for file with cross-section.");
        inputs.insert(inputs.end(), columns.begin() + 1, columns.end());
        const double xs = Parse<double>(columns.at(0));
        current.files.emplace_back(inputs, xs);
    } else {
        return;
    }
    throw param_parsed_exception();
}

} // namespace tuple_skimmer
} // namespace analysis
