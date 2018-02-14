/*! Definition of classes to read TupleSkimmer configuration.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "SkimmerConfig.h"

namespace analysis {
namespace tuple_skimmer {

class SetupEntryReader : public ConfigEntryReaderT<Setup> {
public:
    using ConfigEntryReaderT<Setup>::ConfigEntryReaderT;

    virtual void EndEntry() override
    {
        CheckReadParamCounts("energy_scales", 1, Condition::less_equal);
        CheckReadParamCounts("channels", 1, Condition::less_equal);
        CheckReadParamCounts("tau_ids", 1, Condition::less_equal);
        CheckReadParamCounts("period", 1, Condition::less_equal);
        CheckReadParamCounts("btag_wp", 1, Condition::less_equal);
        CheckReadParamCounts("common_weights", 1, Condition::less_equal);
        CheckReadParamCounts("n_splits", 1, Condition::less_equal);
        CheckReadParamCounts("split_seed", 1, Condition::less_equal);

        CheckReadParamCounts("apply_mass_cut", 1, Condition::less_equal);
        CheckReadParamCounts("apply_charge_cut", 1, Condition::less_equal);
        CheckReadParamCounts("keep_genJets", 1, Condition::less_equal);
        CheckReadParamCounts("keep_genParticles", 1, Condition::less_equal);
        CheckReadParamCounts("keep_MET_cov",1,Condition::less_equal);
        CheckReadParamCounts("tau_iso", 1, Condition::less_equal);
        CheckReadParamCounts("massWindowParams", 0, Condition::greater_equal);


        ConfigEntryReaderT<Setup>::EndEntry();
    }

    virtual void ReadParameter(const std::string& /*param_name*/, const std::string& /*param_value*/,
                               std::istringstream& /*ss*/) override
    {
        ParseEnumList("energy_scales", current.energy_scales);
        ParseEnumList("channels", current.channels);
        ParseEntryList("tau_ids", current.tau_ids);
        ParseEntry("period", current.period);
        ParseEntry("btag_wp", current.btag_wp);
        ParseEntryList("common_weights", current.common_weights);
        ParseEntry("n_splits", current.n_splits);
        ParseEntry("split_seed", current.split_seed);

        ParseEntry("apply_mass_cut", current.apply_mass_cut);
        ParseEntry("apply_charge_cut", current.apply_charge_cut);
        ParseEntry("keep_genJets", current.keep_genJets);
        ParseEntry("keep_genParticles", current.keep_genParticles);
        ParseEntry("keep_MET_cov",current.keep_MET_cov);
        ParseEntry("tau_id_cut", current.tau_id_cut);
        ParseEntry("massWindowParams", current.massWindowParams);
    }
};

class SkimJobEntryReader : public ConfigEntryReaderT<SkimJob> {
public:
    using ConfigEntryReaderT<SkimJob>::ConfigEntryReaderT;

    virtual void EndEntry() override
    {
        CheckReadParamCounts("merged_output", 1, Condition::less_equal);
        CheckReadParamCounts("apply_common_weights", 1, Condition::less_equal);
        CheckReadParamCounts("apply_dm_fix", 1, Condition::less_equal);
        CheckReadParamCounts("weights", 1, Condition::less_equal);

        const size_t n_files = GetReadParamCounts("file");
        const size_t n_files_ex = GetReadParamCounts("file_ex");
        const size_t n_files_xs = GetReadParamCounts("file_xs");

        if(n_files_xs && n_files + n_files_ex)
            throw exception("Cross-section should be specified for all files within a job or for non of them.");

        if(current.ProduceMergedOutput() && n_files_ex)
            throw exception("Entries merged_output and file_ex can't be both specified within single job.");

        ConfigEntryReaderT<SkimJob>::EndEntry();
    }

    virtual void ReadParameter(const std::string& param_name, const std::string& param_value,
                               std::istringstream& /*ss*/) override
    {
        ParseEntry("merged_output", current.merged_output);
        ParseFileDescriptor(param_name, param_value);
        ParseEntry("apply_common_weights", current.apply_common_weights);
        ParseEntry("apply_dm_fix", current.apply_dm_fix);
        ParseEntryList("weights", current.weights);
    }

private:
    void ParseFileDescriptor(const std::string& param_name, const std::string& param_value)
    {
        std::vector<std::string> inputs;
        if(param_name == "file") {
            inputs = SplitValueList(param_value, false);
            current.files.emplace_back(inputs);
        } else if(param_name == "file_ex") {
            auto columns = SplitValueList(param_value, false);
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
};

} // namespace tuple_skimmer
} // namespace analysis
