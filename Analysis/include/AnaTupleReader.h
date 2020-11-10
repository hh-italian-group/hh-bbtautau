/*! Read AnaTuples.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <boost/bimap.hpp>
#include <ROOT/RDataFrame.hxx>
#include "EventAnalyzerDataId.h"
#include "hh-bbtautau/Analysis/include/EventTags.h"

namespace analysis {
namespace bbtautau {

class AnaTupleReader {
public:
    using DataId = EventAnalyzerDataId;
    using Hash = size_t;
    using DataIdBiMap = boost::bimap<DataId, Hash>;
    using DataIdMap = std::map<DataId, std::tuple<double, double>>;
    using NameSet = std::set<std::string>;
    using RDF = ROOT::RDF::RNode;

    AnaTupleReader(const std::string& file_name, Channel channel, NameSet& active_var_names,
                   const std::vector<std::string>& input_friends, const EventTagCreator& event_tagger, const std::string& mdnn_version);
    size_t GetNumberOfEntries() const;
    const DataId& GetDataIdByHash(Hash hash) const;
    const RDF& GetDataFrame() const;
    const std::list<RDF>& GetSkimmedDataFrames() const;
    const std::map<std::string, std::set<std::string>>& GetParametricVariables() const;
    boost::optional<std::string> TryGetBranchType(const std::string& branch_name) const;

private:
    void DefineBranches(const NameSet& active_var_names, bool all, const EventTagCreator& event_tagger, const std::string& mdnn_version);

    static std::vector<std::shared_ptr<TFile>> OpenFiles(const std::string& file_name,
                                                         const std::vector<std::string>& input_friends);
    static std::vector<std::shared_ptr<TTree>> ReadTrees(Channel channel,
                                                         const std::vector<std::shared_ptr<TFile>>& files);

private:
    std::vector<std::shared_ptr<TFile>> files;
    std::vector<std::shared_ptr<TTree>> trees;
    ROOT::RDataFrame dataFrame;
    RDF df;
    std::list<RDF> skimmed_df;
    DataIdBiMap known_data_ids;
    std::map<std::string, std::set<std::string>> parametric_vars;
    std::map<std::string, std::string> branch_types;
};

struct HyperPoint {
    boost::optional<int> spin;
    boost::optional<double> mass;
    boost::optional<double> kl;

    std::string ToString();
};

} // namespace bbtautau
} // namespace analysis
