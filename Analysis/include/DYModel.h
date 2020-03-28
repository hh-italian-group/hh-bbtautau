/*! Definition of DYModel class, the class for DY estimation.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "h-tautau/Analysis/include/EventInfo.h"
#include "hh-bbtautau/Analysis/include/SampleDescriptor.h"
#include "hh-bbtautau/Analysis/include/AnaTuple.h"
#include "hh-bbtautau/Analysis/include/EventAnalyzerDataId.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"

namespace analysis{

class DYModelBase {
public:
    virtual ~DYModelBase(){} //destructor
    virtual void ProcessEvent(const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight,
                      bbtautau::AnaTupleWriter::DataIdMap& dataIds) = 0;

};

class DYModel : public DYModelBase {
public:
    DYModel(const SampleDescriptor& sample,const std::string& working_path);

    virtual void ProcessEvent(const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight,
                      bbtautau::AnaTupleWriter::DataIdMap& dataIds) override;

    template<typename T>
    static size_t Get2WP(T value, std::set<size_t>& wp_set)
    {
        auto prev = wp_set.begin();
        for(auto iter = std::next(prev); iter != wp_set.end() && *iter < value; ++iter) {
            prev = iter;
        }
        return (*prev);
    }

private:
    std::map<std::string,double> scale_factor_maps;
    size_t b_index;
    size_t ht_index;
    size_t njet_index;
    size_t pt_index;
    std::map<std::pair<size_t,size_t>, SampleDescriptorBase::Point> working_points_map;
    DYFitModel fit_method;
    bool ht_found;
    std::set<size_t> ht_wp_set;
    static const std::string& HT_suffix() { static const std::string s = "ht"; return s; }
    bool jet_found;
    std::set<size_t> njet_wp_set;
    static const std::string& NJet_suffix() { static const std::string s = "Jet"; return s; }
    bool pt_found;
    std::set<size_t> pt_wp_set;
    static const std::string& Pt_suffix() { static const std::string s = "Pt"; return s; }

    //static constexpr double b_Flavour = 5;
    //int b_Flavour=5;

    std::map<std::string, std::shared_ptr<TH1D>> pt_weight_histo_map;
    std::map<std::string, double>  fractional_weight_map;
    std::string sampleOrder;

};

class DYModel_HTT : public DYModelBase {
public:
    DYModel_HTT(const SampleDescriptor& sample,const std::string& working_path);

    virtual void ProcessEvent(const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight,
                      bbtautau::AnaTupleWriter::DataIdMap& dataIds) override;

private:
    std::shared_ptr<RooWorkspace> workspace;

};

}
