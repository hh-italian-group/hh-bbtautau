#pragma once
#include <iostream>
#include <algorithm>
#include <boost/format.hpp>
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "Analysis/include/SampleDescriptor.h"
#include "Analysis/include/AnaTuple.h"
#include "Analysis/include/EventAnalyzerDataId.h"

namespace analysis{

template<typename _FirstLeg, typename _SecondLeg>
class DYModel { // simple analyzer definition
public:
    using FirstLeg = _FirstLeg;
    using SecondLeg = _SecondLeg;
    using EventInfo = ::analysis::EventInfo<FirstLeg, SecondLeg>;
    DYModel(const SampleDescriptor& sample)
    {
        // Analyzer initialization (e.g. open input/output files, parse configs...)
        std::cout<<"Name = "<<sample.name<<std::endl;
        std::cout<<"Name suffix = "<<sample.name_suffix<<std::endl;

        const auto& param_names = sample.GetModelParameterNames();
        const auto b_param_iter = param_names.find("b");
        if(b_param_iter == param_names.end())
            throw exception("Unable to find b_parton WP for DY sample");
        b_index = b_param_iter->second;

        std::string ht_suffix = "ht";
        ht_found = sample.name_suffix.find(ht_suffix) != std::string::npos;

        fit_method = sample.fit_method;

        const auto ht_param_iter = param_names.find("ht");
        if(ht_param_iter == param_names.end())
        throw exception("Unable to find ht WP for DY smaple");
        ht_index = ht_param_iter->second;

        for(const auto& sample_wp : sample.working_points) {
            std::cout<<" Working point names = "<<sample_wp.full_name<<std::endl;
            const size_t n_b_partons = static_cast<size_t>(sample_wp.param_values.at(b_index));
            const size_t ht_wp = static_cast<size_t>(sample_wp.param_values.at(ht_index));
            working_points_map[std::pair<size_t,size_t>(n_b_partons,ht_wp)] = sample_wp;
            if(ht_found || fit_method == DYFitModel::NbjetBins_htBins) ht_wp_set.insert(ht_wp);
        }
        if(fit_method == DYFitModel::NbjetBins){
            auto input_file = root_ext::OpenRootFile(sample.norm_sf_file);
            auto scale_factor_histo =  std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*input_file,ToString(fit_method)
                                                                                    +"/scale_factors"));
            int nbins = scale_factor_histo->GetNbinsX();
            for(int i=1; i<=nbins;i++){
                std::string scale_factor_name = scale_factor_histo->GetXaxis()->GetBinLabel(i);
                double value = scale_factor_histo->GetBinContent(i);
                scale_factor_name.erase(2,3);
                scale_factor_maps[scale_factor_name] = value;
            }
        }
        else if(fit_method != DYFitModel::None)
            throw exception("Unable to find the fit method");
    }

    void ProcessEvent(const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight,
                 bbtautau::AnaTupleWriter::DataIdMap& dataIds){
        //static constexpr double pt_cut =18, b_Flavour = 5;
        /*size_t n_genJets = 0;
        for(const auto& b_Candidates : event.GetHiggsBB().GetDaughterMomentums()) {
            for(size_t i=0; i<event->genJets_p4.size(); i++){
                const auto& jet_p4 = event->genJets_p4.at(i);
                const auto& jet_hadronFlavour = event->genJets_hadronFlavour.at(i);
                double deltaR = ROOT::Math::VectorUtil::DeltaR(b_Candidates, jet_p4);
                if (jet_p4.Pt() <= pt_cut || jet_hadronFlavour != b_Flavour || deltaR >= 0.3) continue;
                n_genJets++;
            }
        }*/
        //unsigned int n_bJets = event->lhe_n_b_partons;
        unsigned int n_bJets = event->jets_nTotal_hadronFlavour_b;
        double lheHT = event->lhe_HT;
        int ht_wp = 0;
        std::pair<size_t,size_t> p;
        if(ht_found || fit_method == DYFitModel::NbjetBins_htBins)
            ht_wp = GetHTWP(lheHT);
        p = std::make_pair(std::min((unsigned int)2, n_bJets),0);
        std::map<std::pair<size_t,size_t>,SampleDescriptorBase::Point>::iterator it = working_points_map.find(p);
        if(it == working_points_map.end())
            throw exception("Unable to find WP for DY event with  jets_nTotal_hadronFlavour_b = %1%") % event->jets_nTotal_hadronFlavour_b;
           // throw exception("Unable to find WP for DY event with lhe_n_b_partons = %1%") % event->lhe_n_b_partons;

        auto sample_wp = it->second;
        const auto finalId = anaDataId.Set(sample_wp.full_name);
        double norm_sf = 1;
        if(fit_method == DYFitModel::NbjetBins)
            norm_sf = scale_factor_maps.at(sample_wp.full_name);
        else if(fit_method == DYFitModel::NbjetBins_htBins){
            if(ht_found) norm_sf = scale_factor_maps.at(sample_wp.full_name);
            else norm_sf = scale_factor_maps.at(sample_wp.full_name+ToString(ht_wp)+"ht");
        }
        dataIds[finalId] = std::make_tuple(weight * norm_sf, event.GetMvaScore());

    }

    int GetHTWP(double ht)
    {
        auto prev = ht_wp_set.begin();
        for(auto iter = std::next(prev); iter != ht_wp_set.end() && *iter < ht; ++iter) {
            prev = iter;
        }
        return static_cast<int>(*prev);
    }

private:
    std::map<std::string,double> scale_factor_maps;
    size_t b_index;
    size_t ht_index;
    std::map<std::pair<size_t,size_t>,SampleDescriptorBase::Point> working_points_map;
    DYFitModel fit_method;
    bool ht_found;
    std::set<size_t> ht_wp_set;
};
}
