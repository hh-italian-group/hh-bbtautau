/*! Definition of DYModel class, the class for DY estimation.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

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

class DYModel {
public:
    DYModel(const SampleDescriptor& sample)
    {
        const auto& param_names = sample.GetModelParameterNames();
        const auto b_param_iter = param_names.find("b");
        if(b_param_iter == param_names.end())
            throw exception("Unable to find b_parton WP for DY sample");
        b_index = b_param_iter->second;

        ht_found = sample.name_suffix.find(HT_suffix()) != std::string::npos;

        fit_method = sample.fit_method;

        if(ht_found){
            const auto ht_param_iter = param_names.find(HT_suffix());
            if(ht_param_iter == param_names.end())
                throw exception("Unable to find ht WP for DY smaple");
            ht_index = ht_param_iter->second;
        }

        jet_found = sample.name_suffix.find(NJet_suffix()) != std::string::npos;
        if(jet_found){
            const auto njet_param_iter = param_names.find(NJet_suffix());
            if(njet_param_iter == param_names.end())
                throw exception("Unable to find njet WP for DY smaple");
            njet_index = njet_param_iter->second;
        }

        for(const auto& sample_wp : sample.working_points) {
            const size_t n_b_partons = static_cast<size_t>(sample_wp.param_values.at(b_index));
            if(ht_found){
                const size_t ht_wp = static_cast<size_t>(sample_wp.param_values.at(ht_index));
                working_points_map[std::pair<size_t,size_t>(n_b_partons,ht_wp)] = sample_wp;
                ht_wp_set.insert(ht_wp);
            }
            else if(jet_found){
                const size_t njet_wp = static_cast<size_t>(sample_wp.param_values.at(njet_index));
                working_points_map[std::pair<size_t,size_t>(n_b_partons,njet_wp)] = sample_wp;
                njet_wp_set.insert(njet_wp);
            }
            else working_points_map[std::pair<size_t,size_t>(n_b_partons,0)] = sample_wp;
        }
        if(fit_method == DYFitModel::NbjetBins || fit_method == DYFitModel::NbjetBins_htBins ||
                fit_method == DYFitModel::NbjetBins_NjetBins){
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

    void ProcessEvent(const EventAnalyzerDataId& anaDataId, EventInfoBase& event, double weight,
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
        //unsigned int n_bJets = event->jets_nTotal_hadronFlavour_b;
        //double lheHT = event->lhe_HT;

        auto n_selected_gen_jets =  event->genJets_p4.size();
        size_t n_bJets =  static_cast<size_t>(std::count(event->genJets_hadronFlavour.begin(),
                                                         event->genJets_hadronFlavour.end(), b_Flavour));
        std::pair<size_t,size_t> p(std::min<size_t>(2,n_bJets),0);
        if(ht_found){
            double lheHT_otherjets = event.CalculateGenHT(2);
            size_t ht_wp = Get2WP(lheHT_otherjets,ht_wp_set);
            p.second = ht_wp;
        }
        else if(jet_found){
            size_t njet_wp = Get2WP(n_selected_gen_jets,njet_wp_set);
            p.second = njet_wp;
        }

        std::map<std::pair<size_t,size_t>,SampleDescriptorBase::Point>::iterator it = working_points_map.find(p);
        if(it == working_points_map.end())
            throw exception("Unable to find WP for DY event with  n_selected_bjet = %1%") % n_bJets;
            //throw exception("Unable to find WP for DY event with  jets_nTotal_hadronFlavour_b = %1%") % event->jets_nTotal_hadronFlavour_b;
           // throw exception("Unable to find WP for DY event with lhe_n_b_partons = %1%") % event->lhe_n_b_partons;

        auto sample_wp = it->second;
        const auto finalId = anaDataId.Set(sample_wp.full_name);
        double norm_sf = 1;
        if(fit_method == DYFitModel::NbjetBins)
            norm_sf = scale_factor_maps.at(sample_wp.full_name);
        else if(fit_method == DYFitModel::NbjetBins_htBins){
            if(ht_found) norm_sf = scale_factor_maps.at(sample_wp.full_name);
            else norm_sf = scale_factor_maps.at(sample_wp.full_name+"_"+ToString(it->first.second)+HT_suffix());
        }
        else if(fit_method == DYFitModel::NbjetBins_NjetBins){
            if(jet_found) norm_sf = scale_factor_maps.at(sample_wp.full_name);
            else norm_sf = scale_factor_maps.at(sample_wp.full_name+"_"+ToString(it->first.second)+NJet_suffix());
        }
        dataIds[finalId] = std::make_tuple(weight * norm_sf, event.GetMvaScore());

    }

    size_t Get2WP(double value, std::set<size_t> wp_set) const
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
    std::map<std::pair<size_t,size_t>,SampleDescriptorBase::Point> working_points_map;
    DYFitModel fit_method;
    bool ht_found;
    std::set<size_t> ht_wp_set;
    static const std::string& HT_suffix() { static const std::string s = "ht"; return s; }
    bool jet_found;
    std::set<size_t> njet_wp_set;
    static const std::string& NJet_suffix() { static const std::string s = "Jet"; return s; }

    //static constexpr double b_Flavour = 5;
    int b_Flavour=5;
};
}
