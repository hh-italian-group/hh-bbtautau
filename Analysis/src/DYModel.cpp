/*! Definition of DYModel class, the class for DY estimation.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/DYModel.h"

namespace analysis{

DYModel::DYModel(const SampleDescriptor& sample,const std::string& working_path)
{
    sampleOrder = sample.sampleOrder;
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

    pt_found = sample.name_suffix.find(Pt_suffix()) != std::string::npos;
    if(pt_found){
        const auto pt_param_iter = param_names.find(Pt_suffix());
        if(pt_param_iter == param_names.end())
            throw exception("Unable to find Pt WP for DY smaple");
        pt_index = pt_param_iter->second;
    }

    for(const auto& sample_wp : sample.working_points) {
        const size_t n_b_partons = Parse<size_t>(sample_wp.param_values.at(b_index));
        if(ht_found){
            const size_t ht_wp = static_cast<size_t>(std::stod(sample_wp.param_values.at(ht_index)));
            working_points_map[std::pair<size_t,size_t>(n_b_partons,ht_wp)] = sample_wp;
            ht_wp_set.insert(ht_wp);
        }
        else if(jet_found){
            const size_t njet_wp = static_cast<size_t>(std::stod(sample_wp.param_values.at(njet_index)));
            working_points_map[std::pair<size_t,size_t>(n_b_partons,njet_wp)] = sample_wp;
            njet_wp_set.insert(njet_wp);
        }
        else if(pt_found){
            const size_t pt_wp = static_cast<size_t>(std::stod(sample_wp.param_values.at(pt_index)));
            working_points_map[std::pair<size_t,size_t>(n_b_partons,pt_wp)] = sample_wp;
            pt_wp_set.insert(pt_wp);
        }
        else working_points_map[std::pair<size_t,size_t>(n_b_partons,0)] = sample_wp;
    }
    if(fit_method == DYFitModel::NbjetBins || fit_method == DYFitModel::NbjetBins_htBins ||
            fit_method == DYFitModel::NbjetBins_NjetBins || fit_method == DYFitModel::NbjetBins_ptBins){
        auto input_file = root_ext::OpenRootFile(working_path+"/"+sample.norm_sf_file);
        auto scale_factor_histo =  std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*input_file,ToString(fit_method)
                                                                                +"/scale_factors"));
        int nbins = scale_factor_histo->GetNbinsX();
        for(int i=1; i<=nbins;i++){
            std::string scale_factor_name = scale_factor_histo->GetXaxis()->GetBinLabel(i);
            double value = scale_factor_histo->GetBinContent(i);
            if(scale_factor_name.find("DY") != std::string::npos){
                std::string sf_prefix = "SF_";
                if(scale_factor_name.find("DY_nlo") != std::string::npos) scale_factor_name.insert(7,sf_prefix);
                else if(scale_factor_name.find("DY_lo") != std::string::npos) scale_factor_name.insert(6,sf_prefix);
                scale_factor_maps[scale_factor_name] = value;
            }
        }
    }
    else if(fit_method != DYFitModel::None)
        throw exception("Unable to find the fit method");

    if(sampleOrder == "LO"){
        fractional_weight_map["0Jet"] = 0.813;
        fractional_weight_map["1Jet_0bJet"] = 1.100;
        fractional_weight_map["1Jet_1bJet"] = 1.316;
        fractional_weight_map["2Jet_0bJet"] = 1.044;
        fractional_weight_map["2Jet_1bJet"] = 1.232;
        fractional_weight_map["2Jet_2bJet"] = 1.222;

        auto NLO_weight_file = (root_ext::OpenRootFile(working_path+"/"+sample.NLO_weight_file));
        std::string histo_name = "h_ratio_pt";
        pt_weight_histo_map["0Jet"] = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*NLO_weight_file,histo_name+"0Jet"));
        pt_weight_histo_map["1Jet_0bJet"] = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*NLO_weight_file,histo_name+"1Jet_0bJet"));
        pt_weight_histo_map["1Jet_1bJet"] = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*NLO_weight_file,histo_name+"1Jet_1bJet"));
        pt_weight_histo_map["2Jet_0bJet"] = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*NLO_weight_file,histo_name+"2Jet_0bJet"));
        pt_weight_histo_map["2Jet_1bJet"] = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*NLO_weight_file,histo_name+"2Jet_1bJet"));
        pt_weight_histo_map["2Jet_2bJet"] = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*NLO_weight_file,histo_name+"2Jet_2bJet"));
    }

}

void DYModel::ProcessEvent(const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight,
                           bbtautau::AnaTupleWriter::DataIdMap& dataIds)
{
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

    auto n_selected_gen_jets = event->lhe_n_partons;
    size_t n_bJets = event->lhe_n_b_partons;

    /*if(sampleOrder == "LO"){
        std::string lhe_category = "";
        if(n_selected_gen_jets==0) lhe_category = "0Jet";
        else if (n_selected_gen_jets == 1){
            if(n_bJets==0) lhe_category= "1Jet_0bJet";
            else if (n_bJets==1) lhe_category="1Jet_1bJet";
        }
        else if (n_selected_gen_jets==2){
            if(n_bJets==0) lhe_category="2Jet_0bJet";
            else if(n_bJets==1) lhe_category="2Jet_1bJet";
            else if (n_bJets==2) lhe_category="2Jet_2bJet";
        }

        double fractional_weight = 0;
        double pt_weight =0;
        if (n_selected_gen_jets <= 2){
            fractional_weight = fractional_weight_map[lhe_category];

            for(size_t i=0;i<event->genParticles_p4.size();i++){
                if(event->genParticles_pdg.at(i) != 23) continue;
                double pt = event->genParticles_p4.at(i).Pt();
                pt_weight = pt_weight_histo_map[lhe_category]->GetBinContent(pt_weight_histo_map[lhe_category]->FindBin(pt));
            }
        }
        weight = weight*fractional_weight*pt_weight;
    }*/


    //unsigned int n_bJets = event->jets_nTotal_hadronFlavour_b;
    //double lheHT = event->lhe_HT;

    /*double gen_Zpt=-99.9;
    for(size_t i=0;i<event->genParticles_pdg.size();i++){
        float pdgId = event->genParticles_pdg.at(i);
        std::cout<<i<<" PdgId = "<<pdgId<<std::endl;
        if(pdgId != 23) continue;
        gen_Zpt = event->genParticles_p4.at(i).Pt();
    }
    double gen_weight_NLO = 0;
    double gen_weight_LO = 0;
    if(gen_Zpt != -99.9){
        gen_weight_NLO = (0.876979+gen_Zpt*(4.11598e-03)-(2.35520e-05)*gen_Zpt*gen_Zpt)*
                           (1.10211 * (0.958512 - 0.131835*TMath::Erf((gen_Zpt-14.1972)/10.1525)))*(gen_Zpt<140)
                           +0.891188*(gen_Zpt>=140);
        gen_weight_LO = (8.61313e-01+gen_Zpt*4.46807e-03-1.52324e-05*gen_Zpt*gen_Zpt)*
                                     (1.08683 * (0.95 - 0.0657370*TMath::Erf((gen_Zpt-11.)/5.51582)))*(gen_Zpt<140)
                                     +1.141996*(gen_Zpt>=140);
    }*/


    /*auto n_selected_gen_jets =  event->genJets_p4.size();
    size_t n_bJets =  static_cast<size_t>(std::count(event->genJets_hadronFlavour.begin(),
                                                     event->genJets_hadronFlavour.end(), b_Flavour));*/


    std::pair<size_t,size_t> p(std::min<size_t>(2,n_bJets),0);
    if(ht_found){
        double lheHT_otherjets = event->lhe_HT;
        size_t ht_wp = Get2WP(lheHT_otherjets,ht_wp_set);
        p.second = ht_wp;
    }
    else if(jet_found){
        size_t njet_wp = Get2WP(n_selected_gen_jets,njet_wp_set);
        p.second = njet_wp;
    }
    else if(pt_found){
        double gen_pt = 0.0;
        if(sampleOrder == "LO") gen_pt= 0;
        else if (sampleOrder == "NLO") gen_pt=25;
        bool Z_found =false;
        for(size_t i=0;i<event->genParticles_p4.size();i++){
            if(event->genParticles_pdg.at(i) != 23) continue;
            gen_pt = event->genParticles_p4.at(i).Pt();
            Z_found = true;
            break;
        }
        //int gen_match_muon = static_cast<int>(GenLeptonMatch::Muon);
        if(!Z_found){
            if(event.GetLeg(1)->gen_match() == GenLeptonMatch::Muon && event.GetLeg(2)->gen_match() == GenLeptonMatch::Muon) gen_pt = (event.GetLeg(1)->gen_p4() + event.GetLeg(2)->gen_p4()).Pt();
        }
        size_t pt_wp = Get2WP(gen_pt,pt_wp_set);
        p.second = pt_wp;
    }

    std::map<std::pair<size_t,size_t>,SampleDescriptorBase::Point>::iterator it = working_points_map.find(p);
    if(it == working_points_map.end())
        throw exception("Unable to find WP for DY event with  n_selected_bjet = %1%") % n_bJets;
        //throw exception("Unable to find WP for DY event with  jets_nTotal_hadronFlavour_b = %1%") % event->jets_nTotal_hadronFlavour_b;
       // throw exception("Unable to find WP for DY event with lhe_n_b_partons = %1%") % event->lhe_n_b_partons;

    auto sample_wp = it->second;
    const auto finalId = anaDataId.Set(sample_wp.full_name);
    double norm_sf = 1;

    if(event.HasBjetPair()){
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
        else if(fit_method == DYFitModel::NbjetBins_ptBins){
            if(pt_found) norm_sf = scale_factor_maps.at(sample_wp.full_name);
            else norm_sf = scale_factor_maps.at(sample_wp.full_name+"_"+ToString(it->first.second)+Pt_suffix());
        }
    }
    dataIds[finalId] = std::make_tuple(weight * norm_sf, event.GetMvaScore());
}

DYModel_HTT::DYModel_HTT(const SampleDescriptor& sample,const std::string& working_path) //constructor
{
    auto input_file = root_ext::OpenRootFile(working_path+"/"+sample.norm_sf_file);
    workspace =  std::shared_ptr<RooWorkspace>(root_ext::ReadObject<RooWorkspace>(*input_file,"w"));
}

void DYModel_HTT::ProcessEvent(const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight,
                           bbtautau::AnaTupleWriter::DataIdMap& dataIds)
{
    double norm_sf = 1;
    for(size_t i=0;i<event->genParticles_p4.size();i++){
        if(event->genParticles_pdg.at(i) != 23) continue;
        workspace->var("z_gen_mass")->setVal(event->genParticles_p4.at(i).M());
        workspace->var("z_gen_pt")->setVal(event->genParticles_p4.at(i).Pt());
        norm_sf = workspace->function("zptmass_weight_nom")->getVal();
        break;
    }

    dataIds[anaDataId] = std::make_tuple(weight * norm_sf, event.GetMvaScore());
}

}
