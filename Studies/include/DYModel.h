#pragma once
#include <boost/format.hpp>
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "Analysis/include/SampleDescriptor.h"
#include "Analysis/include/AnaTuple.h"

namespace analysis{
class DYModel { // simple analyzer definition
public:
    using EventInfo = ::analysis::EventInfo<FirstLeg, SecondLeg>;

    std::map<std::string,double> scale_factor_maps;
    static size_t b_index;
    std::map<size_t,Point> working_points_map;
    DYModel(const SampleDescriptor& sample)
    {
        // Analyzer initialization (e.g. open input/output files, parse configs...)
        static auto const find_b_index = [&]() {
            const auto& param_names = sample.GetModelParameterNames();
            const auto b_param_iter = param_names.find("b");
            if(b_param_iter == param_names.end())
                throw exception("Unable to find b_parton WP for DY sample");
            return b_param_iter->second;
        };
        b_index = find_b_index();
        for(const auto& sample_wp : sample.working_points) {
            const size_t n_b_partons = static_cast<size_t>(sample_wp.param_values.at(b_index));
            working_points_map[n_b_partons] = sample_wp;
        }


        std::shared_ptr<TFile> input_file = root_ext::OpenRootFile("hh-bbtautau/McCorrections/data/DY_Scale_factors.root");
        std::shared_ptr<TH1F> scale_factor_histo = std::make_shared<TH1F>(root_ext::ReadObject<TH1F>(*input_file,
                                                                                                     "scale_factors"));
        int nbins = scale_factor_histo->GetNbinsX();
        for(int i=1; i<=nbins;i++){
            std::string scale_factor_name = scale_factor_histo->GetXaxis()->GetBinLabel(i);
            double value = scale_factor_histo->GetBinContent(i);
            scale_factor_name.erase(2,3);
            scale_factor_maps[scale_factor_name] = value;
        }
    }

    ProcessEvent(const EventAnalyzerDataId& anaDataId, EventInfo& event, double weight,
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
        double n_bJets = event->lhe_n_b_partons;
        std::map<size_t,Point>::iterator it = working_points_map.find(min(2,n_bJets));
        bool wp_found = it != working_points_map.end();
        if(wp_found){
            const auto& sample_wp = working_points_map.at(min(2,n_bjets));
            const auto finalId = anaDataId.Set(sample_wp.full_name);
            dataIds[finalId] = std::make_tuple(weight * sample_wp.norm_sf, event.GetMvaScore());
        }
        else
            throw exception("Unable to find WP for DY event with lhe_n_b_partons = %1%") % event->lhe_n_b_partons;

    }

private:

};
}
