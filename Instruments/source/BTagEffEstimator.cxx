/*! Estimate btag efficiencies.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include <vector>
#include <set>
#include <map>

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Core/include/EventTuple.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Cuts/include/btag_Run2.h"
#include "h-tautau/Cuts/include/hh_bbtautau_Run2.h"
#include "AnalysisTools/Core/include/Tools.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "h-tautau/Core/include/TauIdResults.h"
#include "h-tautau/JetTools/include/BTagger.h"
#include "h-tautau/Analysis/include/EventInfo.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, output_file);
    REQ_ARG(analysis::SignalMode, mode);
    REQ_ARG(analysis::BTaggerKind, btagger);
    REQ_ARG(analysis::Period, period);
    OPT_ARG(unsigned, n_threads, 1);
    REQ_ARG(std::vector<std::string>, input_file);
};

namespace analysis {

class BTagData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    const std::vector<double> x_bins = {20,30,40,60,100,150,200,300,1000};
    const std::vector<double> y_bins = {0,0.6,1.2,2.1,2.4};
    TH2D_ENTRY_CUSTOM(h2,x_bins,y_bins)
    TH2D_ENTRY_CUSTOM(eff,x_bins,y_bins)
    TH2D_ENTRY_CUSTOM(validate,x_bins,y_bins)
    TH2D_ENTRY_CUSTOM(val_check,x_bins,y_bins)
};

class BTagEffEstimator {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    BTagEffEstimator(const Arguments& _args) :
        args(_args), outfile(root_ext::CreateRootFile(args.output_file())), signalObjectSelector(args.mode())
    {
        ROOT::EnableThreadSafety();
        if(args.n_threads() > 1)
            ROOT::EnableImplicitMT(args.n_threads());

        anaDataMap["All"+eff] = std::make_shared<BTagData>(outfile,tools::FullPath({"All",effMap[eff]}));
        anaDataMap["All"+valCha] = std::make_shared<BTagData>(outfile,tools::FullPath({"All",valMap[valCha]}));
        for(const auto& tau_sign : tau_signs){
            for (const auto& tau_iso : tau_isos){
                anaDataMap[tau_sign+tau_iso+eff] = std::make_shared<BTagData>(
                            outfile,tools::FullPath({tau_sign, tau_iso, effMap[eff]}));
                anaDataMap[tau_sign+tau_iso+valCha] = std::make_shared<BTagData>(
                            outfile,tools::FullPath({tau_sign,tau_iso,valMap[valCha]}));
            }
            anaDataMap[tau_sign+valIso] = std::make_shared<BTagData>(outfile,tau_sign+"/"+valMap[valIso]);
        }

        for (const auto& tau_iso : tau_isos)
            anaDataMap[valQ+tau_iso] = std::make_shared<BTagData>(outfile,valMap[valQ]+"/"+tau_iso);

    }
    void Run()
    {
        static const std::set<std::string> channels = { ToString(Channel::ETau), ToString(Channel::MuTau),
            ToString(Channel::TauTau)};//, ToString(Channel::MuMu) };
        std::string channel_all = "all";
        std::set<std::string> channel_names = channels;
        channel_names.insert(channel_all);
        const BTagger bTagger(args.period(),args.btagger());
        btag_working_points = {
            {"L", DiscriminatorWP::Loose }, {"M", DiscriminatorWP::Medium}, {"T", DiscriminatorWP::Tight}
        };
        static const std::string btag_wp_all = "all";
        static const std::map<int, std::string> flavours = { { 5, "b" }, { 4, "c" }, { 0, "udsg" } };
        static const std::string flavour_all = "all";
        std::set<std::string> flavour_names = analysis::tools::collect_map_values(flavours);
        flavour_names.insert(flavour_all);
        static const std::string num = "Num", denom = "Denom";

        /*static const std::set<std::string> disabled_branches;
        static const std::set<std::string> enabled_branches = {
            "jets_p4", "SVfit_p4", "extramuon_veto", "extraelec_veto", "q_1", "q_2", "tauId_keys_1", "tauId_values_1",
            "tauId_keys_2", "tauId_values_2", "jets_mva", "jets_csv", "jets_deepCsv_BvsAll", "jets_hadronFlavour",
            "jets_pu_id", "jets_deepFlavour_b", "jets_deepFlavour_bb", "jets_deepFlavour_lepb"
        };*/


        //bool apply_pu_id_cut = args.apply_pu_id_cut() != "no";
        //DiscriminatorWP pu_wp = DiscriminatorWP::Medium;
        //if(apply_pu_id_cut && args.period()=="2016") pu_wp = analysis::Parse<DiscriminatorWP>(args.apply_pu_id_cut());

        for(const auto& channel : channels) {
            const auto leg_types = GetChannelLegTypes(analysis::Parse<Channel>(channel));
            for (const auto& name : args.input_file()){
                std::shared_ptr<TFile> in_file(root_ext::OpenRootFile(name));
                std::shared_ptr<EventTuple> tuple;
                try {
                    tuple = std::make_shared<EventTuple>(channel, in_file.get(), true);//, disabled_branches,
                                                         //enabled_branches);
                } catch(std::exception&) {
                    std::cerr << "WARNING: tree "<<channel<<" not found in file"<<name<< std::endl;
                    continue;
                }

                std::cout << "Processing " << name << "/" << channel << std::endl;

                for(const Event& event : *tuple){
                    auto eventInfo = EventInfo::Create(event, signalObjectSelector, bTagger, DiscriminatorWP::Medium);
                    if(!eventInfo) continue;
                    // const EventEnergyScale es = static_cast<EventEnergyScale>(event.eventEnergyScale);
                   /* if (args.period() == Period::Run2016 && (es != EventEnergyScale::Central || event.jets_p4.size() < 2 || event.extraelec_veto
                            || event.extramuon_veto
                            || std::abs(eventInfo->GetHiggsBB().GetFirstDaughter().GetMomentum().eta()) >= cuts::btag_2016::eta
                            || std::abs(eventInfo->GetHiggsBB().GetSecondDaughter().GetMomentum().eta()) >= cuts::btag_2016::eta)) continue;

                    else if (args.period() == Period::Run2017 && (es != EventEnergyScale::Central || event.jets_p4.size() < 2 || event.extraelec_veto
                            || event.extramuon_veto
                            || std::abs(eventInfo->GetHiggsBB().GetFirstDaughter().GetMomentum().eta()) >= cuts::btag_2017::eta
                            || std::abs(eventInfo->GetHiggsBB().GetSecondDaughter().GetMomentum().eta()) >= cuts::btag_2017::eta)) continue;
                    */
                    if(!eventInfo->HasBjetPair()) continue;

                    /*auto bb = eventInfo->GetHiggsBB().GetFirstDaughter().GetMomentum() + eventInfo->GetHiggsBB().GetSecondDaughter().GetMomentum();

                    if (!cuts::hh_bbtautau_2017::hh_tag::IsInsideMassWindow(eventInfo->GetSVFitResults().momentum.mass(),bb.mass())) continue;*/

                    std::string tau_sign = (eventInfo->GetLeg(1)->charge()+eventInfo->GetLeg(2)->charge()) == 0 ? "OS" : "SS";

                    const bool passTauId = (leg_types.first != LegType::tau || eventInfo->GetLeg(1)->Passed(analysis::TauIdDiscriminator::byDeepTau2017v2p1VSjet,DiscriminatorWP::Medium))
                            && (leg_types.second != LegType::tau || eventInfo->GetLeg(2)->Passed(analysis::TauIdDiscriminator::byDeepTau2017v2p1VSjet,DiscriminatorWP::Medium));


                    std::string tau_iso = passTauId ? "Iso" : "NonIso";

                    for (const auto& jet_candidate : eventInfo->GetCentralJets()){
                        const auto& jet = jet_candidate->GetMomentum();

                        //PU correction
                        /*if(apply_pu_id_cut){
                            if(args.period()=="2016"){
                                double jet_mva = event.jets_mva.at(i);
                                if(!PassJetPuId(jet.Pt(),jet_mva,pu_wp)) continue;
                            }
                                if((event.jets_pu_id.at(i) & (1 << 2)) == 0) continue;
                        }*/

                        const int jet_hadronFlavour = (*jet_candidate)->hadronFlavour();
                        const std::string& jet_flavour = flavours.at(jet_hadronFlavour);

                        //For folder/subfolder structure in Sign and Isolation
                        anaDataMap[tau_sign+tau_iso+eff]->h2(denom, flavour_all, btag_wp_all, channel).
                            Fill(jet.Pt(), std::abs(jet.Eta()));
                        if (channel != "muMu") anaDataMap[tau_sign+tau_iso+eff]->h2(denom, flavour_all, btag_wp_all,
                                channel_all).Fill(jet.Pt(), std::abs(jet.Eta()));
                        anaDataMap[tau_sign+tau_iso+eff]->h2(denom, jet_flavour, btag_wp_all, channel).
                            Fill(jet.Pt(), std::abs(jet.Eta()));
                        if (channel != "muMu") anaDataMap[tau_sign+tau_iso+eff]->h2(denom, jet_flavour, btag_wp_all,
                                channel_all).Fill(jet.Pt(), std::abs(jet.Eta()));

                        //For "ALL" folder
                        anaDataMap["All"+eff]->h2(denom, flavour_all, btag_wp_all, channel).
                            Fill(jet.Pt(), std::abs(jet.Eta()));
                        if (channel != "muMu") anaDataMap["All"+eff]->h2(denom, flavour_all, btag_wp_all,
                                channel_all).Fill(jet.Pt(), std::abs(jet.Eta()));
                        anaDataMap["All"+eff]->h2(denom, jet_flavour, btag_wp_all, channel).
                            Fill(jet.Pt(), std::abs(jet.Eta()));
                        if (channel != "muMu") anaDataMap["All"+eff]->h2(denom, jet_flavour, btag_wp_all,
                                channel_all).Fill(jet.Pt(), std::abs(jet.Eta()));

                        for(const auto& btag_wp : btag_working_points) {
                            if(bTagger.Pass(**jet_candidate, btag_wp.second)) {
                                //For folder/subfolder structure in Sign and Isolation
                                anaDataMap[tau_sign+tau_iso+eff]->h2(num, flavour_all, btag_wp.first, channel).
                                    Fill(jet.Pt(), std::abs(jet.Eta()));
                                if (channel != "muMu") anaDataMap[tau_sign+tau_iso+eff]->h2(num, flavour_all,
                                        btag_wp.first, channel_all).Fill(jet.Pt(), std::abs(jet.Eta()));
                                anaDataMap[tau_sign+tau_iso+eff]->h2(num, jet_flavour,btag_wp.first, channel).
                                    Fill(jet.Pt(), std::abs(jet.Eta()));
                                if (channel != "muMu") anaDataMap[tau_sign+tau_iso+eff]->h2(num, jet_flavour,
                                        btag_wp.first, channel_all).Fill(jet.Pt(), std::abs(jet.Eta()));

                                //For "ALL" folder
                                anaDataMap["All"+eff]->h2(num, flavour_all, btag_wp.first, channel).
                                    Fill(jet.Pt(), std::abs(jet.Eta()));
                                if (channel != "muMu") anaDataMap["All"+eff]->h2(num, flavour_all,
                                        btag_wp.first, channel_all).Fill(jet.Pt(), std::abs(jet.Eta()));
                                anaDataMap["All"+eff]->h2(num, jet_flavour,btag_wp.first, channel).
                                    Fill(jet.Pt(), std::abs(jet.Eta()));
                                if (channel != "muMu") anaDataMap["All"+eff]->h2(num, jet_flavour,
                                        btag_wp.first, channel_all).Fill(jet.Pt(), std::abs(jet.Eta()));

                            }
                        }
                    }//end loop on jets
                }//end loop on events
            }// end loop on files
        }//end loop on channel


        //Create Efficiency Histograms
        for(const auto& channel : channel_names){
            for(const auto& btag_wp : btag_working_points) {
                for(const auto& jet_flavour : flavour_names){
                    for(const auto& tau_sign : tau_signs){
                        for (const auto& tau_iso : tau_isos){
                            //For folder/subfolder structure in Sign and Isolation
                            anaDataMap[tau_sign+tau_iso+eff]->eff(jet_flavour,btag_wp.first,channel).
                                CopyContent(anaDataMap[tau_sign+tau_iso+eff]->h2(num,jet_flavour,btag_wp.first,channel));
                            anaDataMap[tau_sign+tau_iso+eff]->eff(jet_flavour,btag_wp.first,channel).
                                Divide(&anaDataMap[tau_sign+tau_iso+eff]->h2(denom,jet_flavour, btag_wp_all,channel));
                        }
                    }

                    //For "ALL" folder
                    anaDataMap["All"+eff]->eff(jet_flavour,btag_wp.first,channel).
                        CopyContent(anaDataMap["All"+eff]->h2(num,jet_flavour,btag_wp.first,channel));
                    anaDataMap["All"+eff]->eff(jet_flavour,btag_wp.first,channel).
                        Divide(&anaDataMap["All"+eff]->h2(denom,jet_flavour, btag_wp_all,channel));

                }
            }
        }//End of loop to create Efficiency Histograms


        //Make Validation plots

        for(auto channel1 = channel_names.begin(); channel1 != channel_names.end(); ++ channel1){
            for(const auto& btag_wp : btag_working_points) {
                for(const auto& jet_flavour : flavour_names){
                    for(auto tau_sign1 = tau_signs.begin(); tau_sign1 != tau_signs.end(); ++tau_sign1 ){
                        for (auto tau_iso1 = tau_isos.begin(); tau_iso1 != tau_isos.end(); ++tau_iso1){

                            //Make Validation plots between tau signs
                            for(auto tau_sign2 = std::next(tau_sign1); tau_sign2 != tau_signs.end(); ++tau_sign2){

                                anaDataMap[valQ+*tau_iso1]->validate(jet_flavour,btag_wp.first,*channel1,*tau_sign2,
                                        *tau_sign1).CopyContent(anaDataMap[*tau_sign2+*tau_iso1+eff]->
                                        eff(jet_flavour, btag_wp.first,*channel1));
                                anaDataMap[valQ+*tau_iso1]->validate(jet_flavour,btag_wp.first,*channel1,*tau_sign2,
                                        *tau_sign1).Add(&anaDataMap[*tau_sign1+*tau_iso1+eff]->
                                        eff(jet_flavour, btag_wp.first,*channel1),-1.0);
                                checkValidation(anaDataMap[valQ+*tau_iso1]->validate(jet_flavour,btag_wp.first,
                                        *channel1,*tau_sign2,*tau_sign1),tools::FullPath({valMap[valQ],*tau_iso1}),
                                        anaDataMap[valQ+*tau_iso1]->val_check(jet_flavour,btag_wp.first,*channel1,
                                            *tau_sign2,*tau_sign1),
                                        anaDataMap[*tau_sign1+*tau_iso1+eff]->eff(jet_flavour,
                                                                                  btag_wp.first,*channel1),
                                        anaDataMap[*tau_sign2+*tau_iso1+eff]->eff(jet_flavour,
                                                                                  btag_wp.first,*channel1));
                            }

                            //Make Vlaidation plots between Isolations
                            for (auto tau_iso2 = std::next(tau_iso1); tau_iso2 != tau_isos.end(); ++tau_iso2){

                                anaDataMap[*tau_sign1+valIso]->validate(jet_flavour,btag_wp.first,*channel1,*tau_iso2,
                                        *tau_iso1).CopyContent(anaDataMap[*tau_sign1+*tau_iso2+eff]->
                                        eff(jet_flavour, btag_wp.first,*channel1));
                                anaDataMap[*tau_sign1+valIso]->validate(jet_flavour,btag_wp.first,*channel1,*tau_iso2,
                                        *tau_iso1).Add(&anaDataMap[*tau_sign1+*tau_iso1+eff]->
                                        eff(jet_flavour, btag_wp.first,*channel1),-1.0);
                                checkValidation(anaDataMap[*tau_sign1+valIso]->validate(jet_flavour,btag_wp.first,
                                            *channel1,*tau_iso2,*tau_iso1),tools::FullPath({*tau_sign1,valMap[valIso]}),
                                        anaDataMap[*tau_sign1+valIso]->val_check(jet_flavour,btag_wp.first,
                                            *channel1,*tau_iso2,*tau_iso1),
                                        anaDataMap[*tau_sign1+*tau_iso1+eff]->eff(jet_flavour,
                                                                                  btag_wp.first,*channel1),
                                        anaDataMap[*tau_sign1+*tau_iso2+eff]->eff(jet_flavour,
                                                                                  btag_wp.first,*channel1));
                            }

                            //Make Validation plots between Channels In the OS or SS folders
                            for (auto channel2 = std::next(channel1); channel2 != channel_names.end(); ++channel2){

                                anaDataMap[*tau_sign1+*tau_iso1+valCha]->validate(jet_flavour,btag_wp.first,*channel2,
                                        *channel1).CopyContent(anaDataMap[*tau_sign1+*tau_iso1+eff]->
                                        eff(jet_flavour, btag_wp.first,*channel2));
                                anaDataMap[*tau_sign1+*tau_iso1+valCha]->validate(jet_flavour,btag_wp.first,*channel2,
                                        *channel1).Add(&anaDataMap[*tau_sign1+*tau_iso1+eff]->
                                        eff(jet_flavour, btag_wp.first,*channel1),-1.0);
                                checkValidation(anaDataMap[*tau_sign1+*tau_iso1+valCha]->validate(jet_flavour,
                                        btag_wp.first,*channel2,*channel1),tools::FullPath({*tau_sign1,*tau_iso1,
                                        valMap[valCha]}),anaDataMap[*tau_sign1+*tau_iso1+valCha]->val_check(jet_flavour,
                                        btag_wp.first,*channel2,*channel1),
                                        anaDataMap[*tau_sign1+*tau_iso1+eff]->eff(jet_flavour, btag_wp.first,
                                                                                  *channel1),
                                        anaDataMap[*tau_sign1+*tau_iso1+eff]->eff(jet_flavour, btag_wp.first,
                                                                                  *channel2));
                            }
                        }
                    }
                    //Make Validation plot between Channel in "All" folder
                    for (auto channel2 = std::next(channel1); channel2 != channel_names.end(); ++channel2){

                        anaDataMap["All"+valCha]->validate(jet_flavour,btag_wp.first,*channel2,*channel1).
                                CopyContent(anaDataMap["All"+eff]->eff(jet_flavour, btag_wp.first,*channel2));
                        anaDataMap["All"+valCha]->validate(jet_flavour,btag_wp.first,*channel2,*channel1).
                                Add(&anaDataMap["All"+eff]->eff(jet_flavour, btag_wp.first,*channel1),-1.0);
                        checkValidation(
                                anaDataMap["All"+valCha]->validate(jet_flavour,btag_wp.first,*channel2,*channel1),
                                tools::FullPath({"All",valMap[valCha]}),
                                anaDataMap["All"+valCha]->val_check(jet_flavour,btag_wp.first,*channel2,*channel1),
                                anaDataMap["All"+eff]->eff(jet_flavour, btag_wp.first,*channel1),
                                anaDataMap["All"+eff]->eff(jet_flavour, btag_wp.first,*channel2));
                    }
                }
            }
        }//End of creating the validation plots

    }// End of Run Function
private:
    Arguments args;
    std::shared_ptr<TFile> outfile;
    //BTagData anaData;
    std::map<std::string,std::shared_ptr<BTagData>> anaDataMap;
    std::vector<std::string> tau_signs = {"SS","OS"};
    std::vector<std::string> tau_isos = {"NonIso","Iso"};
    std::string eff = "Eff";
    std::map<std::string,std::string> effMap ={ {"Eff", "Efficiency"} };
    std::string valCha = "valCha", valIso = "valIso", valQ = "valQ";
    std::map<std::string,std::string> valMap = { {"valCha","ValidationChannel"} , {"valIso","ValidationIsolation"},
                                                 {"valQ","ValidationCharge"} };
    std::map<std::string, DiscriminatorWP> btag_working_points;
    SignalObjectSelector signalObjectSelector;

    // static bool PassTauIdCut(TauIdResults::BitsContainer id_bits)
    // {
    //     static const std::string disc_name = "byMediumIsolationMVArun2v1DBoldDMwLT";
    //     static const size_t disc_index = TauIdResults::GetBitRefsByName().at(disc_name);
    //
    //     const TauIdResults id_results(id_bits);
    //     return id_results.Result(disc_index);
    // }

    static bool PassJetPuId(double pt, double mva, DiscriminatorWP wp)
    {
        //PU Id cuts
        //https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/RecoJets/JetProducers/python/PileupJetIDCutParams_cfi.py#L31
        static const std::map<DiscriminatorWP, std::vector<std::pair<double,double> > > puId_working_points =
        { { DiscriminatorWP::Tight, { {30,0.26}, {50,0.62} } }, { DiscriminatorWP::Medium, { {30,-0.49}, {50,-0.06} } },
            { DiscriminatorWP::Loose, { {30,-0.96}, {50,-0.92} } } };

        if(!puId_working_points.count(wp)) throw analysis::exception("Unknown working point '%1%'.") % wp;
        const auto& working_point = puId_working_points.at(wp);
        for ( const auto& cut_values: working_point){
            if (pt< cut_values.first) return mva > cut_values.second;
        }
        return true;
    }

    static void checkValidation(const TH2D& validation_histo, const std::string histo_path,TH2D& output_histo,
                                const TH2D& eff_histo1, const TH2D& eff_histo2){
        int nbinsX = validation_histo.GetNbinsX();
        int nbinsY = validation_histo.GetNbinsY();
        std::string histo_name = tools::FullPath({histo_path,validation_histo.GetName()});
        for(int i=1;i<= nbinsX;i++){
            for(int j=1;j<=nbinsY;j++){
                double xPosition = validation_histo.GetXaxis()->GetBinCenter(i);
                double yPosition = validation_histo.GetYaxis()->GetBinCenter(j);
                double value = validation_histo.GetBinContent(i,j);
                double error = validation_histo.GetBinError(i,j);
                double n_sigma = error != 0 ? std::floor(value/error) : 0;
                double eff1 = eff_histo1.GetBinContent(i,j);
                double eff2 = eff_histo2.GetBinContent(i,j);
                bool has_entry = (eff1 != 0) && (eff2 != 0);
                n_sigma = has_entry ? n_sigma : 0 ;
                if (n_sigma < 0) n_sigma = n_sigma+1;
                output_histo.SetBinContent(i,j,n_sigma);
                if(n_sigma !=0) std::cout<<n_sigma<<" sigma: "<<histo_name<<" in (pt,eta) bin ("<<
                    xPosition<<","<<yPosition<<"):"<<value<<" +- "<<error<<std::endl;
            }
        }
    }

};

} // namespace analysis

PROGRAM_MAIN(analysis::BTagEffEstimator, Arguments) // definition of the main program function
