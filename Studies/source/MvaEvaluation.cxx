/*! Study for Mva Training
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "hh-bbtautau/Studies/include/MvaConfiguration.h"
#include "hh-bbtautau/Analysis/include/MvaVariables.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Ranking.h"
#include "TMVA/Config.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodBDT.h"
#include "TMVA/IMethod.h"
#include "TMVA/ClassifierFactory.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/ResultsClassification.h"
#include "TMVA/ROCCurve.h"
#include "TMath.h"
#include <fstream>
#include <random>
#include "hh-bbtautau/Analysis/include/MvaConfigurationReader.h"
#include "hh-bbtautau/Analysis/include/MvaReader.h"
#include "hh-bbtautau/Studies/include/MvaTuple.h"
#include "hh-bbtautau/Analysis/include/MvaVariablesStudy.h"
#include "hh-bbtautau/Studies/include/MvaMethods.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(std::string, file_xml);
    REQ_ARG(std::string, method_name);
    REQ_ARG(Long64_t, number_events);
    REQ_ARG(int, min);
    REQ_ARG(int, max);
    REQ_ARG(uint_fast32_t, number_sets);
    REQ_ARG(uint_fast32_t, seed);
    REQ_ARG(uint_fast32_t, seed2);
    REQ_ARG(bool, isLegacy);
    REQ_ARG(bool, isLow);
    OPT_ARG(std::vector<std::string>, variables, {});
};

namespace analysis {
namespace mva_study{

class MVAEvaluation {
public:
    using Event = ::ntuple::Event;
    using EventTuple = ::ntuple::EventTuple;
    using DataVector = std::vector<double>;
    using MassData = std::map<SampleId, DataVector>;


    MVAEvaluation(const Arguments& _args): args(_args),
        outfile(root_ext::CreateRootFile(args.output_file()+".root")), seed_split(0,29), test_vs_training(0, args.number_sets()-static_cast<uint_fast32_t>(1))
    {
        MvaSetupCollection setups;
        SampleEntryListCollection samples_list;

        ConfigReader configReader;
        MvaConfigReader setupReader(setups);
        configReader.AddEntryReader("SETUP", setupReader, true);
        SampleConfigReader sampleReader(samples_list);
        configReader.AddEntryReader("FILES", sampleReader, false);
        configReader.ReadConfig(args.cfg_file());
        if (args.isLegacy()) mva_setup = setups.at("LowF");
        else mva_setup = setups.at("Low");

        if(!samples_list.count("Samples"))
            throw exception("Samples don't found");
        samples = samples_list.at("Samples").files;
        std::cout<<"SAMPLES: "<<samples.size()<<std::endl;

        for (const auto& v : args.variables())
            enabled_vars.insert(v);
        std::cout<<"VARS: "<<enabled_vars.size()<<std::endl;

        reader = std::make_shared<MvaReader>();

    }

    std::map<SampleId, std::map<size_t,std::shared_ptr<TH1D>>> CreateOutputHistos(std::map<SampleId, std::map<size_t, std::vector<double>>> data)
    {
        std::map<SampleId, std::map<size_t,std::shared_ptr<TH1D>>> histo;
        for(const auto& sample : data){
            const auto& type = sample.second;
            for(const auto& entry : type){
                if (sample.first.IsSignal()){
                    for (const auto& value : entry.second){
                        if (!histo.count(sample.first))
                             histo[sample.first][tot] = std::make_shared<TH1D>(("Signal_"+std::to_string(sample.first.mass)+"_output").c_str(),("Signal"+std::to_string(sample.first.mass)+"_output").c_str(), nbin, bin_min, bin_max);
                        histo.at(sample.first).at(tot)->Fill(value);
                        if (!histo[sample.first].count(entry.first))
                            histo[sample.first][entry.first] = std::make_shared<TH1D>(("Signal_"+std::to_string(sample.first.mass)+"_output_type"+std::to_string(entry.first)).c_str(),("Signal_"+std::to_string(sample.first.mass)+"_output_type"+std::to_string(entry.first)).c_str(), nbin, bin_min, bin_max);
                        histo[sample.first][entry.first]->Fill(value);
                    }
                }
                if (sample.first.IsBackground()){
                    for (const auto& value : entry.second){
                        if (!histo.count(sample.first))
                            histo[sample.first][tot] = std::make_shared<TH1D>(("Bkg_"+std::to_string(sample.first.mass)+"_output").c_str(),("Bkg"+std::to_string(sample.first.mass)+"_output").c_str(), nbin, bin_min, bin_max);
                        histo.at(sample.first).at(tot)->Fill(value);
                        if (!histo[sample.first].count(entry.first))
                            histo[sample.first][entry.first] = std::make_shared<TH1D>(("Bkg_"+std::to_string(sample.first.mass)+"_output_type"+std::to_string(entry.first)).c_str(),("Bkg_"+std::to_string(sample.first.mass)+"_output_type"+std::to_string(entry.first)).c_str(), nbin, bin_min, bin_max);
                        histo[sample.first][entry.first]->Fill(value);
                    }
                }
            }
        }
        auto directory = root_ext::GetDirectory(*outfile.get(), "Evaluation");
        for(const auto& h : histo){
            for (auto his : h.second){
                root_ext::WriteObject(*his.second,directory);
            }
        }
        return histo;
    }

    std::vector<int> CreateMassRange(const int& min, const int& max) const
    {
        std::set<int> masses;
        for(const auto& sample : samples) {
            if(sample.id.IsSignal() && sample.id.mass>=min && sample.id.mass<=max)
                masses.insert(sample.id.mass);
        }
        return std::vector<int>(masses.begin(), masses.end());
    }

    void Run()
    {
        std::map<SampleId, std::map<size_t, std::vector<double>>> data;
        std::cout<<"Variabili iniziali: "<<enabled_vars.size()<<std::endl;
        const auto mass_range = CreateMassRange(args.min(), args.max());
        std::cout<<"SAMPLES nel RANGE:"<<mass_range.size() <<std::endl;
        ::analysis::Range<int> range(args.min(), args.max());
        auto vars = reader->AddRange(range, args.method_name(), args.file_xml(), enabled_vars, args.isLegacy(), args.isLow());
        std::mt19937_64 seed_gen(args.seed());

        for(size_t j = 0; j<mva_setup.channels.size(); j++){
            std::cout << mva_setup.channels[j] << std::endl;
            for(const SampleEntry& entry : samples)
            {
                if ( entry.id.IsSignal() && (entry.id.mass<args.min() || entry.id.mass>args.max())) continue;
                if ( entry.id.IsBackground() && Parse<Channel>(entry.channel) != mva_setup.channels.at(j) )
                    continue;
                auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
                EventTuple tuple(ToString(mva_setup.channels[j]), input_file.get(), true, {} , GetMvaBranches());
                Long64_t tot_entries = 0;
                std::mt19937_64 gen2(args.seed2());
                for(Long64_t current_entry = 0; tot_entries < args.number_events() && current_entry < tuple.GetEntries(); ++current_entry) {
                    uint_fast32_t seed1 = seed_split(gen2);
                    if (seed1>15) continue;
                    tuple.GetEntry(current_entry);
                    const Event& event = tuple.data();
                    if (static_cast<EventEnergyScale>(event.eventEnergyScale) != EventEnergyScale::Central || (event.q_1+event.q_2) != 0 || event.jets_p4.size() < 2
                        || event.extraelec_veto == true || event.extramuon_veto == true || event.jets_p4[0].eta() > cuts::btag_2016::eta
                        || event.jets_p4[1].eta() > cuts::btag_2016::eta)
                        continue;
                    auto bb = event.jets_p4[0] + event.jets_p4[1];
                    if (!cuts::hh_bbtautau_2016::hh_tag::IsInsideEllipse(event.SVfit_p4.mass(), bb.mass()))
                        continue;
                    tot_entries++;

                    std::uniform_int_distribution<uint_fast32_t> seed_distr(100000, std::numeric_limits<uint_fast32_t>::max());
                    uint_fast32_t test_split = test_vs_training(seed_gen);
                    gen.seed(seed_distr(seed_gen));
                    if (entry.id.IsBackground()) {
                        for (const auto mass : mass_range){
                            const SampleId sample_bkg(SampleType::Bkg_TTbar, mass);
                            vars->AddEvent(event, sample_bkg, entry.weight);
                            data[sample_bkg][test_split].push_back(reader->Evaluate(event, mass, args.method_name()));
                        }
                    }
                    else{
                        vars->AddEvent(event, entry.id, entry.weight);
                        data[entry.id][test_split].push_back(reader->Evaluate(event, entry.id.mass, args.method_name()));
                    }
                }
                std::cout << " channel " << mva_setup.channels[j] << "    " << entry.filename << " number of events: " << tot_entries << std::endl;
            }
        }

        std::map<SampleId, std::map<size_t,std::shared_ptr<TH1D>>> outputBDT = CreateOutputHistos(data);
        auto histo_roc = std::make_shared<TH2D>("ROC","ROC", static_cast<int>((args.max()-args.min())/9),args.min()-5, args.max()+5, 100, 0,1);
        std::map<int, double> roc;
        auto reader = vars->GetReader();
        auto method = dynamic_cast<TMVA::MethodBase*>(reader->FindMVA(args.method_name()));
        for (const auto& mass : mass_range){
            SampleId sample_sgn(SampleType::Sgn_Res, mass);
            SampleId sample_bkg(SampleType::Bkg_TTbar, mass);
            roc[mass] = method->GetROCIntegral(outputBDT.at(sample_sgn).at(tot).get(), outputBDT.at(sample_bkg).at(tot).get());
            histo_roc->Fill(mass, roc[mass]);
            std::cout<<mass<<"  ROC: "<<roc[mass]<<std::endl;
        }

        root_ext::WriteObject(*histo_roc, outfile.get());
        auto directory = root_ext::GetDirectory(*outfile.get(), "Kolmogorov");
        std::cout<<"kolmogorov"<<std::endl;
        std::map<SampleId, double> kolmogorov = Kolmogorov(data, directory);

        auto directory_sb = root_ext::GetDirectory(*outfile.get(), "Significance");
        std::cout<<"Significativity"<<std::endl;
        std::map<int, std::pair<double, PhysicalValue>> sign = EstimateSignificativity(mass_range, outputBDT, directory_sb, false);

        MvaTuple mva_tuple(outfile.get(), false);
        for (const auto& sample : kolmogorov){
            mva_tuple().KS_mass.push_back(sample.first.mass);
            mva_tuple().KS_value.push_back(sample.second);
            if (sample.first.IsSignal())
                mva_tuple().KS_type.push_back(1);
            else if (sample.first.IsBackground())
                mva_tuple().KS_type.push_back(-1);
        }

        for (const auto& value : roc){
            mva_tuple().roc_value.push_back(value.second);
            mva_tuple().roc_mass.push_back(value.first);
        }

        for (const auto& value : sign){
            mva_tuple().significance.push_back(value.second.second.GetValue());
            mva_tuple().significance_err.push_back(value.second.second.GetStatisticalError());
            mva_tuple().significance_mass.push_back(value.first);
        }
        mva_tuple.Fill();
        mva_tuple.Write();

    }
private:
    Arguments args;
    SampleEntryCollection samples;
    MvaSetup mva_setup;
    std::shared_ptr<TFile> outfile;
    std::mt19937_64 gen, gen2;
    MvaVariables::VarNameSet enabled_vars;
    std::shared_ptr<MvaReader> reader;
    std::uniform_int_distribution<uint_fast32_t>  seed_split, test_vs_training;
};

}
}

PROGRAM_MAIN(analysis::mva_study::MVAEvaluation, Arguments) // definition of the main program function


// ./run.sh MvaEvaluation ~/Desktop/tuples prova hh-bbtautau/Studies/config/mva_config.cfg ~/Desktop/xml/n/TMVAClassification_500t_PU_mass_newvars_HIGH_oldvars.weights.xml BDT::500t_PU_mass_newvars_HIGH_oldvars 2000 400 900 2 12345678 1234567 true false dphi_l1MET dphi_htautauMET dphi_hbbMET dphi_hbbhtautau dR_l1l2 dR_b1b2 MT_l1 MT_l2
// ./run.sh MvaEvaluation ~/Desktop/tuples prova_low hh-bbtautau/Studies/config/mva_config.cfg ~/Desktop/xml/n/TMVAClassification_500t_PU_mass_newvars_LOW.weights.xml BDT::500t_PU_mass_newvars_LOW 20000 250 350 2 12345678 1234567 true true dphi_l1MET dphi_htautauMET dphi_hbbMET dphi_hbbhtautau dR_l1l2Pt_htautau dR_b1b2Pt_hb MT_l1 MT_l2
// ./run.sh MvaEvaluation ~/Desktop/tuples prova_low hh-bbtautau/Studies/config/mva_config.cfg ~/workspace/hh-analysis/mydataloader/weight/myFactoryLowMassprova_Grad_shrinkage_0.1_12345678.weights.xml BDT::Grad_shrinkage_0.1_12345678 20000 250 320 2 12345678 1234567 false false HT_otherjets MT2 MT_htautau MT_l1 MT_tot abs_dphi_l1MET abs_dphi_htautauMET dR_b1b2_boosted dR_l1l2Pt_htautau dR_l1l2_boosted dphi_b1b2 dphi_hbbhtautau dphi_l2MET mass_l1l2MET mass_top1 mass_top2 p_zeta p_zetavisible pt_hbb pt_htautau pt_l1 pt_l1l2MET pt_l2
