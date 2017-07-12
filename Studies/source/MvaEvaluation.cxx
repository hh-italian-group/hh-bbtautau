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
#include "hh-bbtautau/Studies/include/MvaVariablesStudy.h"
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
    REQ_ARG(std::string, range);
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
        mva_setup = setups.at(args.range());

        if(!samples_list.count("Samples"))
            throw exception("Samples don't found");
        samples = samples_list.at("Samples").files;
        std::cout<<"SAMPLES: "<<samples.size()<<std::endl;


        enabled_vars.insert(mva_setup.variables.begin(), mva_setup.variables.end());
        if(mva_setup.use_mass_var) {
            enabled_vars.insert("mass");
            enabled_vars.insert("channel");
        }
        std::cout<<"VARS: "<<enabled_vars.size()<<std::endl;

    }

    void CreateOutputHistos(std::map<SampleId, std::map<size_t, std::vector<double>>> data, BDTData::Entry& outputBDT)
    {
        for(const auto& sample : data){
            for(const auto& entry : sample.second){
                std::vector<BDTData::Hist*> outs = { &outputBDT(sample.first, tot), &outputBDT(sample.first, entry.first) };
                for (const auto& value : entry.second){
                    for(auto out : outs)
                        out->Fill(value);
                }
            }
        }
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
        auto vars = reader.AddRange(range, args.method_name(), args.file_xml(), enabled_vars, args.isLegacy(), args.isLow());
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
                            data[sample_bkg][test_split].push_back(reader.Evaluate(event, mass, args.method_name()));
                        }
                    }
                    else{
                        vars->AddEvent(event, entry.id, entry.weight);
                        data[entry.id][test_split].push_back(reader.Evaluate(event, entry.id.mass, args.method_name()));
                    }
                }
                std::cout << " channel " << mva_setup.channels[j] << "    " << entry.filename << " number of events: " << tot_entries << std::endl;
            }
        }

        BDTData outputBDT(outfile, "Evaluation");
        CreateOutputHistos(data, outputBDT.bdt_out);
        auto histo_roc = std::make_shared<TH2D>("ROC","ROC", 66, 245, 905, 100, 0,1);
        std::map<int, double> roc;
        auto reader_method = vars->GetReader();
        auto method = dynamic_cast<TMVA::MethodBase*>(reader_method->FindMVA(args.method_name()));
        for (const auto& mass : mass_range){
            SampleId sample_sgn(SampleType::Sgn_Res, mass);
            SampleId sample_bkg(SampleType::Bkg_TTbar, mass);
            roc[mass] = method->GetROCIntegral(&outputBDT.bdt_out(sample_sgn,tot), &outputBDT.bdt_out(sample_bkg,tot));
            histo_roc->Fill(mass, roc[mass]);
            std::cout<<mass<<"  ROC: "<<roc[mass]<<std::endl;
        }

        root_ext::WriteObject(*histo_roc, outfile.get());
        auto directory = root_ext::GetDirectory(*outfile.get(), "Kolmogorov");
        std::cout<<"kolmogorov"<<std::endl;
        std::map<SampleId, double> kolmogorov = Kolmogorov(data, directory);

        auto directory_sb = root_ext::GetDirectory(*outfile.get(), "Significance");
        std::cout<<"Significativity"<<std::endl;
        std::map<int, std::pair<double, PhysicalValue>> sign = EstimateSignificativity(mass_range, outputBDT.bdt_out, directory_sb, false);

        MvaTuple mva_tuple(outfile.get(), false);
        for (const auto& sample : kolmogorov){
            mva_tuple().KS_mass.push_back(sample.first.mass);
            mva_tuple().KS_value.push_back(sample.second);
            mva_tuple().KS_type.push_back(static_cast<int>(sample.first.sampleType));
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
    MvaReader reader;
    std::uniform_int_distribution<uint_fast32_t>  seed_split, test_vs_training;
};

}
}

PROGRAM_MAIN(analysis::mva_study::MVAEvaluation, Arguments) // definition of the main program function




