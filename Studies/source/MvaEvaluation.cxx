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
    REQ_ARG(size_t, number_sets);
    REQ_ARG(uint_fast32_t, seed);
    REQ_ARG(bool, isLegacy);
    REQ_ARG(bool, isLow);
    REQ_ARG(std::string, range);
    REQ_ARG(std::string, channel);
    REQ_ARG(int, spin);
    OPT_ARG(bool, blind, 1);
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
        outfile(root_ext::CreateRootFile(args.output_file()+".root")), test_vs_training(0, args.number_sets()-1)
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
            else if (sample.id.IsSM())
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

        std::vector<ChannelSpin> set{{"muTau",0},{"eTau",0},{"tauTau",0},{"muTau",2},{"eTau",2},{"tauTau",2},{"muTau",1},{"eTau",1},{"tauTau",1},{"muTau",-1},{"eTau",-1},{"tauTau",-1}};



        for (const auto& s : set){
            std::cout << s.first << s.second <<std::endl;
            if (s.first != args.channel()) continue;
            for(const SampleEntry& entry : samples)
            {
                if ( entry.id.IsSignal() && (entry.id.mass<args.min() || entry.id.mass>args.max())) continue;
                if ( entry.spin != s.second) continue;
                if ( !entry.id.IsBackground() && entry.spin!=args.spin()) continue;

                auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
                auto tuple = ntuple::CreateEventTuple(s.first, input_file.get(), true, ntuple::TreeState::Skimmed);
                ntuple::SummaryTuple sumtuple("summary",input_file.get(), true);

                Long64_t tot_entries = 0;
                for (Long64_t  current_entry = 0; current_entry < tuple->GetEntries(); current_entry++) {

                    tuple->GetEntry(current_entry);
                    const Event& event = tuple->data();
                    if(tot_entries >= args.number_events()) break;
                    sumtuple.GetEntry(current_entry);
                    if (args.blind())
                        if (event.split_id >= (sumtuple.data().n_splits/2)) continue;
                    if (!args.blind())
                        if (event.split_id < (sumtuple.data().n_splits/2)) continue;
                    tot_entries++;
                    std::uniform_int_distribution<uint_fast32_t> seed_distr(100000, std::numeric_limits<uint_fast32_t>::max());
                    size_t test_split = test_vs_training(seed_gen);
                    gen.seed(seed_distr(seed_gen));
                    if (entry.id.IsBackground()) {
                        for (const auto mass : mass_range){
                            const SampleId sample_bkg(SampleType::Bkg_TTbar, mass);
                            vars->AddEvent(event, sample_bkg, entry.spin, s.first, entry.weight);
                            double eval = reader.Evaluate(event, mass, args.method_name(), entry.spin, s.first);
                            data[sample_bkg][test_split].push_back(eval);
                            data[bkg][test_split].push_back(eval);
                            data[sample_bkg][test_train].push_back(eval);
                            data[bkg][test_train].push_back(eval);
                        }
                    }
                    else{
                        vars->AddEvent(event, entry.id, entry.spin, s.first, entry.weight);
                        double eval = reader.Evaluate(event, entry.id.mass, args.method_name(), entry.spin, s.first);
                        data[entry.id][test_split].push_back(eval);
                        data[mass_tot][test_split].push_back(eval);
                        data[entry.id][test_train].push_back(eval);
                        data[mass_tot][test_train].push_back(eval);
                    }
                }
                std::cout << " channel " << s.first << "    " << entry.filename << " number of events: " << tot_entries << std::endl;
            }
        }

        BDTData outputBDT(outfile, "Evaluation");
        CreateOutputHistos(data, outputBDT.bdt_out);
        BDTData difference(outfile, "Difference");
        auto histo_roc = std::make_shared<TGraph>();
        histo_roc->SetTitle("ROC");
        histo_roc->SetName("ROC");
        std::map<int, double> roc;
        auto reader_method = vars->GetReader();
        auto method = dynamic_cast<TMVA::MethodBase*>(reader_method->FindMVA(args.method_name()));
        int i = 0;
        auto directory_roc = root_ext::GetDirectory(*outfile.get(), "ROCCurve");
        for (const auto& mass : mass_range){
            std::cout<<mass<<std::endl;
            SampleId sample_sgn(SampleType::Sgn_Res, mass);
            SampleId sample_bkg(SampleType::Bkg_TTbar, mass);
            roc[mass] = method->GetROCIntegral(&outputBDT.bdt_out(sample_sgn,tot), &outputBDT.bdt_out(sample_bkg,tot));
            histo_roc->SetPoint(i, mass, roc[mass]);
            std::cout<<mass<<"  ROC: "<<roc[mass]<<std::endl;
            i++;
            auto directory_roc_mass = root_ext::GetDirectory(*directory_roc, std::to_string(mass));
            std::vector<float> mvaS, mvaB;
            for (const auto& eval : data[sample_sgn][test_train])
                mvaS.push_back(static_cast<float>(eval));
            for (const auto& eval : data[sample_bkg][test_train])
                mvaB.push_back(static_cast<float>(eval));
            TMVA::ROCCurve roccurve(mvaS, mvaB);
            auto graph = roccurve.GetROCCurve();
            root_ext::WriteObject(*graph, directory_roc_mass);

        }
        root_ext::WriteObject(*histo_roc, outfile.get());

        auto roctot = method->GetROCIntegral(&outputBDT.bdt_out(mass_tot,tot), &outputBDT.bdt_out(bkg,tot));
        std::cout<<roctot<<std::endl;

        std::vector<float> mvaS, mvaB;
        for (const auto& eval : data[mass_tot][test_train])
            mvaS.push_back(static_cast<float>(eval));
        for (const auto& eval : data[bkg][test_train])
            mvaB.push_back(static_cast<float>(eval));

        TMVA::ROCCurve roccurve(mvaS, mvaB);
        auto graph = roccurve.GetROCCurve();
        root_ext::WriteObject(*graph, directory_roc);

        auto directory = root_ext::GetDirectory(*outfile.get(), "Kolmogorov");
        std::cout<<"kolmogorov"<<std::endl;
        std::map<SampleId, double> kolmogorov = Kolmogorov(data, outputBDT.bdt_out,  difference.difference ,directory);
        auto directory_chi = root_ext::GetDirectory(*outfile.get(), "Chi2");
        std::cout<<"chi2"<<std::endl;
        std::map<SampleId, double> chi = ChiSquare(data, outputBDT.bdt_out, directory_chi);

        auto directory_sb = root_ext::GetDirectory(*outfile.get(), "Significance");
        std::cout<<"Significativity"<<std::endl;
        std::map<int, std::pair<double, PhysicalValue>> sign = EstimateSignificativity(mass_range, outputBDT.bdt_out, directory_sb, false);

        MvaTuple mva_tuple(outfile.get(), false);
        for (const auto& sample : kolmogorov){
            mva_tuple().KS_mass.push_back(sample.first.mass);
            mva_tuple().KS_value.push_back(sample.second);
            mva_tuple().KS_type.push_back(static_cast<int>(sample.first.sampleType));
        }

        for (const auto& sample : chi){
            mva_tuple().chi_mass.push_back(sample.first.mass);
            mva_tuple().chi_value.push_back(sample.second);
            mva_tuple().chi_type.push_back(static_cast<int>(sample.first.sampleType));
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
    std::uniform_int_distribution<size_t>  test_vs_training;
};

}
}

PROGRAM_MAIN(analysis::mva_study::MVAEvaluation, Arguments) // definition of the main program function




