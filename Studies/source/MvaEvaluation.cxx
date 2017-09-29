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
    OPT_ARG(std::string, error_file, "");
    OPT_ARG(bool, all_data, 1);
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
        outfile(root_ext::CreateRootFile(args.output_file()+".root")), gen(args.seed()), test_vs_training(0, args.number_sets()-1)
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
        samples = samples_list.at("Sample").files;
        std::cout<<"SAMPLES: "<<samples.size()<<std::endl;


        enabled_vars.insert(mva_setup.variables.begin(), mva_setup.variables.end());
        if(mva_setup.use_mass_var) {
            enabled_vars.insert("mass");
            enabled_vars.insert("channel");
            enabled_vars.insert("spin");
        }
        std::cout<<"VARS: "<<enabled_vars.size()<<std::endl;

    }

    void CreateOutputHistos(std::map<ChannelSampleIdSpin, std::map<size_t, std::vector<double>>> data, BDTData::Entry& outputBDT)
    {
        for(const auto& sample : data){
            std::cout<<sample.first.channel<<"    "<<sample.first.sample_id.sampleType<<"    "<<sample.first.sample_id.mass<<"    "<<sample.first.spin<<std::endl;
            for(const auto& entry : sample.second){

                std::vector<BDTData::Hist*> outs = { &outputBDT(sample.first.channel, sample.first.sample_id, sample.first.spin, entry.first)};
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
        std::map<ChannelSampleIdSpin, std::map<size_t, std::vector<double>>> data;
        std::cout<<"Variabili iniziali: "<<enabled_vars.size()<<std::endl;
        const auto mass_range = CreateMassRange(args.min(), args.max());
        std::cout<<"SAMPLES nel RANGE:"<<mass_range.size() <<std::endl;
        ::analysis::Range<int> range(args.min(), args.max());
        auto vars = reader.AddRange(range, args.method_name(), args.file_xml(), enabled_vars, args.isLegacy(), args.isLow());

        std::vector<ChannelSpin> set{{"muTau",0},{"eTau",0}, {"tauTau",0},{"muTau",2},{"eTau",2}, {"tauTau",2},
                                     {"tauTau",SM_spin}, {"muTau",SM_spin},{"eTau",SM_spin},
                                     {"muTau",bkg_spin},{"eTau",bkg_spin}, {"tauTau",bkg_spin}};

        for (const auto& s : set){
            std::cout << s.channel << s.spin <<std::endl;
            if (s.channel != args.channel()) continue;
            for(const SampleEntry& entry : samples)
            {
                if ( entry.id.IsSignal() && (entry.id.mass<args.min() || entry.id.mass>args.max())) continue;
                if ( entry.spin != s.spin) continue;
                if ( !entry.id.IsBackground() && entry.spin!=args.spin()) continue;

                auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
                auto tuple = ntuple::CreateEventTuple(s.channel, input_file.get(), true, ntuple::TreeState::Skimmed);
                auto sumtuple = ntuple::CreateSummaryTuple("summary",input_file.get(), true, ntuple::TreeState::Skimmed);
                auto mergesummary = ntuple::MergeSummaryTuple(*sumtuple.get());

                Long64_t tot_entries = 0;
                for(const Event& event : *tuple) {
                    if(tot_entries >= args.number_events()) break;
                    if (!args.all_data() && args.blind()!=(event.split_id >= mergesummary.n_splits/2)) continue;
                    tot_entries++;
                    size_t test_split = test_vs_training(gen);
                    if (entry.id.IsBackground()) {
                        for (const auto mass : mass_range){
                            SampleId sample_bkg(SampleType::Bkg_TTbar, mass);
                            ChannelSampleIdSpin id_ch_sample_spin{args.channel(), sample_bkg, args.spin()};
                            ChannelSampleIdSpin id_ch_bkg_spin{args.channel(), bkg, args.spin()};
                            double eval = reader.Evaluate(event, mass, args.method_name(), args.spin(), args.channel());
                            data[id_ch_sample_spin][test_split].push_back(eval);
                            data[id_ch_bkg_spin][test_split].push_back(eval);
                            data[id_ch_sample_spin][test_train].push_back(eval);
                            data[id_ch_bkg_spin][test_train].push_back(eval);
                        }
                    }
                    else{
                        double eval = reader.Evaluate(event, entry.id.mass, args.method_name(), args.spin(), args.channel());
                        ChannelSampleIdSpin id_ch_sample_spin{args.channel(), entry.id, args.spin()};
                        ChannelSampleIdSpin id_ch_tot_spin{args.channel(), mass_tot, args.spin()};
                        data[id_ch_sample_spin][test_split].push_back(eval);
                        data[id_ch_tot_spin][test_split].push_back(eval);
                        data[id_ch_sample_spin][test_train].push_back(eval);
                        data[id_ch_tot_spin][test_train].push_back(eval);
                    }
                }
                std::cout << " channel " << s.channel << "    " << entry.filename << " number of events: " << tot_entries << std::endl;
            }
        }

        std::map<ChannelSampleIdSpin, PhysicalValue> err_training, err_testing;
        if (args.error_file().size()){
            std::shared_ptr<TFile> in_file(root_ext::OpenRootFile(args.error_file()));
            MvaTuple myTree("mva_result", in_file.get(), true);
            for(const MvaResults& results : myTree) {
                err_training = GetRocTrainingIntegralMap(results);
                err_testing = GetRocTestingIntegralMap(results);
            }
        }

        BDTData outputBDT(outfile, "Evaluation");
        CreateOutputHistos(data, outputBDT.bdt_out);

//        BDTData difference(outfile, "Difference");
        auto histo_roc = std::make_shared<TGraphErrors>();
        histo_roc->SetTitle("ROC");
        histo_roc->SetName("ROC");
        std::map<int, double> roc, roc_testing, roc_training;
        auto reader_method = vars->GetReader();
        auto method = dynamic_cast<TMVA::MethodBase*>(reader_method->FindMVA(args.method_name()));
        int i = 0;
        auto directory_roc = root_ext::GetDirectory(*outfile.get(), "ROCCurve");
        for (const auto& mass : mass_range){
            std::cout<<mass<<std::endl;
            SampleId sample_sgn(SampleType::Sgn_Res, mass);
            SampleId sample_bkg(SampleType::Bkg_TTbar, mass);
            ChannelSampleIdSpin id_sgn(args.channel(), sample_sgn, args.spin());
            ChannelSampleIdSpin id_bkg(args.channel(), sample_bkg, args.spin());
            roc[mass] = method->GetROCIntegral(&outputBDT.bdt_out(args.channel(), sample_sgn, args.spin(), test_train), &outputBDT.bdt_out(args.channel(),sample_bkg,args.spin(), test_train));
            roc_testing[mass] = method->GetROCIntegral(&outputBDT.bdt_out(args.channel(), sample_sgn, args.spin(), 0), &outputBDT.bdt_out(args.channel(),sample_bkg,args.spin(),0));
            roc_training[mass] = method->GetROCIntegral(&outputBDT.bdt_out(args.channel(), sample_sgn, args.spin(), 1), &outputBDT.bdt_out(args.channel(),sample_bkg,args.spin(),1));
            histo_roc->SetPoint(i, mass, roc[mass]);
            histo_roc->SetPointError(i,0, err_testing[id_sgn].GetFullError());

            std::cout<<mass<<"  ROC: "<<roc[mass] << " +- "<<err_testing[id_sgn].GetFullError()<<std::endl;
            i++;
            auto directory_roc_mass = root_ext::GetDirectory(*directory_roc, std::to_string(mass));
            std::vector<float> mvaS, mvaB;
            for (const auto& eval : data[id_sgn][test_train])
                mvaS.push_back(static_cast<float>(eval));
            for (const auto& eval : data[id_bkg][test_train])
                mvaB.push_back(static_cast<float>(eval));
            TMVA::ROCCurve roccurve(mvaS, mvaB);
            auto graph = roccurve.GetROCCurve();
            root_ext::WriteObject(*graph, directory_roc_mass);

        }
        root_ext::WriteObject(*histo_roc, outfile.get());

        auto roctot = method->GetROCIntegral(&outputBDT.bdt_out(args.channel(), mass_tot, args.spin(), test_train),
                                             &outputBDT.bdt_out(args.channel(),  bkg, args.spin(), test_train));
//        auto roctot_testing = method->GetROCIntegral(&outputBDT.bdt_out(args.channel(), mass_tot, args.spin(), 0),
//                                                     &outputBDT.bdt_out(args.channel(),  bkg, args.spin(), 0));
//        auto roctot_training = method->GetROCIntegral(&outputBDT.bdt_out(args.channel(), mass_tot, args.spin(), 1),
//                                                      &outputBDT.bdt_out(args.channel(),  bkg, args.spin(), 1));

        std::cout<<roctot<<std::endl;

        std::vector<float> mvaS, mvaB;
        ChannelSampleIdSpin id_mass_tot{args.channel(), mass_tot,args.spin()};
        ChannelSampleIdSpin id_bkg{args.channel(),bkg,args.spin()};
        for (const auto& eval : data[id_mass_tot][test_train])
            mvaS.push_back(static_cast<float>(eval));
        for (const auto& eval : data[id_bkg][test_train])
            mvaB.push_back(static_cast<float>(eval));

        TMVA::ROCCurve roccurve(mvaS, mvaB);
        auto graph = roccurve.GetROCCurve();
        root_ext::WriteObject(*graph, directory_roc);

        auto directory = root_ext::GetDirectory(*outfile.get(), "Kolmogorov");
        std::cout<<"kolmogorov"<<std::endl;
        auto kolmogorov = KolmogorovTest(data, outputBDT.bdt_out,directory,true);
        auto directory_chi = root_ext::GetDirectory(*outfile.get(), "Chi2");
        std::cout<<"chi2"<<std::endl;
        auto chi = ChiSquareTest(data, outputBDT.bdt_out, directory_chi,true);

        auto directory_sb = root_ext::GetDirectory(*outfile.get(), "Significance");
        std::cout<<"Significativity"<<std::endl;
        auto sign = EstimateSignificativity(args.channel(), args.spin(), mass_range, outputBDT.bdt_out, directory_sb, false);

        MvaTuple mva_tuple(outfile.get(), false);
        for (const auto& sample : kolmogorov){
            mva_tuple().KS_mass.push_back(sample.first.sample_id.mass);
            mva_tuple().KS_value.push_back(sample.second);
            mva_tuple().KS_type.push_back(static_cast<int>(sample.first.sample_id.sampleType));
            mva_tuple().KS_spin.push_back(static_cast<int>(sample.first.spin));
            mva_tuple().KS_channel.push_back(sample.first.channel);
        }

        for (const auto& sample : chi){
            mva_tuple().chi_mass.push_back(sample.first.sample_id.mass);
            mva_tuple().chi_value.push_back(sample.second);
            mva_tuple().chi_type.push_back(static_cast<int>(sample.first.sample_id.sampleType));
            mva_tuple().chi_spin.push_back(static_cast<int>(sample.first.spin));
            mva_tuple().chi_channel.push_back(sample.first.channel);
        }

        for (const auto& value : roc){
            mva_tuple().roc_testing_value.push_back(value.second);
            mva_tuple().roc_testing_mass.push_back(value.first);
            mva_tuple().roc_testing_channel.push_back(args.channel());
            mva_tuple().roc_testing_spin.push_back(args.spin());
            mva_tuple().roc_testing_type.push_back(static_cast<int>(SampleType::Sgn_Res));
        }

        mva_tuple().roc_testing_value.push_back(roctot);
        mva_tuple().roc_testing_mass.push_back(mass_tot.mass);
        mva_tuple().roc_testing_channel.push_back(args.channel());
        mva_tuple().roc_testing_spin.push_back(args.spin());
        mva_tuple().roc_testing_type.push_back(static_cast<int>(SampleType::Sgn_Res));

        for (const auto& value : sign){
            mva_tuple().optimal_cut.push_back(value.second.cut);
            mva_tuple().significance.push_back(value.second.significance.GetValue());
            mva_tuple().significance_err.push_back(value.second.significance.GetStatisticalError());
            mva_tuple().significance_mass.push_back(value.first.sample_id.mass);
            mva_tuple().significance_type.push_back(static_cast<int>(value.first.sample_id.sampleType));
            mva_tuple().significance_spin.push_back(value.first.spin);
            mva_tuple().significance_channel.push_back(value.first.channel);
        }
        mva_tuple.Fill();
        mva_tuple.Write();

    }
private:
    Arguments args;
    SampleEntryCollection samples;
    MvaSetup mva_setup;
    std::shared_ptr<TFile> outfile;
    std::mt19937_64 gen;
    MvaVariables::VarNameSet enabled_vars;
    MvaReader reader;
    std::uniform_int_distribution<size_t>  test_vs_training;
};

}
}

PROGRAM_MAIN(analysis::mva_study::MVAEvaluation, Arguments) // definition of the main program function




