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
#include "TMVA/ROCCurve.h"
#include "TMVA/ClassifierFactory.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/ResultsClassification.h"
#include "TMVA/ROCCurve.h"
#include "TMath.h"
#include <fstream>
#include <random>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include "hh-bbtautau/Analysis/include/MvaConfigurationReader.h"
#include "hh-bbtautau/Studies/include/MvaTuple.h"
#include "hh-bbtautau/Studies/include/MvaVariablesStudy.h"
#include "hh-bbtautau/Studies/include/MvaMethods.h"
#include "hh-bbtautau/McCorrections/include/HH_nonResonant_weight.h"

struct Arguments { // list of all program arguments
        REQ_ARG(std::string, input_path);
        REQ_ARG(std::string, output_file);
        REQ_ARG(std::string, cfg_file);
        REQ_ARG(Long64_t, number_events);
        REQ_ARG(std::string, range);
        REQ_ARG(Long64_t, number_variables);
        REQ_ARG(uint_fast32_t, which_test);
        REQ_ARG(std::string, suffix);
        OPT_ARG(size_t, number_sets, 2);
        OPT_ARG(uint_fast32_t, seed, std::numeric_limits<uint_fast32_t>::max());
        OPT_ARG(std::string, save, "");
        OPT_ARG(bool, all_data, 1);
        OPT_ARG(bool, blind, 1);
        OPT_ARG(uint_fast32_t, subdivisions, 2);
        OPT_ARG(bool, is_SM, false);
        OPT_ARG(bool, is_BSM, false);
        OPT_ARG(std::string, coeffFile, "");
        OPT_ARG(std::string, input_histo, "");
};

namespace analysis {
namespace mva_study{

class MvaVariablesTMVA : public MvaVariables {
public:
    using DataVector = std::vector<double>;
    using DataVectorF = std::vector<float>;
    struct SampleData{
        double sampleweight;
        std::deque<std::pair<DataVector, double>> data;
    };
    using MassData_pair = std::map<ChannelSampleIdSpin, SampleData>;
    static constexpr size_t max_n_vars = 1000;
    DataVector variable;
    DataVectorF variable_float;
    std::map<size_t, MassData_pair> data_pair;
    std::map<std::string, size_t> name_indices;
    std::vector<std::string> names;
    std::shared_ptr<TMVA::DataLoader> loader;
    std::shared_ptr<TMVA::Reader> reader;
    MvaVariablesTMVA(size_t _number_set = 2, uint_fast32_t _seed = std::numeric_limits<uint_fast32_t>::max(),
                            const VarNameSet& _enabled_vars = {}) :
        MvaVariables(_number_set, _seed, _enabled_vars), variable_float(max_n_vars),
        loader(new TMVA::DataLoader("mydataloader")), reader(new TMVA::Reader)
    {}

    virtual void SetValue(const std::string& name, double value, char type = 'F') override
    {
        if (!name_indices.count(name)){
            variable.push_back(0);
            name_indices[name] = variable.size()-1;
            names.push_back(name);
            loader->AddVariable(name, type);
            reader->AddVariable(name, &variable_float.at(variable.size()-1));
        }
        variable.at(name_indices.at(name)) = value;
    }

    virtual void AddEventVariables(size_t istraining, const SampleId& mass, double weight, double sampleweight, int spin, std::string channel) override
    {
        ChannelSampleIdSpin id{channel, mass, spin};
        data_pair[istraining][id].data.emplace_back(variable, weight);
        data_pair[istraining][id].sampleweight = sampleweight;
    }

    void UploadEvents()
    {
        std::map< ChannelSampleIdSpin, double> weights;
        for(const auto& entry : data_pair[0]) {
            const  ChannelSampleIdSpin id = entry.first;
            double tot_weight = 0;
            for (const auto& val : data_pair){
                for (const auto& x : val.second.at(id).data){
                    tot_weight+= x.second;
                }
            }
            weights[id] = entry.second.sampleweight / (tot_weight);
        }
        std::cout<<"pesi"<<std::endl;
        for(const auto& data_entry : data_pair) {
            const bool istraining = data_entry.first;
            const TMVA::Types::ETreeType treetype = istraining ? TMVA::Types::kTraining : TMVA::Types::kTesting;
            for(const auto& m_entry : data_entry.second) {
                const  ChannelSampleIdSpin id = m_entry.first;
                const SampleId m = id.sample_id;
                const std::string samplename = m.IsSignal() ? "Signal" : "Background";
                for(const auto& vars : m_entry.second.data) {
                    loader->AddEvent(samplename, treetype, vars.first , vars.second*weights.at(id));
                }
            }
        }
    }

    void SaveEvents(const std::string& file_name)
    {
        static const std::string sep = ",";
        std::stringstream ss;
        std::ofstream save(file_name, std::ofstream::out);
        ss << "index" << sep << "IsSignal" << sep;
        for (const auto& var: names){
            ss << var << sep ;
        }
        ss << "Weight" << std::endl;
        size_t n = 0;
        for(const auto& data_entry : data_pair) {
            for(const auto& m_entry : data_entry.second) {
                const SampleId m = m_entry.first.sample_id;
                for(const auto& vars : m_entry.second.data) {
                    ss << n << sep <<  m.IsSignal() << sep;
                    for (const auto& var : vars.first)
                        ss << var << sep;
                    ss << vars.second <<std::endl;
                    ++n;
                }
            }
        }
        boost::iostreams::filtering_streambuf< boost::iostreams::input> in;
        in.push( boost::iostreams::gzip_compressor());
        in.push( ss );
        boost::iostreams::copy(in, save);
    }


    double Evaluation(const std::string& method_name, const std::vector<double>& event)
    {
        if(event.size() != variable.size())
            throw exception("Invalid event size.");
        for(size_t n = 0; n < variable.size(); ++n)
            variable_float[n] = static_cast<float>(event[n]);
        return reader->EvaluateMVA(method_name);
    }

    std::vector<std::pair<double,double>> EvaluateForAllEvents(const std::string& method_name, size_t istraining, const std::string channel, const int& spin,
                                                               const SampleId& sample)
    {
        ChannelSampleIdSpin id{channel, sample, spin};
        if(!data_pair.count(istraining) || !data_pair.at(istraining).count(id))
            throw exception("Sample not found.");
        const auto& events = data_pair.at(istraining).at(id).data;
        std::vector<std::pair<double,double>> result;
        result.reserve(events.size());
        for(const auto& event : events)
            result.emplace_back(Evaluation(method_name, event.first), event.second);
        return result;
    }

    virtual std::shared_ptr<TMVA::Reader> GetReader() override { return reader;}
};


class MVATraining {
public:
    using Event = ::ntuple::Event;
    using EventTuple = ::ntuple::EventTuple;
    using SummaryTuple = ntuple::SummaryTuple;

    std::vector<ChannelSpin> set_SM{{"tauTau",SM_spin}, {"muTau",SM_spin},{"eTau",SM_spin},
                                 {"muTau",bkg_spin},{"eTau",bkg_spin}, {"tauTau",bkg_spin}};
    std::vector<ChannelSpin> set_R{{"tauTau",0}, {"muTau",0},{"eTau",0},
                                   {"tauTau",2}, {"muTau",2},{"eTau",2},
                                 {"muTau",bkg_spin},{"eTau",bkg_spin}, {"tauTau",bkg_spin}};
    std::vector<ChannelSpin> set;

    MVATraining(const Arguments& _args): args(_args),
        outfile(root_ext::CreateRootFile(args.output_file()+"_"+std::to_string(args.which_test())+".root"))
    {
        MvaSetupCollection setups;
        SampleEntryListCollection samples_list;

        ConfigReader configReader;
        MvaConfigReader setupReader(setups);
        configReader.AddEntryReader("SETUP", setupReader, true);
        SampleConfigReader sampleReader(samples_list);
        configReader.AddEntryReader("FILES", sampleReader, false);
        configReader.ReadConfig(args.cfg_file());

        if(!samples_list.count(args.range()))
            throw exception("Samples at '%1%' mass don't found") %args.range();
        samples = samples_list.at(args.range()).files;

        if(!setups.count(args.range()+std::to_string(args.number_variables())))
            throw exception("Setups at '%1%' range don't found") %args.range();
        mva_setup = setups.at(args.range()+std::to_string(args.number_variables()));
        std::cout<<args.range()+std::to_string(args.number_variables())<<std::endl;
        enabled_vars.insert(mva_setup.variables.begin(), mva_setup.variables.end());
        std::cout<<enabled_vars.size()<<std::endl;

        if(mva_setup.use_mass_var) {
            enabled_vars.insert("mass");
            enabled_vars.insert("channel");
            enabled_vars.insert("spin");
        }

        std::cout<<"quante variabii? "<<enabled_vars.size()<<std::endl;
        std::mt19937_64 seed_gen(args.seed());
        std::uniform_int_distribution<uint_fast32_t> seed_distr(100000, std::numeric_limits<uint_fast32_t>::max());
        uint_fast32_t seed = seed_distr(seed_gen);
        gen.seed(seed_distr(seed_gen));
        gen2.seed(seed_distr(seed_gen));
        vars = std::make_shared<MvaVariablesTMVA>(args.number_sets(), seed, enabled_vars);

        set = args.is_SM() ? set_SM : set_R;

        if (args.is_BSM()) {
            reweight5D = std::make_shared<NonResHH_EFT::WeightProvider>(args.coeffFile());
//            benchmarks = NonResHH_EFT::BenchmarkCollection::Benchmarks();
            const std::vector<std::vector<double>> kl_scan = {
                { -20, -10, -5, -1, 0, 1, 2, 5, 10, 30 }, { 1 }, { 0 }, { 0 }, { 0 }
            };
            benchmarks = NonResHH_EFT::BenchmarkCollection::CreateScanBenchmarks(kl_scan);
        }
    }

    void TrainAllMethods(const TMVA::Factory& factory)
    {
       std::map< TString, TMVA::Factory::MVector*> fMethodsMap = factory.fMethodsMap;
       std::map<TString,TMVA::Factory::MVector*>::iterator itrMap;
       for(itrMap = fMethodsMap.begin(); itrMap != fMethodsMap.end(); itrMap++)
       {
          TMVA::Factory::MVector *methods=itrMap->second;
          TMVA::Factory::MVector::iterator itrMethod;

          for( itrMethod = methods->begin(); itrMethod != methods->end(); itrMethod++ ) {
              TMVA::Event::SetIsTraining(kTRUE);
              TMVA::MethodBase* mva = dynamic_cast<TMVA::MethodBase*>(*itrMethod);
              if (mva==0)
                  throw exception("Method not found.");
              if (mva->DataInfo().GetNClasses() < 2 )
                  throw exception("You want to do classification training, but specified less than two classes.");
              if (mva->Data()->GetNTrainingEvents() < 10) {
                throw exception("Method '%1%' not trained (training tree has less entries ['%2%'] than required [10]).")
                          % mva->GetMethodName() % mva->Data()->GetNTrainingEvents();
              }
              mva->TrainMethod();
          }
       }
    }

    void RankingImportanceVariables(const std::vector<std::pair<std::string, double>>& importance,
                                    std::map<std::string, std::shared_ptr<TH1D>>& histo_rank) const
    {
       for (const auto& pair: importance){
            if (!histo_rank.count(pair.first)){
                histo_rank[pair.first] = std::make_shared<TH1D>((pair.first+"_importance").c_str(),(pair.first+"_importance").c_str(),100, 0, 0.1);
                histo_rank.at(pair.first)->SetCanExtend(TH1::kAllAxes);
                histo_rank.at(pair.first)->SetXTitle("importance");
            }
            histo_rank.at(pair.first)->Fill(pair.second);
        }
    }

    std::map<std::string, size_t> RankingPoisitionVariables(std::vector<std::pair<std::string, double>> importance,
                                                            std::map<std::string, std::shared_ptr<TH1D>>& histo_rank) const
    {
        std::map<std::string, size_t> position;

        std::sort(importance.begin(), importance.end(), [](auto el1, auto el2){
            return el1.second > el2.second;
        } );
        for (size_t i = 1; i <= importance.size(); i++){
            position[importance[i-1].first] = i;
        }

        for (const auto& pair: position){
           if (!histo_rank.count(pair.first)){
               histo_rank[pair.first] = std::make_shared<TH1D>((pair.first+"_position").c_str(),(pair.first+"_position").c_str(), position.size(), 0.5, position.size() + 0.5);
               histo_rank[pair.first]->SetXTitle("position");
           }
            histo_rank[pair.first]->Fill(pair.second);
        }
        return position;
    }

    void EvaluateMethod(std::map<ChannelSampleIdSpin, std::map<size_t, std::vector<double>>>& evaluation,
                        BDTData::Entry& outputBDT, const std::string& method_name)
    {
        std::string weightfile = "mydataloader/weights/myFactory"+args.output_file()+"_"+method_name+".weights.xml";
        vars->reader->BookMVA(method_name, weightfile);

        for(const auto& type_entry : vars->data_pair){
            const auto& sample_entry = type_entry.second;
            std::cout<<type_entry.first<<std::endl;
            for(const auto& entry : sample_entry){

                auto current_evaluation = vars->EvaluateForAllEvents(method_name, type_entry.first, entry.first.channel,
                                                                     entry.first.spin, entry.first.sample_id);

                ChannelSampleIdSpin allch_sample_spin{all_channel, entry.first.sample_id, entry.first.spin};

                ChannelSampleIdSpin allch_sample_allspin{all_channel, entry.first.sample_id, spin_tot};
                ChannelSampleIdSpin ch_sample_allspin{entry.first.channel, entry.first.sample_id, spin_tot};
                if (entry.first.sample_id.IsBackground()){
                    allch_sample_allspin.spin = bkg_spin;
                    ch_sample_allspin.spin = bkg_spin;
                }

                for (const auto& val: current_evaluation){
                    evaluation[entry.first][type_entry.first].push_back(val.first);
                    evaluation[allch_sample_spin][type_entry.first].push_back(val.first);
                    evaluation[allch_sample_allspin][type_entry.first].push_back(val.first);
                    evaluation[ch_sample_allspin][type_entry.first].push_back(val.first);
                }

                ChannelSampleIdSpin allch_allsample_allspin{all_channel, mass_tot, spin_tot};
                ChannelSampleIdSpin allch_allsample_spin{all_channel, mass_tot, entry.first.spin};

                ChannelSampleIdSpin ch_allsample_allspin{entry.first.channel, mass_tot, spin_tot};
                ChannelSampleIdSpin ch_allsample_spin{entry.first.channel, mass_tot, entry.first.spin};

                if (entry.first.sample_id.IsBackground()){
                    allch_allsample_allspin.sample_id = bkg;
                    allch_allsample_allspin.spin = bkg_spin;

                    ch_allsample_allspin.sample_id = bkg;
                    ch_allsample_allspin.spin = bkg_spin;

                    allch_allsample_spin.sample_id = bkg;
                    allch_allsample_spin.spin = bkg_spin;

                    ch_allsample_spin.sample_id = bkg;
                    ch_allsample_spin.spin = bkg_spin;
                }

                auto& eval_aaa = evaluation[allch_allsample_allspin][type_entry.first];
                eval_aaa.insert(eval_aaa.end(), evaluation[entry.first][type_entry.first].begin(), evaluation[entry.first][type_entry.first].end());

                auto& eval_caa = evaluation[ch_allsample_allspin][type_entry.first];
                eval_caa.insert(eval_caa.end(), evaluation[entry.first][type_entry.first].begin(), evaluation[entry.first][type_entry.first].end());

                if (entry.first.sample_id.IsSignal()){
                    auto& eval_aas = evaluation[allch_allsample_spin][type_entry.first];
                    eval_aas.insert(eval_aas.end(), evaluation[entry.first][type_entry.first].begin(), evaluation[entry.first][type_entry.first].end());
                    auto& eval_cas = evaluation[ch_allsample_spin][type_entry.first];
                    eval_cas.insert(eval_cas.end(), evaluation[entry.first][type_entry.first].begin(), evaluation[entry.first][type_entry.first].end());
                }

                std::set<BDTData::Hist*> outs = { &outputBDT(entry.first.channel, entry.first.sample_id, entry.first.spin, type_entry.first),
                                                  &outputBDT(entry.first.channel, entry.first.sample_id, entry.first.spin, tot),
                                                  &outputBDT(all_channel, entry.first.sample_id, entry.first.spin, type_entry.first),
                                                  &outputBDT(all_channel, entry.first.sample_id, entry.first.spin, tot),
                                                  &outputBDT(entry.first.channel, entry.first.sample_id, spin_tot, type_entry.first),
                                                  &outputBDT(entry.first.channel, entry.first.sample_id, spin_tot, tot),
                                                  &outputBDT(all_channel, entry.first.sample_id, spin_tot, type_entry.first),
                                                  &outputBDT(all_channel, entry.first.sample_id, spin_tot, tot),

                                                  &outputBDT(ch_allsample_spin.channel, ch_allsample_spin.sample_id, ch_allsample_spin.spin, type_entry.first),
                                                  &outputBDT(ch_allsample_spin.channel, ch_allsample_spin.sample_id, ch_allsample_spin.spin, tot),
                                                  &outputBDT(allch_allsample_spin.channel, allch_allsample_spin.sample_id, allch_allsample_spin.spin, type_entry.first),
                                                  &outputBDT(allch_allsample_spin.channel, allch_allsample_spin.sample_id, allch_allsample_spin.spin, tot),
                                                  &outputBDT(ch_allsample_allspin.channel, ch_allsample_allspin.sample_id, ch_allsample_allspin.spin, type_entry.first),
                                                  &outputBDT(ch_allsample_allspin.channel, ch_allsample_allspin.sample_id, ch_allsample_allspin.spin, tot),
                                                  &outputBDT(allch_allsample_allspin.channel, allch_allsample_allspin.sample_id, allch_allsample_allspin.spin, type_entry.first),
                                                  &outputBDT(allch_allsample_allspin.channel, allch_allsample_allspin.sample_id, allch_allsample_allspin.spin, tot)
                };

                for (const auto& value : current_evaluation){
                    for(auto out : outs)
                        out->Fill(value.first, value.second);
                }
            }
        }
    }

    std::vector<int> CreateMassRange()
    {
        std::set<int> masses;
        for(const auto& sample : samples) {
            if(sample.id.IsSignal()){
                masses.insert(sample.id.mass);
            }
        }
        return std::vector<int>(masses.begin(), masses.end());
    }

    std::vector<std::pair<int,int>> CreateMassSpinRange(const Range<int>& range)
    {
        std::set<std::pair<int,int>> mass_spin;
        for(const SampleEntry& entry : samples)
        {
            if ( entry.id.IsBackground()) continue;
            if ( entry.id.IsSignal() && !range.Contains(entry.id.mass)) continue;
            std::pair<int,int> pair(entry.id.mass, entry.spin);
            mass_spin.insert(pair);
        }

        return std::vector<std::pair<int,int>>(mass_spin.begin(), mass_spin.end());
    }

    void Run()
    {
        std::cout<<"Variabili iniziali: "<<enabled_vars.size()<<std::endl;
        const auto range = mva_setup.mass_range;
        auto mass_range = CreateMassRange();
        std::cout<< "quante masse? "<<mass_range.size() <<std::endl;
        auto mass_spin = CreateMassSpinRange(range);
        std::cout<< "quante coppie massa-spin? "<<mass_spin.size() <<std::endl;

        std::uniform_int_distribution<size_t> it(0, mass_spin.size() - 1);
        std::uniform_int_distribution<size_t> bp(0, benchmarks.size() - 1);
        std::cout<<bkg<<std::endl;

        for (const auto& s : set){
            std::cout << s.channel << s.spin <<std::endl;
            for(const SampleEntry& entry : samples)
            {
                if ( entry.id.IsSignal() && !range.Contains(entry.id.mass) ) continue;
                if ( entry.spin != s.spin) continue;

                std::pair<int,int> pair_check(entry.id.mass, s.spin);
                if ( !entry.id.IsBackground() && std::find(mass_spin.begin(),  mass_spin.end(), pair_check) ==  mass_spin.end())
                    throw exception("Mass and spin couple is not in the possible values set");

                auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
                auto tuple = ntuple::CreateEventTuple(s.channel, input_file.get(), true, ntuple::TreeState::Skimmed);
                auto sumtuple = ntuple::CreateSummaryTuple("summary",input_file.get(), true, ntuple::TreeState::Skimmed);
                auto mergesummary = ntuple::MergeSummaryTuple(*sumtuple.get());

                if(entry.id.IsSM() && args.is_BSM()) {
                    reweight5D->AddFile(*input_file);
                    reweight5D->CreatePdfs();
                }

                Long64_t tot_entries = 0;
                for(Event event : *tuple) {
                    if(tot_entries >= args.number_events()) break;
                    LorentzVectorE_Float bb = event.jets_p4[0] + event.jets_p4[1];
                    if (args.suffix() == "_ANcut"){
                        if (!cuts::hh_bbtautau_2016::hh_tag::m_hh_window().IsInside(event.SVfit_p4.mass(),bb.mass())) continue;
                    }
//                    if (event.p4_1.pt()<40 || event.p4_2.pt()<40) continue;
                    int which_set=0;
                    if(!args.all_data()){
                        if (args.blind()){
                            if (event.split_id >= (mergesummary.n_splits/2)) continue;
                            auto step = (mergesummary.n_splits/2)/args.subdivisions();
                            if (event.split_id >= step*args.which_test() && event.split_id < (step*(args.which_test()+1))) {
                                which_set = 0;
                            }
                            else which_set = 1;
                        }
                        if (!args.blind()){
                            if (event.split_id < (mergesummary.n_splits/2)) continue;
                            auto step = (mergesummary.n_splits/2)/args.subdivisions();
                            if (event.split_id >= (mergesummary.n_splits/2+step*args.which_test()) && event.split_id < (mergesummary.n_splits/2+(step*(args.which_test()+1))) ) {
                                which_set = 0;
                            }
                            else which_set = 1;
                        }
                    }
                    else {
                        if (event.split_id >= (mergesummary.n_splits/2)) which_set = 0;
                        else which_set = 1;
                    }
                    tot_entries++;

                    auto eventInfoPtr =  analysis::MakeEventInfo(Parse<Channel>(s.channel), event) ;
                    EventInfoBase& eventbase = *eventInfoPtr;
                    if (args.suffix() == "_newcut"){
                        if (!cuts::hh_bbtautau_2016::hh_tag::new_m_hh_window().IsInside(eventbase.GetHiggsTTMomentum(false).M(),bb.mass())) continue;
                    }

                    if (!args.is_SM() && args.is_BSM())
                        throw exception("Impossible to reweight events in different benchmark scenario if you don't use SM sample");                    

                    if (entry.id.IsBackground()) {
                        std::pair<int,int> pair_mass_spin;
//                        if (entry.filename == "TT.root") weight_bkg = 831.76/mergesummary.totalShapeWeight; //To Fix
//                        if (entry.filename == "DYJetsToLL_M-50.root") weight_bkg = 5765.4/mergesummary.totalShapeWeight;
                        if (args.is_SM() && args.is_BSM()){
                            size_t which_benchmark = bp(gen2);
                            const auto& benchmark = benchmarks.at_index(which_benchmark);
                            pair_mass_spin = std::make_pair(SampleId::SM().mass, benchmark.point.kl);
                            const SampleId sample_bkg(SampleType::Bkg_TTbar, pair_mass_spin.first);
                            vars->AddEvent(eventbase, sample_bkg, pair_mass_spin.second, entry.weight, which_set);
                        }
                        else {
                            pair_mass_spin = mass_spin.at(it(gen));
                            const SampleId sample_bkg(SampleType::Bkg_TTbar, pair_mass_spin.first);
                            vars->AddEvent(eventbase, sample_bkg, pair_mass_spin.second,entry.weight, which_set);
                        }
                    }
                    else {
                        if (args.is_SM() && args.is_BSM()){
                            size_t which_benchmark = bp(gen2);
                            const auto& benchmark = benchmarks.at_index(which_benchmark);
                            double benchmarkWeight = reweight5D->Get(event, benchmark.point);
                            event.weight_total = eventbase->weight_total/eventbase->weight_bsm_to_sm*benchmarkWeight;
                            vars->AddEvent(eventbase, entry.id, static_cast<int>(benchmark.point.kl), entry.weight , which_set);
                        }
                        else vars->AddEvent(eventbase, entry.id, entry.spin, entry.weight , which_set);
                    }
                }
                std::cout << " channel " << s.channel << "    " << entry.filename << " number of events: " << tot_entries << std::endl;
            }
        }

        if (args.save().size()) {
            vars->SaveEvents(args.save());
            return;
        }
        std::cout<<"upload"<<std::endl;
        vars->UploadEvents();
        std::cout<<"prepare training and test tree"<<std::endl;
        vars->loader->PrepareTrainingAndTestTree( "","", "SplitMode=Random" );

        MvaOptionCollection options = mva_setup.CreateOptionCollection(true);
        const Grid_ND grid(options.GetPositionLimits());
        std::map<std::string, std::string> methods;
        for(const auto& point : grid){
            methods[options.GetName(point)+"_"+std::to_string(args.which_test())] = options.GetConfigString(point);
        }
        std::cout << methods.size() << " metodi" << std::endl;

        auto directory = root_ext::GetDirectory(*outfile.get(), "Evaluation");
        auto directory_roc = root_ext::GetDirectory(*outfile.get(), "ROCCurve");
        auto directory_ks = root_ext::GetDirectory(*outfile.get(), "Kolmogorov");
        auto directory_chi = root_ext::GetDirectory(*outfile.get(), "Chi2");
        auto directory_sb = root_ext::GetDirectory(*outfile.get(), "Significance");

//        auto difference = std::make_shared<BDTData>(outfile.get());

        std::map<std::string, std::shared_ptr<TH1D>> histo_rank_importance;
        std::map<std::string, std::shared_ptr<TH1D>> histo_rank_position;

        MvaTuple mva_tuple(outfile.get(), false);

        for(const auto& point : grid) {
            for (const auto& val : options.GetOptionNames()){
                mva_tuple().param_names.push_back(val.first);
                mva_tuple().param_positions.push_back(val.second);
                mva_tuple().param_values.push_back(options.GetNumericValue(point, val.first));
            }
        }

        int num = 1;
        for(const auto& m : methods){
            std::map<ChannelSampleIdSpin, double> roc_testing;
            std::map<ChannelSampleIdSpin, double> roc_training;
            std::vector<std::pair<std::string, double>> importance;
            std::map<ChannelSampleIdSpin, std::map<size_t, std::vector<double>>> evaluation;
            std::shared_ptr<BDTData> outputBDT;

            mva_tuple().name = m.first;

            std::cout<<"Quale metodo? "<<num<<std::endl;
            auto factory = std::make_shared<TMVA::Factory>("myFactory"+args.output_file(), outfile.get(),"!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification");
            factory->BookMethod(vars->loader.get(), TMVA::Types::kBDT, m.first, m.second);
            std::cout<<"Booked"<<std::endl;

            TrainAllMethods(*factory);
            std::cout<<"Trained"<<std::endl;
            outputBDT = std::make_shared<BDTData>(directory, m.first);

            EvaluateMethod(evaluation, outputBDT->bdt_out, m.first);
            std::cout<<"Evaluated"<<std::endl;

            auto method = dynamic_cast<TMVA::MethodBDT*>(factory->GetMethod(vars->loader->GetName(), m.first));
            for (size_t i = 0; i<vars->names.size(); i++){
                importance.emplace_back(vars->names[i], method->GetVariableImportance(static_cast<UInt_t>(i)));
            }
            for (const auto& ch_spin : set){
                ChannelSampleIdSpin id_sgn_allch_sp{all_channel, mass_tot, ch_spin.spin};
                ChannelSampleIdSpin id_sgn_ch_sp{ch_spin.channel, mass_tot, ch_spin.spin};
                ChannelSampleIdSpin id_bkg_allch_sp{all_channel, bkg, bkg_spin};
                ChannelSampleIdSpin id_bkg_ch_sp{ch_spin.channel, bkg, bkg_spin};

                if (args.range()!="SM"){
                    if (ch_spin.spin == SM_spin) continue;
                }
                else if (ch_spin.spin == SM_spin){
                    id_sgn_allch_sp.sample_id.sampleType = SampleType::Sgn_NonRes;
                    id_sgn_allch_sp.sample_id.mass = 0;
                    id_sgn_ch_sp.sample_id.sampleType = SampleType::Sgn_NonRes;
                    id_sgn_ch_sp.sample_id.mass = 0;
                }
                else continue;

                if (args.is_BSM()){
                    std::cout<<"è BSM"<<std::endl;
                    SampleId sample_sgn(SampleType::Sgn_NonRes, SampleId::SM().mass);
                    SampleId sample_bkg(SampleType::Bkg_TTbar, SampleId::SM().mass);
                    id_sgn_allch_sp.sample_id = sample_sgn;
                    id_sgn_ch_sp.sample_id = sample_sgn;
                    id_bkg_allch_sp.sample_id = sample_bkg;
                    id_bkg_ch_sp.sample_id = sample_bkg;
                    if (!vars->data_pair[0].count(id_sgn_ch_sp)) continue;

                    for(const auto& id : benchmarks){
                        std::cout<<"kl: "<<id.second.point.kl<<std::endl;
                        roc_testing[id_sgn_allch_sp] = method->GetROCIntegral(&outputBDT->bdt_out(id_sgn_allch_sp.channel, id_sgn_allch_sp.sample_id, id.second.point.kl, 0),
                                                                                          &outputBDT->bdt_out(id_bkg_allch_sp.channel, id_bkg_allch_sp.sample_id, id.second.point.kl, 0));
                        roc_training[id_sgn_allch_sp] = method->GetROCIntegral(&outputBDT->bdt_out(id_sgn_allch_sp.channel, id_sgn_allch_sp.sample_id, id.second.point.kl, 1),
                                                                                          &outputBDT->bdt_out(id_bkg_allch_sp.channel, id_bkg_allch_sp.sample_id, id.second.point.kl, 1));
                        std::cout<<"channel: "<< id_sgn_allch_sp.channel<<"  sampleid: "<< id_sgn_allch_sp.sample_id<<" kl: "<<id.second.point.kl <<" ---  ROC testing:  "<<roc_testing[id_sgn_allch_sp]<<std::endl;

                        roc_testing[id_sgn_ch_sp] = method->GetROCIntegral(&outputBDT->bdt_out(id_sgn_ch_sp.channel, id_sgn_ch_sp.sample_id, id.second.point.kl, 0),
                                                                                          &outputBDT->bdt_out(id_bkg_ch_sp.channel, id_bkg_ch_sp.sample_id, id.second.point.kl, 0));
                        roc_training[id_sgn_ch_sp] = method->GetROCIntegral(&outputBDT->bdt_out(id_sgn_ch_sp.channel, id_sgn_ch_sp.sample_id, id_sgn_ch_sp.spin, 1),
                                                                                          &outputBDT->bdt_out(id_bkg_ch_sp.channel, id_bkg_ch_sp.sample_id, id.second.point.kl, 1));
                        std::cout<<"channel: "<< id_sgn_ch_sp.channel<<"  sampleid: "<< id_sgn_ch_sp.sample_id<<" kl: "<<id.second.point.kl <<" ---  ROC testing:  "<<roc_testing[id_sgn_ch_sp]<<std::endl;
                    }
                }

                if (ch_spin.spin == bkg_spin) {
                    id_sgn_ch_sp.spin = spin_tot;
                    id_sgn_allch_sp.spin = spin_tot;
                }
                if (!roc_testing.count(id_sgn_allch_sp)){
                    roc_testing[id_sgn_allch_sp] = method->GetROCIntegral(&outputBDT->bdt_out(id_sgn_allch_sp.channel, id_sgn_allch_sp.sample_id, id_sgn_allch_sp.spin, 0),
                                                                                   &outputBDT->bdt_out(id_bkg_allch_sp.channel, id_bkg_allch_sp.sample_id, id_bkg_allch_sp.spin, 0));
                    roc_training[id_sgn_allch_sp] = method->GetROCIntegral(&outputBDT->bdt_out(id_sgn_allch_sp.channel, id_sgn_allch_sp.sample_id, id_sgn_allch_sp.spin, 1),
                                                                                    &outputBDT->bdt_out(id_bkg_allch_sp.channel, id_bkg_allch_sp.sample_id, id_bkg_allch_sp.spin, 1));
                    std::cout<<"channel: "<< id_sgn_allch_sp.channel<<"  sampleid: "<< id_sgn_allch_sp.sample_id<<" spin: "<<id_sgn_allch_sp.spin <<" ---  ROC testing:  "<<roc_testing[id_sgn_allch_sp]<<std::endl;
                }
                 if (!roc_testing.count(id_sgn_ch_sp)){
                     roc_testing[id_sgn_ch_sp] = method->GetROCIntegral(&outputBDT->bdt_out(id_sgn_ch_sp.channel, id_sgn_ch_sp.sample_id, id_sgn_ch_sp.spin, 0),
                                                                                 &outputBDT->bdt_out(id_bkg_ch_sp.channel, id_bkg_ch_sp.sample_id, id_bkg_ch_sp.spin, 0));
                     roc_training[id_sgn_ch_sp] = method->GetROCIntegral(&outputBDT->bdt_out(id_sgn_ch_sp.channel, id_sgn_ch_sp.sample_id, id_sgn_ch_sp.spin, 1),
                                                                                  &outputBDT->bdt_out(id_bkg_ch_sp.channel, id_bkg_ch_sp.sample_id, id_bkg_ch_sp.spin, 1));
                 }



                if ( ch_spin.spin == SM_spin || ch_spin.spin == bkg_spin) continue;

                for(const auto& sample : mass_range){
                    SampleId sample_sgn(SampleType::Sgn_Res, sample);
                    SampleId sample_bkg(SampleType::Bkg_TTbar, sample);

                    id_sgn_allch_sp.sample_id = sample_sgn;
                    id_sgn_ch_sp.sample_id = sample_sgn;
                    id_bkg_allch_sp.sample_id = sample_bkg;
                    id_bkg_ch_sp.sample_id = sample_bkg;

                    if (!vars->data_pair[0].count(id_sgn_ch_sp)) continue;
                    if (!roc_testing.count(id_sgn_allch_sp)){
                        roc_testing[id_sgn_allch_sp] = method->GetROCIntegral(&outputBDT->bdt_out(id_sgn_allch_sp.channel, id_sgn_allch_sp.sample_id, id_sgn_allch_sp.spin, 0),
                                                                                          &outputBDT->bdt_out(id_bkg_allch_sp.channel, id_bkg_allch_sp.sample_id, id_bkg_allch_sp.spin, 0));
                        roc_training[id_sgn_allch_sp] = method->GetROCIntegral(&outputBDT->bdt_out(id_sgn_allch_sp.channel, id_sgn_allch_sp.sample_id, id_sgn_allch_sp.spin, 1),
                                                                                          &outputBDT->bdt_out(id_bkg_allch_sp.channel, id_bkg_allch_sp.sample_id, id_bkg_allch_sp.spin, 1));
                        std::cout<<"channel: "<< id_sgn_allch_sp.channel<<"  sampleid: "<< id_sgn_allch_sp.sample_id<<" spin: "<<id_sgn_allch_sp.spin <<" ---  ROC testing:  "<<roc_testing[id_sgn_allch_sp]<<std::endl;

                    }
                    if (!roc_testing.count(id_sgn_ch_sp)){
                        roc_testing[id_sgn_ch_sp] = method->GetROCIntegral(&outputBDT->bdt_out(id_sgn_ch_sp.channel, id_sgn_ch_sp.sample_id, id_sgn_ch_sp.spin, 0),
                                                                                          &outputBDT->bdt_out(id_bkg_ch_sp.channel, id_bkg_ch_sp.sample_id, id_bkg_ch_sp.spin, 0));
                        roc_training[id_sgn_ch_sp] = method->GetROCIntegral(&outputBDT->bdt_out(id_sgn_ch_sp.channel, id_sgn_ch_sp.sample_id, id_sgn_ch_sp.spin, 1),
                                                                                          &outputBDT->bdt_out(id_bkg_ch_sp.channel, id_bkg_ch_sp.sample_id, id_bkg_ch_sp.spin, 1));
                    }

                }
            }

            for (const auto& value : roc_testing){
                mva_tuple().roc_testing_value.push_back(value.second);
                mva_tuple().roc_testing_channel.push_back(value.first.channel);
                mva_tuple().err_roc_testing.push_back(0.);
                mva_tuple().roc_testing_mass.push_back(value.first.sample_id.mass);
                mva_tuple().roc_testing_type.push_back(static_cast<int>(value.first.sample_id.sampleType));
                mva_tuple().roc_testing_spin.push_back(value.first.spin);
            }

            for (const auto& value : roc_training){
                mva_tuple().roc_training_value.push_back(value.second);
                mva_tuple().roc_training_channel.push_back(value.first.channel);
                mva_tuple().err_roc_training.push_back(0.);
                mva_tuple().roc_training_mass.push_back(value.first.sample_id.mass);
                mva_tuple().roc_training_type.push_back(static_cast<int>(value.first.sample_id.sampleType));
                mva_tuple().roc_training_spin.push_back(value.first.spin);
            }

            auto directory_roc_method = root_ext::GetDirectory(*directory_roc, m.first);
            std::vector<float> mvaS, mvaB;
            int spin = args.range() == "SM" ? 1 : spin_tot;
            ChannelSampleIdSpin id_sgn{all_channel, mass_tot, spin};
            ChannelSampleIdSpin id_bkg{all_channel, bkg, spin};
            for (auto& eval : evaluation[id_sgn][0])
                mvaS.push_back(static_cast<float>(eval));
            for (auto& eval : evaluation[id_bkg][0])
                mvaB.push_back(static_cast<float>(eval));
            TMVA::ROCCurve roccurve(mvaS, mvaB);
            auto graph = roccurve.GetROCCurve();
            root_ext::WriteObject(*graph, directory_roc_method);
            num++;

            std::cout<<"importance"<<std::endl;
            RankingImportanceVariables(importance, histo_rank_importance);
            std::cout<<"position"<<std::endl;
            auto position = RankingPoisitionVariables(importance, histo_rank_position);

            for (const auto& var: importance){
                mva_tuple().importance.push_back(var.second);
                mva_tuple().var_name.push_back(var.first);
                mva_tuple().position.push_back(position.at(var.first));
            }

            std::map<ChannelSampleIdSpin,double> kolmogorov;
            std::map<ChannelSampleIdSpin,double> chi2;
            std::map<ChannelSampleIdSpin, OptimalSignificance> sign;

            auto directory_ks_method = root_ext::GetDirectory(*directory_ks, m.first);
            auto directory_chi_method = root_ext::GetDirectory(*directory_chi, m.first);
            auto directory_sb_method = root_ext::GetDirectory(*directory_sb, m.first);
            std::cout<<"----"<<m.first<<"----"<<std::endl;
            std::cout<<"Kolmogorov"<<std::endl;
            kolmogorov = KolmogorovTest(evaluation, outputBDT->bdt_out, directory_ks_method, true);
            for (const auto& sample : kolmogorov){
                mva_tuple().KS_value.push_back(sample.second);
                mva_tuple().KS_channel.push_back(sample.first.channel);
                mva_tuple().KS_mass.push_back(sample.first.sample_id.mass);
                mva_tuple().KS_type.push_back(static_cast<int>(sample.first.sample_id.sampleType));
                mva_tuple().KS_spin.push_back(sample.first.spin);
            }

            std::cout<<"Chi"<<std::endl;
            chi2 = ChiSquareTest(evaluation, outputBDT->bdt_out,  directory_chi_method, true);
            for (const auto& sample : chi2){
                mva_tuple().chi_value.push_back(sample.second);
                mva_tuple().chi_channel.push_back(sample.first.channel);
                mva_tuple().chi_mass.push_back(sample.first.sample_id.mass);
                mva_tuple().chi_type.push_back(static_cast<int>(sample.first.sample_id.sampleType));
                mva_tuple().chi_spin.push_back(sample.first.spin);
            }

            std::cout<<"Significance"<<std::endl;
            for (const auto& ch_spin : set){
                if (ch_spin.spin == bkg_spin)
                    sign  = EstimateSignificativity(ch_spin.channel, spin_tot, mass_range, outputBDT->bdt_out, directory_sb_method, true);
                else
                    sign  = EstimateSignificativity(ch_spin.channel, ch_spin.spin, mass_range, outputBDT->bdt_out, directory_sb_method, true);
            }
            for (const auto& entry : sign){
                mva_tuple().optimal_cut.push_back(entry.second.cut);
                mva_tuple().significance.push_back(entry.second.significance.GetValue());
                mva_tuple().significance_err.push_back(entry.second.significance.GetFullError());
                mva_tuple().significance_channel.push_back(entry.first.channel);
                mva_tuple().significance_mass.push_back(entry.first.sample_id.mass);
                mva_tuple().significance_type.push_back(static_cast<int>(entry.first.sample_id.sampleType));
                mva_tuple().significance_spin.push_back(entry.first.spin);
            }
            mva_tuple.Fill();
            outputBDT->bdt_out("all_channel","TT",bkg_spin,1).Scale(1/outputBDT->bdt_out("all_channel","TT",bkg_spin,1).Integral());
            outputBDT->bdt_out("all_channel","TT",bkg_spin,0).Scale(1/outputBDT->bdt_out("all_channel","TT",bkg_spin,0).Integral());
            outputBDT->bdt_out("all_channel","TT",bkg_spin,tot).Scale(1/outputBDT->bdt_out("all_channel","TT",bkg_spin,tot).Integral());

            outputBDT->bdt_out("all_channel","Mtot",3,1).Scale(1/outputBDT->bdt_out("all_channel","Mtot",3,1).Integral());
            outputBDT->bdt_out("all_channel","Mtot",3,0).Scale(1/outputBDT->bdt_out("all_channel","Mtot",3,0).Integral());
            outputBDT->bdt_out("all_channel","Mtot",3,tot).Scale(1/outputBDT->bdt_out("all_channel","Mtot",3,tot).Integral());
        }

        auto directory_ranking = root_ext::GetDirectory(*outfile.get(), "Ranking");
        for (const auto& var: histo_rank_importance)
            root_ext::WriteObject(*var.second, directory_ranking);

        for (const auto& var: histo_rank_position)
            root_ext::WriteObject(*var.second, directory_ranking);

        mva_tuple.Write();


    }
private:
    Arguments args;
    SampleEntryCollection samples;
    MvaSetup mva_setup;
    std::shared_ptr<TFile> outfile;
    std::mt19937_64 gen, gen2;
    MvaVariables::VarNameSet enabled_vars;
    std::shared_ptr<MvaVariablesTMVA> vars;
    std::uniform_int_distribution<uint_fast32_t> test_vs_training;
    std::shared_ptr<NonResHH_EFT::WeightProvider> reweight5D;
    NonResHH_EFT::BenchmarkCollection benchmarks;
};

}
}

PROGRAM_MAIN(analysis::mva_study::MVATraining, Arguments) // definition of the main program function
