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


struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(Long64_t, number_events);
    REQ_ARG(std::string, range);
    REQ_ARG(Long64_t, number_variables);
    REQ_ARG(int, which_test);
    OPT_ARG(size_t, number_sets, 2);
    OPT_ARG(uint_fast32_t, seed, std::numeric_limits<uint_fast32_t>::max());
    OPT_ARG(uint_fast32_t, seed2, 1234567);
    OPT_ARG(std::string, channel, "");
    OPT_ARG(std::string, save, "");
};

namespace analysis {
namespace mva_study{

class MvaVariablesTMVA : public MvaVariables {
public:
    using DataVector = std::vector<double>;
    using DataVectorF = std::vector<float>;
    using MassData_pair = std::map<SampleId, std::deque<std::pair<DataVector, double>>>;
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

    virtual void AddEventVariables(size_t istraining, const SampleId& mass, double weight) override
    {
        data_pair[istraining][mass].emplace_back(variable, weight);
    }

    void UploadEvents()
    {
        std::map<SampleId, double> weights;
        for(const auto& entry : data_pair[0]) {
            const SampleId m = entry.first;
            double tot_weight = 0;
            for (const auto& val : data_pair){
                for (const auto& x : val.second.at(m)){
                    tot_weight+= x.second;
                }
            }
            weights[m] = 1. / (tot_weight);
            std::cout<<m<<" "<<weights[m]<<std::endl;
        }

        for(const auto& data_entry : data_pair) {
            const bool istraining = data_entry.first;
            const TMVA::Types::ETreeType treetype = istraining ? TMVA::Types::kTraining : TMVA::Types::kTesting;
            for(const auto& m_entry : data_entry.second) {
                const SampleId m = m_entry.first;
                const std::string samplename = m.IsSignal() ? "Signal" : "Background";
                for(const auto& vars : m_entry.second) {
                    loader->AddEvent(samplename, treetype, vars.first , weights.at(m));
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
                const SampleId m = m_entry.first;
                for(const auto& vars : m_entry.second) {
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

    std::vector<std::pair<double,double>> EvaluateForAllEvents(const std::string& method_name, size_t istraining, const SampleId& sample)
    {
        if(!data_pair.count(istraining) || !data_pair.at(istraining).count(sample))
            throw exception("Sample not found.");
        const auto& events = data_pair.at(istraining).at(sample);
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

    MVATraining(const Arguments& _args): args(_args),
        outfile(root_ext::CreateRootFile(args.output_file()+"_"+std::to_string(args.which_test())+".root")), gen2(args.seed2()),  seed_split(0,29)
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
        }

        if (args.channel().size()){
            channels.push_back(args.channel());
        }
        else
            for (const auto& channel : mva_setup.channels)
                channels.push_back(ToString(channel));

        std::cout<<enabled_vars.size()<<std::endl;
        std::mt19937_64 seed_gen(args.seed());
        std::uniform_int_distribution<uint_fast32_t> seed_distr(100000, std::numeric_limits<uint_fast32_t>::max());
        uint_fast32_t seed = seed_distr(seed_gen);
        gen.seed(seed_distr(seed_gen));
        vars = std::make_shared<MvaVariablesTMVA>(args.number_sets(), seed, enabled_vars);
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

    void RankingImportanceVariables(const std::map<std::string, std::vector<std::pair<std::string, double>>>& importance) const
    {
        std::map<std::string, std::shared_ptr<TH1D>> histo_rank;
        auto directory = root_ext::GetDirectory(*outfile.get(), "Ranking");
        for (const auto& method: importance){
            for (const auto& pair: method.second){
                if (!histo_rank.count(pair.first)){
                    histo_rank[pair.first] = std::make_shared<TH1D>((pair.first+"_importance").c_str(),(pair.first+"_importance").c_str(),100, 0, 0.1);
                    histo_rank.at(pair.first)->SetCanExtend(TH1::kAllAxes);
                    histo_rank.at(pair.first)->SetXTitle("importance");
                }
                histo_rank.at(pair.first)->Fill(pair.second);
            }
        }
        for (const auto& var: histo_rank)
            root_ext::WriteObject(*var.second, directory);
    }

    std::map<std::string, std::map<std::string, size_t>> RankingPoisitionVariables(const std::map<std::string, std::vector<std::pair<std::string, double>>>& importance) const
    {
        std::map<std::string, std::map<std::string, size_t>> position;
        for(auto method: importance){
            std::sort(method.second.begin(), method.second.end(), [](auto el1, auto el2){
                return el1.second > el2.second;
            } );
            for (size_t i = 1; i <= method.second.size(); i++){
                position[method.first][method.second[i-1].first] = i;
            }
        }
        std::map<std::string, std::shared_ptr<TH1D>> histo_rank;
        auto directory = root_ext::GetDirectory(*outfile.get(), "Ranking");
        for (const auto& method: position){
            for (const auto& pair: method.second){
               if (!histo_rank.count(pair.first)){
                   histo_rank[pair.first] = std::make_shared<TH1D>((pair.first+"_position").c_str(),(pair.first+"_position").c_str(), method.second.size(), 0.5, method.second.size() + 0.5);
                   histo_rank[pair.first]->SetXTitle("position");
               }
                histo_rank[pair.first]->Fill(pair.second);
            }
        }
        for (const auto& var: histo_rank)
            root_ext::WriteObject(*var.second, directory);

        return position;
    }

    void EvaluateMethod(std::map<SampleId, std::map<size_t, std::vector<double>>>& evaluation, BDTData::Entry& outputBDT, const std::string& method_name)
    {
        std::string weightfile = "mydataloader/weights/myFactory"+args.output_file()+"_"+method_name+".weights.xml";
        vars->reader->BookMVA(method_name, weightfile);
        for(const auto& type_entry : vars->data_pair){
            const auto& sample_entry = type_entry.second;
            std::cout<<type_entry.first<<std::endl;
            for(const auto& entry : sample_entry){
                auto current_evaluation = vars->EvaluateForAllEvents(method_name, type_entry.first, entry.first);
                for (const auto& val: current_evaluation)
                    evaluation[entry.first][type_entry.first].push_back(val.first);

                const SampleId tot_sample = entry.first.IsSignal() ? mass_tot : bkg;
                auto& eval = evaluation[tot_sample][type_entry.first];
                eval.insert(eval.end(), evaluation[entry.first][type_entry.first].begin(), evaluation[entry.first][type_entry.first].end());

                std::vector<BDTData::Hist*> outs = { &outputBDT(tot_sample, tot), &outputBDT(tot_sample, type_entry.first),
                                                     &outputBDT(entry.first, tot), &outputBDT(entry.first, type_entry.first) };
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
            if(sample.id.IsSignal())
                masses.insert(sample.id.mass);
        }
        return std::vector<int>(masses.begin(), masses.end());
    }

    void Run()
    {
        std::cout<<"Variabili iniziali: "<<enabled_vars.size()<<std::endl;
        const auto range = mva_setup.mass_range;
        auto mass_range = CreateMassRange();
        std::cout<< mass_range.size() <<std::endl;
        std::uniform_int_distribution<size_t> it(0, mass_range.size() - 1);
        std::uniform_int_distribution<int> split(0, 3);
        std::cout<<bkg<<std::endl;

        for(size_t j = 0; j<channels.size(); j++){
            std::cout << channels[j] << std::endl;
            for(const SampleEntry& entry : samples)
            {
                if ( entry.id.IsSignal() && !range.Contains(entry.id.mass) ) continue;
                if ( entry.id.IsBackground() && entry.channel != channels.at(j) )
                    continue;
                auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
                auto tuple = ntuple::CreateEventTuple(ToString(channels[j]), input_file.get(), true, ntuple::TreeState::Skimmed);
                Long64_t tot_entries = 0;
                for(const Event& event : *tuple) {
                    if(tot_entries >= args.number_events()) break;
                    if (event.split_id >= 15) continue;
                    tot_entries++;
                    int set = split(gen2);
                    int which_set = set == args.which_test() ? 0 : 1;
                    if (entry.id.IsBackground()) {
                        const int mass_background = mass_range.at(it(gen));
                        const SampleId sample_bkg(SampleType::Bkg_TTbar, mass_background);
                        vars->AddEvent(event, sample_bkg, entry.weight, which_set);
                    }
                    else vars->AddEvent(event, entry.id, entry.weight, which_set);
                }
                std::cout << " channel " << channels[j] << "    " << entry.filename << " number of events: " << tot_entries << std::endl;
            }
        }
        if (args.save().size()) {
            vars->SaveEvents(args.save());
            return;
        }
        vars->UploadEvents();
        vars->loader->PrepareTrainingAndTestTree( "","", "SplitMode=Random" );

        MvaOptionCollection options = mva_setup.CreateOptionCollection(true);
        const Grid_ND grid(options.GetPositionLimits());
        std::map<std::string, std::string> methods;
        for(const auto& point : grid){
            methods[options.GetName(point)+"_"+std::to_string(args.which_test())] = options.GetConfigString(point);
        }
        std::cout << methods.size() << " metodi" << std::endl;

        std::map<std::string, double> ROCintegral;
        std::map<std::string, std::map<int, double>> roc;
        std::map<std::string, std::vector<std::pair<std::string, double>>> importance;
        std::map<std::string, std::map<SampleId, std::map<size_t, std::vector<double>>>> evaluation;
        std::map<std::string, std::shared_ptr<BDTData>> outputBDT;
        auto directory = root_ext::GetDirectory(*outfile.get(), "Evaluation");
        auto directory_roc = root_ext::GetDirectory(*outfile.get(), "ROCCurve");
        for(const auto& m : methods){
            auto factory = std::make_shared<TMVA::Factory>("myFactory"+args.output_file(), outfile.get(),"!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification");
            factory->BookMethod(vars->loader.get(), TMVA::Types::kBDT, m.first, m.second);
            std::cout<<"Booked"<<std::endl;
            TrainAllMethods(*factory);
            std::cout<<"Trained"<<std::endl;
            outputBDT[m.first] = std::make_shared<BDTData>(directory, m.first);
            EvaluateMethod(evaluation[m.first], outputBDT[m.first]->bdt_out, m.first);
            std::cout<<"Evaluated"<<std::endl;
            auto method = dynamic_cast<TMVA::MethodBDT*>(factory->GetMethod(vars->loader->GetName(), m.first));
            for (size_t i = 0; i<vars->names.size(); i++){
                importance[m.first].emplace_back(vars->names[i], method->GetVariableImportance(static_cast<UInt_t>(i)));
            }
            ROCintegral[m.first] = method->GetROCIntegral(&outputBDT.at(m.first)->bdt_out(mass_tot, 0), &outputBDT.at(m.first)->bdt_out(bkg, 0));
            std::cout<<"ROC "<<ROCintegral[m.first]<<std::endl;
            for(const auto& sample : mass_range){
                SampleId sample_sgn(SampleType::Sgn_Res, sample);
                SampleId sample_bkg(SampleType::Bkg_TTbar, sample);
                roc[m.first][sample] = method->GetROCIntegral(&outputBDT.at(m.first)->bdt_out(sample_sgn, 0), &outputBDT.at(m.first)->bdt_out(sample_bkg, 0));
                std::cout<<sample<<"    "<<roc[m.first][sample]<<std::endl;
            }
            auto directory_roc_method = root_ext::GetDirectory(*directory_roc, m.first);
            std::vector<float> mvaS, mvaB;
            for (auto& eval : evaluation[m.first][mass_tot][0])
                mvaS.push_back(static_cast<float>(eval));
            for (auto& eval : evaluation[m.first][bkg][0])
                mvaB.push_back(static_cast<float>(eval));

            TMVA::ROCCurve roccurve(mvaS, mvaB);
            auto graph = roccurve.GetROCCurve();
            root_ext::WriteObject(*graph, directory_roc_method);
        }
        std::cout<<"importance"<<std::endl;
        RankingImportanceVariables(importance);
        std::cout<<"position"<<std::endl;
        auto position = RankingPoisitionVariables(importance);
        std::map<std::string, std::map<SampleId,double>> kolmogorov;
        std::map<std::string, std::map<SampleId,double>> chi2;
        std::map<std::string, std::map<int, std::pair<double, PhysicalValue>>> sign;
        auto directory_ks = root_ext::GetDirectory(*outfile.get(), "Kolmogorov");
        auto directory_chi = root_ext::GetDirectory(*outfile.get(), "Chi2");
        auto directory_sb = root_ext::GetDirectory(*outfile.get(), "Significance");
        auto difference = std::make_shared<BDTData>(outfile.get());
        for (const auto& m: methods){
            auto directory_ks_method = root_ext::GetDirectory(*directory_ks, m.first);
            auto directory_chi_method = root_ext::GetDirectory(*directory_chi, m.first);
            auto directory_sb_method = root_ext::GetDirectory(*directory_sb, m.first);
            std::cout<<"----"<<m.first<<"----"<<std::endl;
            std::cout<<"Kolmogorov"<<std::endl;
            kolmogorov[m.first] = Kolmogorov(evaluation[m.first], outputBDT.at(m.first)->bdt_out, difference->difference, directory_ks_method);
            std::cout<<"Chi"<<std::endl;
            chi2[m.first] = ChiSquare(evaluation[m.first], outputBDT.at(m.first)->bdt_out,  directory_chi_method);
            std::cout<<"Significance"<<std::endl;
            sign[m.first]  = EstimateSignificativity(mass_range, outputBDT.at(m.first)->bdt_out, directory_sb_method, true);
        }

        MvaTuple mva_tuple(outfile.get(), false);
        std::cout<<"options"<<std::endl;
        for(const auto& point : grid) {
            auto name = options.GetName(point)+"_"+std::to_string(args.which_test());
            mva_tuple().name = name;

            for (const auto val : options.GetOptionNames()){
                mva_tuple().param_names.push_back(val.first);
                mva_tuple().param_positions.push_back(val.second);
                mva_tuple().param_values.push_back(options.GetNumericValue(point, val.first));
            }

            for (const auto entry : sign[name]){
                mva_tuple().optimal_cut.push_back(entry.second.first);
                mva_tuple().significance.push_back(entry.second.second.GetValue());
                mva_tuple().significance_err.push_back(entry.second.second.GetFullError());
                mva_tuple().significance_mass.push_back(entry.first);
            }

            for (const auto& sample : kolmogorov.at(name)){
                mva_tuple().KS_mass.push_back(sample.first.mass);
                mva_tuple().KS_value.push_back(sample.second);
                mva_tuple().KS_type.push_back(static_cast<int>(sample.first.sampleType));
            }

            for (const auto& sample : chi2.at(name)){
                mva_tuple().chi_mass.push_back(sample.first.mass);
                mva_tuple().chi_value.push_back(sample.second);
                mva_tuple().chi_type.push_back(static_cast<int>(sample.first.sampleType));
            }

            mva_tuple().ROCIntegral = ROCintegral[name];
            for (const auto& value : roc[name]){
                mva_tuple().roc_value.push_back(value.second);
                mva_tuple().roc_mass.push_back(value.first);
            }

            for (const auto& var: importance.at(name)){
                mva_tuple().importance.push_back(var.second);
                mva_tuple().var_name.push_back(var.first);
                mva_tuple().position.push_back(position.at(name).at(var.first));
            }
            mva_tuple.Fill();
        }
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
    std::uniform_int_distribution<uint_fast32_t> test_vs_training, seed_split;
    std::vector<std::string> channels;
};

}
}

PROGRAM_MAIN(analysis::mva_study::MVATraining, Arguments) // definition of the main program function
