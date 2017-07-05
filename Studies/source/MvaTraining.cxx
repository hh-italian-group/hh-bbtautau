/*! Study for Mva Training
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "hh-bbtautau/Analysis/include/MvaConfiguration.h"
#include "hh-bbtautau/Analysis/include/MvaVariables.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
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
#include "hh-bbtautau/Analysis/include/MvaTuple.h"
#include "hh-bbtautau/Analysis/include/MvaVariablesStudy.h"
#include "hh-bbtautau/Analysis/include/MvaMethods.h"


struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(Long64_t, number_events);
    REQ_ARG(std::string, range);
    OPT_ARG(size_t, number_sets, 2);
    OPT_ARG(uint_fast32_t, seed, std::numeric_limits<uint_fast32_t>::max());
    OPT_ARG(uint_fast32_t, seed2, 1234567);
};

namespace analysis {
namespace mva_study{

const SampleId mass_tot = SampleId::MassTot();
const SampleId bkg = SampleId::Bkg();
constexpr int nbin = 220;
constexpr double bin_min = -1.1;
constexpr double bin_max = 1.1;

class MvaVariablesTMVA : public MvaVariables {
public:
    using DataVector = std::vector<double>;
    using DataVectorF = std::vector<float>;
    using MassData = std::map<SampleId, std::deque<DataVector>>;
    static constexpr size_t max_n_vars = 1000;
    DataVector variable;
    DataVectorF variable_float;
    std::map<size_t, MassData> data;
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

    virtual void AddEventVariables(size_t istraining, const SampleId& mass, double /*weight*/) override
    {
        data[istraining][mass].push_back(variable);
    }

    void UploadEvents()
    {
        std::map<SampleId, double> weights;
        for(const auto& entry : data[0]) {
            const SampleId m = entry.first;
            if(m.IsBackground()) weights[m] = 1;
            else weights[m] = 1. / (data[0].at(m).size() + data[1][m].size() );
        }

        for(const auto& data_entry : data) {
            const bool istraining = data_entry.first;
            const TMVA::Types::ETreeType treetype = istraining ? TMVA::Types::kTraining : TMVA::Types::kTesting;
            for(const auto& m_entry : data_entry.second) {
                const SampleId m = m_entry.first;
                const std::string samplename = m.IsSignal() ? "Signal" : "Background";
                for(const auto& vars : m_entry.second) {
                    loader->AddEvent(samplename, treetype, vars, weights.at(m));
                }
            }
        }
    }

    double Evaluate(const std::string& method_name, const std::vector<double>& event)
    {
        if(event.size() != variable.size())
            throw exception("Invalid event size.");
        for(size_t n = 0; n < variable.size(); ++n)
            variable_float[n] = static_cast<float>(event[n]);
        return reader->EvaluateMVA(method_name);
    }

    std::vector<double> EvaluateForAllEvents(const std::string& method_name, size_t istraining, const SampleId& sample)
    {
        if(!data.count(istraining) || !data.at(istraining).count(sample))
            throw exception("Sample not found.");
        const auto& events = data.at(istraining).at(sample);
        std::vector<double> result;
        result.reserve(events.size());
        for(const auto& event : events)
            result.push_back(Evaluate(method_name, event));
        return result;
    }
};


class MVATraining {
public:
    using Event = ::ntuple::Event;
    using EventTuple = ::ntuple::EventTuple;

    MVATraining(const Arguments& _args): args(_args),
        outfile(root_ext::CreateRootFile(args.output_file()+"_"+std::to_string(args.seed())+".root")),
        gen2(args.seed2()), seed_split(0,29)
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

        if(!setups.count(args.range()))
            throw exception("Setups at '%1%' range don't found") %args.range();
        mva_setup = setups.at(args.range());

        enabled_vars.insert(mva_setup.variables.begin(), mva_setup.variables.end());
        if(mva_setup.use_mass_var) {
            enabled_vars.insert("mass");
            enabled_vars.insert("channel");
        }
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

    void EvaluateMethod(std::map<std::string, std::map<SampleId, std::map<size_t, std::vector<double>>>>& evaluation,
                        std::map<std::string, std::map<SampleId, std::shared_ptr<TH1D>>>& outputBDT, const std::string& method_name)
    {
        std::string weightfile = "mydataloader/weights/myFactory_"+method_name+".weights.xml";
        vars->reader->BookMVA(method_name, weightfile);
        outputBDT[method_name][mass_tot] = std::make_shared<TH1D>(("Signal_output_"+method_name).c_str(),("Signal_output_"+method_name).c_str(), nbin, bin_min, bin_max);
        outputBDT[method_name][bkg] = std::make_shared<TH1D>(("Bkg_output_"+method_name).c_str(),("Bkg_output_"+method_name).c_str(), nbin, bin_min, bin_max);
        for(const auto& type_entry : vars->data){
            const auto& sample_entry = type_entry.second;
            auto& evaluation_mass_tot = evaluation[method_name][mass_tot][type_entry.first];
            auto& evaluation_bkg_tot = evaluation[method_name][bkg][type_entry.first];
            for(const auto& entry : sample_entry){
                auto& current_evaluation = evaluation[method_name][entry.first][type_entry.first];
                current_evaluation = vars->EvaluateForAllEvents(method_name, type_entry.first, entry.first);
                if (entry.first.IsSignal()){
                    evaluation_mass_tot.insert(evaluation_mass_tot.end(), current_evaluation.begin(), current_evaluation.end());
                    for (const auto& value : current_evaluation)
                        outputBDT.at(method_name).at(mass_tot)->Fill(value);
                }
                if (entry.first.IsBackground()){
                    evaluation_bkg_tot.insert(evaluation_bkg_tot.end(), current_evaluation.begin(), current_evaluation.end());
                    for (const auto& value : current_evaluation)
                        outputBDT.at(method_name).at(bkg)->Fill(value);
                }
            }
        }
        auto directory = root_ext::GetDirectory(*outfile.get(), "Evaluation");
        auto directory_methody = root_ext::GetDirectory(*directory, (method_name).c_str());
        root_ext::WriteObject(*outputBDT[method_name][mass_tot], directory_methody);
        root_ext::WriteObject(*outputBDT[method_name][bkg], directory_methody);
    }

    std::map<std::string, std::map<SampleId,double>> Kolmogorov(const std::map<std::string, std::map<SampleId, std::map<size_t, std::vector<double>>>>& evaluation) const
    {
        std::map<std::string, std::map<SampleId,double>> kolmogorov;
        std::shared_ptr<TH1D> histo_kolmogorov;
        histo_kolmogorov = std::make_shared<TH1D>("kolmogorov", "kolmogorov", 50, 0, 1.01);
        histo_kolmogorov->SetXTitle("KS");
        auto directory = root_ext::GetDirectory(*outfile.get(), "Kolmogorov");
        for(const auto& method : evaluation){
            for (const auto& sample : method.second){
                std::map<size_t, std::vector<double>> ks_vector;
                for (auto tvt : sample.second){
                    auto type = tvt.first;
                    std::sort(tvt.second.begin(), tvt.second.end());
                    ks_vector[type] = std::move(tvt.second);
                }
                double ks = TMath::KolmogorovTest(static_cast<int>(ks_vector.at(0).size()), ks_vector.at(0).data(),
                                                  static_cast<int>(ks_vector.at(1).size()), ks_vector.at(1).data(), "");
                kolmogorov[method.first][sample.first] = ks;
                histo_kolmogorov->Fill(ks);
            }
        }
        root_ext::WriteObject(*histo_kolmogorov, directory);
        return kolmogorov;
    }

    std::map<std::string, std::pair<double, double>> EstimateSignificativity(const std::map<std::string, std::map<SampleId, std::shared_ptr<TH1D>>>& outputBDT)
    {
        std::map<std::string, std::pair<double, double>> sign;
        auto directory = root_ext::GetDirectory(*outfile.get(), "Significance");
        for(const auto& method : outputBDT){
            std::vector<std::pair<double,double>> cuts;
            auto histo_sign = CreatePlot((method.first+"_significance").c_str(),(method.first+"_significance").c_str(),"output BDT","S/(sqrt(B))");
            for(int i = 0; i<=nbin; ++i){
                auto output = method.second.at(mass_tot)->GetBinCenter(i);
                auto s = method.second.at(mass_tot)->Integral(i, nbin+1);
                auto b = method.second.at(bkg)->Integral(i, nbin+1);
                double significance = 0;
                if (b != 0) significance = s/pow(b, 0.5);
                cuts.emplace_back(output, significance);
                histo_sign->SetPoint(i, output, significance);
            }
            std::sort(cuts.begin(), cuts.end(), [](std::pair<double, double> el1, std::pair<double, double> el2){
                return el1.second > el2.second;
            });
            sign[method.first] = cuts.front();
            root_ext::WriteObject(*histo_sign, directory);
        }
        return sign;
    }

    std::vector<int> CreateMassRange() const
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
        const auto mass_range = CreateMassRange();
        std::cout<< mass_range.size() <<std::endl;
        std::uniform_int_distribution<size_t> it(0, mass_range.size() - 1);

        for(size_t j = 0; j<mva_setup.channels.size(); j++){
            std::cout << mva_setup.channels[j] << std::endl;
            for(const SampleEntry& entry : samples)
            {
                if ( entry.id.IsSignal() && !range.Contains(entry.id.mass) ) continue;
                if ( entry.id.IsBackground() && Parse<Channel>(entry.channel) != mva_setup.channels.at(j) )
                    continue;
                auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
                EventTuple tuple(ToString(mva_setup.channels[j]), input_file.get(), true, {} , GetMvaBranches());
                Long64_t tot_entries = 0;
                for(Long64_t current_entry = 0; tot_entries < args.number_events() && current_entry < tuple.GetEntries(); ++current_entry) {
                    uint_fast32_t seed = seed_split(gen2);
                    if (seed>15) continue;
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
                    if (entry.id.IsBackground()) {
                        const int mass_background = mass_range.at(it(gen));
                        const SampleId sample_bkg(SampleType::Bkg_TTbar, mass_background);
                        vars->AddEvent(event, sample_bkg, entry.weight);
                    }
                    else vars->AddEvent(event, entry.id, entry.weight);
                }
                std::cout << " channel " << mva_setup.channels[j] << "    " << entry.filename << " number of events: " << tot_entries << std::endl;
            }
        }
        vars->UploadEvents();
        vars->loader->PrepareTrainingAndTestTree( "","", "SplitMode=Random" );

        MvaOptionCollection options = mva_setup.CreateOptionCollection(true);
        const Grid_ND grid(options.GetPositionLimits());
        std::map<std::string, std::string> methods;
        for(const auto& point : grid){
            methods[options.GetName(point)+"_"+std::to_string(args.seed())] = options.GetConfigString(point);
        }
        std::cout << methods.size() << " metodi" << std::endl;

        std::map<std::string, double> ROCintegral;
        std::map<std::string, std::vector<std::pair<std::string, double>>> importance;
        std::map<std::string, std::map<SampleId, std::map<size_t, std::vector<double>>>> evaluation;
        std::map<std::string, std::map<SampleId, std::shared_ptr<TH1D>>> outputBDT;
        for(const auto& m : methods){
            auto factory = std::make_shared<TMVA::Factory>("myFactory", outfile.get(),"!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification");
            factory->BookMethod(vars->loader.get(), TMVA::Types::kBDT, m.first, m.second);
            TrainAllMethods(*factory);
            EvaluateMethod(evaluation, outputBDT, m.first);
            auto method = dynamic_cast<TMVA::MethodBDT*>(factory->GetMethod(vars->loader->GetName(), m.first));
            for (size_t i = 0; i<vars->names.size(); i++){
                importance[m.first].emplace_back(vars->names[i], method->GetVariableImportance(static_cast<UInt_t>(i)));
            }
            ROCintegral[m.first] = method->GetROCIntegral(outputBDT.at(m.first).at(mass_tot).get(), outputBDT.at(m.first).at(bkg).get());
        }



        std::cout<<"importance"<<std::endl;
        RankingImportanceVariables(importance);
        std::cout<<"position"<<std::endl;
        auto position = RankingPoisitionVariables(importance);
        std::cout<<"kolmogorov"<<std::endl;
        std::map<std::string, std::map<SampleId,double>> kolmogorov = Kolmogorov(evaluation);
        std::cout<<"Significativity"<<std::endl;
        std::map<std::string, std::pair<double, double>> sign = EstimateSignificativity(outputBDT);

        ntuple::MvaTuple mva_tuple(outfile.get(), false);
        std::cout<<"options"<<std::endl;
        for(const auto& point : grid) {
            auto name = options.GetName(point)+"_"+std::to_string(args.seed());
//            mva_tuple().name = name;
            mva_tuple().NTrees = static_cast<UInt_t>(options.GetNumericValue(point, "NTrees"));
            mva_tuple().shrinkage = options.GetNumericValue(point, "shrinkage");
            mva_tuple().BaggedSampleFraction = options.GetNumericValue(point, "BaggedSampleFraction");
            mva_tuple().MaxDepth = options.GetNumericValue(point, "MaxDepth");
            mva_tuple().MinNodeSize = options.GetNumericValue(point, "MinNodeSize");
            mva_tuple().cut = sign[name].first;
            mva_tuple().significance = sign[name].second;
            mva_tuple().ROCIntegral = ROCintegral[name];
            for (const auto& sample : kolmogorov.at(name)){
                mva_tuple().KS_mass.push_back(sample.first.mass);
                mva_tuple().KS_value.push_back(sample.second);
                if (sample.first.IsSignal())
                    mva_tuple().KS_type.push_back(1);
                else if (sample.first.IsBackground())
                    mva_tuple().KS_type.push_back(-1);
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
    std::uniform_int_distribution<uint_fast32_t> test_vs_training;
    std::uniform_int_distribution<uint_fast32_t> seed_split;
};

}
}

PROGRAM_MAIN(analysis::mva_study::MVATraining, Arguments) // definition of the main program function
