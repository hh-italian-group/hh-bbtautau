/*! Study for Mva Training
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "hh-bbtautau/Analysis/include/Lester_mt2_bisect.h"
#include "hh-bbtautau/Analysis/include/MvaConfiguration.h"
#include "hh-bbtautau/Analysis/include/MvaVariables.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
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
#include "TMath.h"
#include <fstream>
#include <random>
#include "hh-bbtautau/Analysis/include/MvaConfigurationReader.h"

#include "AnalysisTools/Core/include/SmartTree.h"

#define MVA_DATA() \
    VAR(UInt_t, NTrees) \
    VAR(double, shrinkage) \
    VAR(double, BaggedSampleFraction) \
    VAR(double, MaxDepth) \
    VAR(double, MinNodeSize) \
    VAR(double, ROCIntegral) \
    VAR(std::vector<double>, KS_value) \
    VAR(std::vector<double>, KS_mass) \
    VAR(std::vector<double>, position) \
    VAR(std::vector<double>, importance) \
    VAR(std::vector<std::string>, var_name) \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(ntuple, MvaResults, MvaTuple, MVA_DATA, "mva_result")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(ntuple, MvaTuple, MVA_DATA)
#undef VAR
#undef MVA_DATA

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(Long64_t, number_events);
    OPT_ARG(size_t, number_sets, 2);
    OPT_ARG(uint_fast32_t, seed, std::numeric_limits<uint_fast32_t>::max());
};

namespace analysis {
namespace mva_study{

const SampleId mass_tot(SampleType::Sgn_Res, 2000);
const SampleId bkg(SampleType::Bkg_TTbar, 0);
class MvaVariablesTMVA : public MvaVariables {
public:
    using DataVector = std::vector<double>;
    using DataVectorF = std::vector<float>;
    using MassData = std::map<SampleId, std::deque<DataVector>>;
    DataVector variable;
    DataVectorF variable_float;
    std::map<size_t, MassData> data;
    std::shared_ptr<TMVA::DataLoader> loader;
    std::shared_ptr<TMVA::Reader> reader;
    MvaVariablesTMVA(size_t _number_set = 2, uint_fast32_t _seed = std::numeric_limits<uint_fast32_t>::max(),
                            const VarNameSet& _enabled_vars = {}) :
        MvaVariables(_number_set, _seed, _enabled_vars), loader(new TMVA::DataLoader("mydataloader")), reader(new TMVA::Reader) {}

    virtual void SetValue(const std::string& name, double value, char type = 'F') override
    {
        if (!names.count(name)){
            variable.push_back(0);
            names[name] = variable.size()-1;
            loader->AddVariable(name, type);
            variable_float.resize(variable.size());
            reader->AddVariable(name, &variable_float.at(names.at(name)));
        }
        variable.at(names.at(name)) = value;

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

void TrainAllMethods(std::shared_ptr<TMVA::Factory> factory)
{
   std::map< TString, TMVA::Factory::MVector*> fMethodsMap = factory->fMethodsMap;
   std::map<TString,TMVA::Factory::MVector*>::iterator itrMap;
   for(itrMap = fMethodsMap.begin(); itrMap != fMethodsMap.end(); itrMap++)
   {
      TMVA::Factory::MVector *methods=itrMap->second;
      TMVA::Factory::MVector::iterator itrMethod;

      for( itrMethod = methods->begin(); itrMethod != methods->end(); itrMethod++ ) {
          TMVA::Event::SetIsTraining(kTRUE);
          TMVA::MethodBase* mva = dynamic_cast<TMVA::MethodBase*>(*itrMethod);
          if(mva==0) continue;
          if (mva->DataInfo().GetNClasses() < 2 )
              std::cout << TMVA::kFATAL << "You want to do classification training, but specified less than two classes." << std::endl;
          if (mva->Data()->GetNTrainingEvents() < 10) {
            std::cout << TMVA::kWARNING << "Method " << mva->GetMethodName() << " not trained (training tree has less entries ["
             << mva->Data()->GetNTrainingEvents() << "] than required [10]" << std::endl;
            continue;
          }
          mva->TrainMethod();
      }
   }
}

class MVATraining {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    MVATraining(const Arguments& _args): args(_args),
        outfile(root_ext::CreateRootFile(args.output_file()+"_"+std::to_string(args.seed())+".root")),
        factory(new TMVA::Factory ("myFactory", outfile.get(),"!V:!Silent:Color:DrawProgressBar:Transformations=I,D:AnalysisType=Classification")),
        gen(args.seed()), test_vs_training(100000, std::numeric_limits<uint_fast32_t>::max())
    {
        MvaSetupCollection setups;
        SampleEntryListCollection samples_list;

        ConfigReader configReader;
        MvaConfigReader setupReader(setups);
        configReader.AddEntryReader("SETUP", setupReader, true);
        SampleConfigReader sampleReader(samples_list);
        configReader.AddEntryReader("FILES", sampleReader, false);
        configReader.ReadConfig(args.cfg_file());

        samples = samples_list.at("inputs").files;
        mva_setup = setups.at("LowMassRange");

        enabled_vars.insert(mva_setup.variables.begin(), mva_setup.variables.end());
        if(mva_setup.use_mass_var) {
            enabled_vars.insert("mass");
            enabled_vars.insert("channel");
        }
        uint_fast32_t seed = test_vs_training(gen);
        vars = std::make_shared<MvaVariablesTMVA>(args.number_sets(), seed, enabled_vars);
    }

    std::map<std::string, std::vector<std::pair<std::string, double>>> RankingImportanceVariables(Grid_ND grid, MvaOptionCollection options){
        std::map<std::string, std::vector<std::pair<std::string, double>>> importance;

        for(const auto& point : grid) {
            auto name = options.GetName(point)+"_"+std::to_string(args.seed());
            auto method = dynamic_cast<TMVA::MethodBDT*>(factory->GetMethod(vars->loader->GetName(), name));
            if(!method)
                throw exception("method '%1%' not found.") % name;
            if (mva_setup.use_mass_var){
//                importance[name].emplace_back("channel", method->GetVariableImportance(0));
                importance[name].emplace_back("mass", method->GetVariableImportance(1));
            }
            for (size_t i = 0; i<mva_setup.variables.size(); i++){
                importance[name].emplace_back(mva_setup.variables[i], method->GetVariableImportance(static_cast<UInt_t>(i)+2));
            }
        }

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

        return importance;
    }

    std::map<std::string, std::map<std::string, double>> RankingPoisitionVariables(std::map<std::string, std::vector<std::pair<std::string, double>>> importance){
        std::map<std::string, std::map<std::string, double>> position;

        for(auto method: importance){
            std::sort(method.second.begin(), method.second.end(), [](auto el1, auto el2){
                return el1.second > el2.second;
            } );
            for (size_t i = 0; i < method.second.size(); i++){
                position[method.first][method.second[i].first] = i;
            }
        }

        std::map<std::string, std::shared_ptr<TH1D>> histo_rank;
        auto directory = root_ext::GetDirectory(*outfile.get(), "Ranking");
        for (const auto& method: position){
            for (const auto& pair: method.second){
               if (!histo_rank.count(pair.first)){
                   histo_rank[pair.first] = std::make_shared<TH1D>((pair.first+"_position").c_str(),(pair.first+"_position").c_str(), method.second.size(), 0, method.second.size());
                    histo_rank[pair.first]->SetXTitle("position");
               }
                histo_rank[pair.first]->Fill(pair.second);
            }
        }

        for (const auto& var: histo_rank)
            root_ext::WriteObject(*var.second, directory);

        return position;
    }

    std::map<std::string, std::map<SampleId,double>> Kolmogorov(std::map<std::string, std::map<SampleId, std::map<size_t, std::vector<double>>>> evaluation)
    {
        std::map<std::string, std::map<SampleId,double>> kolmogorov;
        for(const auto& method : evaluation){
            for (const auto& sample : method.second){
                std::map<int, double*> ks_vector;
                std::map<int, size_t> dim;
                for (auto tvt : sample.second){
                    auto type = static_cast<int>(tvt.first);
                    std::sort(tvt.second.begin(), tvt.second.end());
                    std::cout<<tvt.second.front()<<"    "<<sample.first.mass<<"    "<<tvt.second.size()<<" "<<tvt.second.back()<<std::endl;
                    dim[type] = tvt.second.size();
                    ks_vector[type] = tvt.second.data();
                }
                double ks = TMath::KolmogorovTest(static_cast<int>(dim.at(0)), ks_vector.at(0), static_cast<int>(dim.at(1)), ks_vector.at(1), "M");
                kolmogorov[method.first][sample.first] = ks;
            }
        }
        return kolmogorov;
    }

    void Run()
    {
        std::cout<<"Variabili iniziali: "<<enabled_vars.size()<<std::endl;
        std::vector<int> vec_mass;
        auto range = mva_setup.mass_range;
        size_t begin = 0, end = 0;
        for (size_t i=0; i<samples.size(); i++){
            vec_mass.push_back(samples[i].id.mass);
            if (samples[i].id.mass == range.min()) begin = i;
            if (samples[i].id.mass == range.max()) end = i;
        }

        for(size_t j = 0; j<mva_setup.channels.size(); j++){

            std::cout << ToString(mva_setup.channels[j]) << std::endl;

            for(const SampleEntry& entry : samples)
            {
                if ( entry.id.IsSignal() && !range.Contains(entry.id.mass) ) continue;

                if ( entry.id.IsBackground() && Parse<Channel>(entry.channel) != mva_setup.channels.at(j) )
                    continue;

                auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
                EventTuple tuple(ToString(mva_setup.channels[j]), input_file.get(), true, {} , GetMvaBranches());

                Long64_t tot_entries = 0, current_entry = 0;
                int mass_background = 0;
                SampleId sample_bkg;

                if (entry.id.IsBackground()) {
                    std::uniform_int_distribution<size_t> it(begin, end);
                    mass_background = vec_mass[static_cast<size_t>(it(gen))];
                    std::cout<<mass_background<<std::endl;
                    sample_bkg.sampleType = SampleType::Bkg_TTbar;
                    sample_bkg.mass = mass_background;
                }

                while(tot_entries < args.number_events() && current_entry < tuple.GetEntries()) {
                    tuple.GetEntry(current_entry);
                    const Event& event = tuple.data();
                    if (event.eventEnergyScale != 0 || (event.q_1+event.q_2) != 0 || event.jets_p4.size() < 2
                        || event.extraelec_veto == true || event.extramuon_veto == true || event.jets_p4[0].eta() > cuts::btag_2016::eta
                        || event.jets_p4[1].eta() > cuts::btag_2016::eta){
                        current_entry++;
                        continue;
                    }
                    LorentzVectorE_Float bb = event.jets_p4[0] + event.jets_p4[1];
                    double ellipse_cut = pow(event.SVfit_p4.mass()-116,2)/pow(35.,2) + pow(bb.mass()-111,2)/pow(45.,2);
                    if (ellipse_cut>1){
                        current_entry++;
                        continue;
                    }

                    current_entry++;
                    tot_entries++;

                    if (entry.id.IsBackground()) {
                        vars->AddEvent(event, sample_bkg, entry.weight);
                    }
                    else vars->AddEvent(event, entry.id, entry.weight);
                }
                std::cout << " channel " << ToString(mva_setup.channels[j]) << "    " << entry << " number of events: " << tot_entries << std::endl;
            }
        }
        vars->UploadEvents();
        vars->loader->PrepareTrainingAndTestTree( "","", "SplitMode=Random" );
        MvaOptionCollection options = mva_setup.CreateOptionCollection(true);
        const Grid_ND grid(options.GetPositionLimits());
        std::cout<<"grid"<<std::endl;
        std::map<std::string, std::string> methods;
        for(const auto& point : grid){
            methods[options.GetName(point)+"_"+std::to_string(args.seed())] = options.GetConfigString(point);
            std::cout << options.GetName(point)+"_"+std::to_string(args.seed()) << "    " << options.GetConfigString(point) << std::endl;
        }

        for(const auto& m : methods)
            factory->BookMethod(vars->loader.get(), TMVA::Types::kBDT, m.first, m.second);

        TrainAllMethods(factory);
//        factory->TrainAllMethods();
//        factory->TestAllMethods();
//        factory->EvaluateAllMethods();

        for(const auto& m : methods){
            TString weightfile = ("mydataloader/weights/myFactory_"+m.first+".weights.xml");
            vars->reader->BookMVA(m.first, weightfile);
        }

        auto importance = RankingImportanceVariables(grid, options);
        std::cout<<"importance"<<std::endl;
        auto position = RankingPoisitionVariables(importance);
        std::cout<<"position"<<std::endl;


        std::map<std::string, std::map<SampleId, std::map<size_t, std::vector<double>>>> evaluation;
        for(const auto& m : methods){
            for(const auto& type_entry : vars->data){
                auto sample_entry = type_entry.second;
                for(const auto& entry : sample_entry){
                    evaluation[m.first][entry.first][type_entry.first] = vars->EvaluateForAllEvents(m.first, type_entry.first, entry.first);
                }
            }
        }

        std::map<std::string, std::map<SampleId,double>> kolmogorov = Kolmogorov(evaluation);

        ntuple::MvaTuple mva_tuple(outfile.get(), false);
        auto n_trees_option = options.at("NTrees");
        auto shrinkage_option = options.at("shrinkage");
        auto BaggedSampleFraction_option = options.at("BaggedSampleFraction");
        auto MaxDepth_option = options.at("MaxDepth");
        auto MinNodeSize_option = options.at("MinNodeSize");
        std::cout<<"options"<<std::endl;

        for(const auto& point : grid) {
            auto name = options.GetName(point)+"_"+std::to_string(args.seed());
            size_t pos = point.at(options.GetIndex("NTrees"));
            mva_tuple().NTrees = static_cast<UInt_t>(n_trees_option->GetNumericValue(pos));
            pos = point.at(options.GetIndex("shrinkage"));
            mva_tuple().shrinkage = shrinkage_option->GetNumericValue(pos);

            pos = point.at(options.GetIndex("BaggedSampleFraction"));
            mva_tuple().BaggedSampleFraction = BaggedSampleFraction_option->GetNumericValue(pos);
            pos = point.at(options.GetIndex("MaxDepth"));
            mva_tuple().MaxDepth = MaxDepth_option->GetNumericValue(pos);
            pos = point.at(options.GetIndex("MinNodeSize"));
            mva_tuple().MinNodeSize = MinNodeSize_option->GetNumericValue(pos);
            mva_tuple().ROCIntegral = static_cast<double>(factory->GetROCIntegral(vars->loader.get(), name));

            std::vector<double> ks;
            std::vector<double> mass;
            for (const auto& sample : kolmogorov.at(name)){
                ks.push_back(sample.second);
                mass.push_back(static_cast<double>(sample.first.mass));
            }
            mva_tuple().KS_mass = mass;
            mva_tuple().KS_value = ks;

            std::vector<double> posit;
            std::vector<double> impo;
            std::vector<std::string> variables;
            for (const auto& var: importance.at(name)){
                impo.push_back(var.second);
                variables.push_back(var.first);
                posit.push_back(position.at(name).at(var.first));
            }
            mva_tuple().position = posit;
            mva_tuple().importance = impo;
            mva_tuple().var_name = variables;
            mva_tuple.Fill();
        }
        mva_tuple.Write();

    }
private:
    Arguments args;
    SampleEntryCollection samples;
    MvaSetup mva_setup;
    std::shared_ptr<TFile> outfile;
    std::shared_ptr<TMVA::Factory> factory;
    std::mt19937 gen;
    MvaVariables::VarNameSet enabled_vars;
    std::shared_ptr<MvaVariablesTMVA> vars;
    std::uniform_int_distribution<uint_fast32_t> test_vs_training;
};

}
}

PROGRAM_MAIN(analysis::mva_study::MVATraining, Arguments) // definition of the main program function
//./run.sh MvaTraining --input_path ~/Desktop/tuples --output_file myfile.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name muTau --number_events 10000000

