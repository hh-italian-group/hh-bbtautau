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
    VAR(Float_t, shrinkage) \
    VAR(Float_t, BaggedSampleFraction) \
    VAR(UInt_t, MaxDepth) \
    VAR(Float_t, MinNodeSize) \
    VAR(Float_t, ROCIntegral) \
    VAR(double, KS_Sgn) \
    VAR(double, KS_Bkg) \
    VAR(std::vector<double>, KS_value) \
    VAR(std::vector<double>, KS_mass) \
    VAR(double, KS_260) \
    VAR(double, KS_270) \
    VAR(double, KS_300) \
    VAR(double, KS_320) \
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
    OPT_ARG(int, number_sets, 2);
    OPT_ARG(uint_fast32_t, seed, std::numeric_limits<uint_fast32_t>::max());
};

namespace analysis {
namespace mva_study{



const SampleId mass_tot(SampleType::Sgn_Res, 2000);
const SampleId bkg(SampleType::Bkg_TTbar, 0);
class MvaVariablesTMVA_loader : public MvaVariables{

public:
    using DataVector = std::vector<double>;
    using MassData = std::map<SampleId, std::deque<DataVector>>;
    DataVector variable_loader;
    std::map<size_t, MassData> data;
    std::shared_ptr<TMVA::DataLoader> loader;
    MvaVariablesTMVA_loader(std::shared_ptr<TMVA::DataLoader> dloader,  size_t _number_set = 2,
                            uint_fast32_t _seed = std::numeric_limits<uint_fast32_t>::max(), const VarNameSet& _enabled_vars = {}) :
        MvaVariables(_number_set, _seed, _enabled_vars), loader(dloader) {}

    virtual void SetValue(const std::string& name, double value) override
    {
        if (!names.count(name)){
            variable_loader.push_back(0);
            names[name] = variable_loader.size()-1;
            loader->AddVariable(name);
        }
        variable_loader.at(names.at(name)) = value;
    }

    virtual void AddEventVariables(size_t istraining, const SampleId& mass, double weight) override
    {
        data[istraining][mass].push_back(variable_loader);
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
};

class MvaVariablesTMVA_reader : public MvaVariables{

public:
    using DataVector = std::vector<float>;
    using MassData = std::map<SampleId, std::deque<DataVector>>;
    DataVector variables;
    std::map<std::string, size_t> reader_names;
    std::map<size_t, MassData> data;
    std::shared_ptr<TMVA::Reader> reader;
    MvaVariablesTMVA_reader(std::shared_ptr<TMVA::Reader> _reader, size_t _number_set = 2, uint_fast32_t _seed = std::numeric_limits<uint_fast32_t>::max(),
                            const VarNameSet& _enabled_vars = {}) :
        MvaVariables(_number_set, _seed, _enabled_vars),  reader(_reader) {}

    virtual void SetValue(const std::string& name, double value) override
    {
        if (!reader_names.count(name)){
            variables.push_back(0);
            reader_names[name] = variables.size()-1;
            reader->AddVariable(name, &variables[reader_names.at(name)]);
        }
        variables.at(reader_names.at(name)) = value;
    }

    virtual void AddEventVariables(size_t istraining, const SampleId& mass, double weight) override
    {
        data[istraining][mass].push_back(variables);
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
        dataloader(new TMVA::DataLoader ("mydataloader")), reader(new TMVA::Reader("!V:!!Silent:Color")),
        gen(args.seed()), test_vs_training(100000, std::numeric_limits<uint_fast32_t>::max()), which_set(0, args.number_sets()-1)

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
        uint_fast32_t seed = test_vs_training(gen);
        vars = std::make_shared<MvaVariablesTMVA_loader>(dataloader, args.number_sets(), seed, enabled_vars);
        r_vars = std::make_shared<MvaVariablesTMVA_reader>(reader, args.number_sets(), seed, enabled_vars);
    }

    void RankingVariables(Grid_ND grid, MvaOptionCollection options){
        std::map<std::string, std::map<std::string, double>> var_ranking;
        for(const auto& point : grid) {
            auto name = options.GetName(point)+"_"+std::to_string(args.seed());
            auto method = dynamic_cast<TMVA::MethodBDT*>(factory->GetMethod(dataloader->GetName(), name));
            if(!method)
                throw exception("method '%1%' not found.") % name;
            if (mva_setup.use_mass_var){
                var_ranking["channel"][name] = method->GetVariableImportance(0);               
                var_ranking["mass"][name] = method->GetVariableImportance(1);
            }
            for (size_t i = 0; i<mva_setup.variables.size(); i++){
                var_ranking[mva_setup.variables[i]][name] = method->GetVariableImportance(i+2);
            }
        }
        std::map<std::string, std::shared_ptr<TH1D>> histo_rank;
        auto directory = root_ext::GetDirectory(*outfile.get(), "Ranking");
        for (const auto& r: var_ranking){
            histo_rank[r.first] = std::make_shared<TH1D>((r.first).c_str(),(r.first).c_str(),100, 0, 0.1);
            histo_rank.at(r.first)->SetCanExtend(TH1::kAllAxes);
            for (const auto& entry: r.second){
                histo_rank.at(r.first)->Fill(entry.second);
            }
            root_ext::WriteObject(*histo_rank.at(r.first), directory);
        }
    }

    void Run()
    {
        std::cout<<"Variabili iniziali: "<<enabled_vars.size()<<std::endl;
        std::vector<int> vec_mass;
        auto range = mva_setup.mass_range;
        int begin, end;
        for (int i=0; i<samples.size(); i++){
            vec_mass.push_back(samples[i].id.mass);
            if (samples[i].id.mass == range.min()) begin = i;
            if (samples[i].id.mass == range.max()) end = i;
        }

        for(size_t j = 0; j<mva_setup.channels.size(); j++){
            for(const SampleEntry& entry : samples)
            {
                if (!entry.id.IsBackground() && !range.Contains(entry.id.mass)) continue;
                if ( entry.channel.size() && Parse<Channel>(entry.channel) != mva_setup.channels.at(j) ) continue;
                int mass_background;
                if (entry.id.IsBackground()) {
                    std::uniform_int_distribution<> it(begin, end);
                    mass_background = vec_mass[it(gen)];
                    std::cout<<mass_background<<std::endl;
                }
                auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
                EventTuple tuple(ToString(mva_setup.channels[j]), input_file.get(), true, {} ,GetMvaBranches());
                Long64_t tot_entries = 0, current_entry = 0;

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

                    if (mva_setup.use_mass_var){
                        if (!vars->names.count("channel")){
                            vars->variable_loader.push_back(0);
                            vars->names["channel"] = vars->variable_loader.size()-1;
                            vars->loader->AddVariable("channel");
                        }
                        vars->variable_loader.at(vars->names.at("channel")) = j;

                        if (!vars->names.count("mass")){
                            vars->variable_loader.push_back(0);
                            vars->names["mass"] = vars->variable_loader.size()-1;
                            vars->loader->AddVariable("mass");
                        }
                        if (!entry.id.IsBackground())
                            vars->variable_loader.at(vars->names.at("mass")) = entry.id.mass;
                        else
                            vars->variable_loader.at(vars->names.at("mass")) = mass_background;

                        float channel = static_cast<float>(j);
                        if (!r_vars->reader_names.count("channel")){
                            r_vars->variables.push_back(0);
                            r_vars->reader_names["channel"] = r_vars->variables.size()-1;
                            r_vars->reader->AddVariable("channel", &r_vars->variables[r_vars->reader_names.at("channel")]);
                        }
                        r_vars->variables.at(r_vars->reader_names.at("channel")) = channel;

                        float mass;
                        if (!r_vars->reader_names.count("mass")){
                            r_vars->variables.push_back(0);
                            r_vars->reader_names["mass"] = r_vars->variables.size()-1;
                            r_vars->reader->AddVariable("mass", &r_vars->variables[r_vars->reader_names.at("mass")]);
                        }
                        if (!entry.id.IsBackground()){
                            mass = static_cast<float>(entry.id.mass);
                            r_vars->variables.at(r_vars->reader_names.at("mass")) = mass;
                        }
                        else{
                            mass = static_cast<float>(mass_background);
                            r_vars->variables.at(r_vars->reader_names.at("mass")) = mass;
                        }

                    }
                    vars->AddEvent(event, entry.id, entry.weight);
                    r_vars->AddEvent(event, entry.id, entry.weight);
                }
                std::cout << entry << " number of events: " << tot_entries << std::endl;
            }
        }

        vars->UploadEvents();
        dataloader->PrepareTrainingAndTestTree( "","", "SplitMode=Random" );

        MvaOptionCollection options = mva_setup.CreateOptionCollection(true);
        const Grid_ND grid(options.GetPositionLimits());
        std::map<std::string, std::string> methods;
        for(const auto& point : grid)
            methods[options.GetName(point)+"_"+std::to_string(args.seed())] = options.GetConfigString(point);

        for(const auto& m : methods){
            factory->BookMethod(dataloader.get(), TMVA::Types::kBDT, m.first, m.second);
            TString weightfile = ("mydataloader/weights/myFactory_"+m.first+".weights.xml");
            reader->BookMVA(m.first, weightfile);
        }

        std::cout<<methods.size()<<" METODI!!!"<<std::endl;
        TrainAllMethods(factory);
        RankingVariables(grid, options);

        ntuple::MvaTuple mva_tuple(outfile.get(), false);
        auto n_trees_option = options.at("NTrees");
        auto shrinkage_option = options.at("shrinkage");
//        auto BaggedSampleFraction_option = options.at("BaggedSampleFraction");
//        auto MaxDepth_option = options.at("MaxDepth");
//        auto MinNodeSize_option = options.at("MinNodeSize");

        for(const auto& point : grid) {
            auto name = options.GetName(point)+"_"+std::to_string(args.seed());
            size_t pos = point.at(options.GetIndex("NTrees"));
            mva_tuple().NTrees = n_trees_option->GetNumericValue(pos);
            pos = point.at(options.GetIndex("shrinkage"));
            mva_tuple().shrinkage = shrinkage_option->GetNumericValue(pos);

//            pos = point.at(options.GetIndex("BaggedSampleFraction"));
//            mva_tuple().BaggedSampleFraction = BaggedSampleFraction_option->GetNumericValue(pos);
//            pos = point.at(options.GetIndex("MaxDepth"));
//            mva_tuple().MaxDepth = MaxDepth_option->GetNumericValue(pos);
//            pos = point.at(options.GetIndex("MinNodeSize"));
//            mva_tuple().MinNodeSize = MinNodeSize_option->GetNumericValue(pos);
            mva_tuple().ROCIntegral = factory->GetROCIntegral(dataloader.get(), name);
            mva_tuple.Fill();
        }

        std::map<SampleId, std::map<std::string, std::map<int,std::vector<double>>>> evaluation;
        for(size_t j = 0; j<mva_setup.channels.size(); j++){
            for(const SampleEntry& entry : samples)
            {
                if (!entry.id.IsBackground() && !range.Contains(entry.id.mass)) continue;
                if ( entry.channel.size() && Parse<Channel>(entry.channel) != mva_setup.channels.at(j) ) continue;
                int mass_background;
                if (entry.id.IsBackground()) {
                    std::uniform_int_distribution<> it(begin, end);
                    mass_background = vec_mass[it(gen)];
                    std::cout<<mass_background<<std::endl;
                }
                auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
                EventTuple tuple(ToString(mva_setup.channels[j]), input_file.get(), true, {} ,GetMvaBranches());
                Long64_t tot_entries = 0, current_entry = 0;
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

                    uint_fast32_t seed = test_vs_training(gen);
                    std::mt19937 generator(seed);
                    for(const auto& m : methods){
                        Double_t mvaValue = reader->EvaluateMVA(m.first);
                        evaluation[entry.id][m.first][which_set(generator)].push_back(mvaValue);
                        evaluation[mass_tot][m.first][which_set(generator)].push_back(mvaValue);
                    }
                }
                std::cout << entry << " number of events: " << tot_entries << std::endl;
            }
        }

        std::map<SampleId, std::map<std::string, double>> kolmogorov;

        for(const auto& sample : evaluation){
            auto mass = sample.first;
            for (const auto& bdt : sample.second){
                auto method = bdt.first;
                std::map<int, double*> ks_vector;
                std::map<int, int> dim;
                for (auto tvt : bdt.second){
                    auto type = tvt.first;
                    std::sort(tvt.second.begin(), tvt.second.end());
                    dim[type] = tvt.second.size();
                    ks_vector[type] = tvt.second.data();
                }
                double ks = TMath::KolmogorovTest(dim.at(0), ks_vector.at(0), dim.at(1), ks_vector.at(1), "");
                std::cout<<ks<<std::endl;
                kolmogorov[mass][method] = ks;
            }
        }

        for(const auto& point : grid) {
            auto name = options.GetName(point)+"_"+std::to_string(args.seed());
            mva_tuple().KS_Bkg = kolmogorov.at(bkg).at(name);
            mva_tuple().KS_Sgn = kolmogorov.at(mass_tot).at(name);
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
    std::shared_ptr<TMVA::DataLoader> dataloader;
    std::shared_ptr<TMVA::Reader> reader;
    std::mt19937 gen;
    MvaVariables::VarNameSet enabled_vars;
    std::shared_ptr<MvaVariablesTMVA_loader> vars;
    std::shared_ptr<MvaVariablesTMVA_reader> r_vars;
    std::uniform_int_distribution<uint_fast32_t> test_vs_training;
    std::uniform_int_distribution<size_t> which_set;
};

}
}

PROGRAM_MAIN(analysis::mva_study::MVATraining, Arguments) // definition of the main program function
//./run.sh MvaTraining --input_path ~/Desktop/tuples --output_file myfile.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name muTau --number_events 10000000

