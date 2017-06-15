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

namespace {
Range<int> sm(0,200), low_mass(250, 320), medium_mass(340,400), high_mass(450,900);
std::vector<Range<int>> ranges{sm, low_mass, medium_mass, high_mass};
}

class MvaVariablesTMVA : public MvaVariables{

public:
    using DataVector = std::vector<double>;
    using MassData = std::map<SampleId, std::deque<DataVector>>;
    DataVector variables;
    std::map<size_t, MassData> data;
    std::shared_ptr<TMVA::DataLoader> loader;
    MvaVariablesTMVA(std::shared_ptr<TMVA::DataLoader> dloader, size_t _number_set = 2, uint_fast32_t _seed = std::numeric_limits<uint_fast32_t>::max(), const VarNameSet& _enabled_vars = {}) :
        MvaVariables(_number_set, _seed, _enabled_vars), loader(dloader){}

    virtual void SetValue(const std::string& name, double value) override
    {
        if (!names.count(name)){
            variables.push_back(0);
            names[name] = variables.size()-1;
            loader->AddVariable(name);
        }
        variables.at(names.at(name)) = value;
    }

    virtual void AddEventVariables(size_t istraining, const SampleId& mass, double weight) override
    {
        data[istraining][mass].push_back(variables);
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
            std::cout<<istraining<<std::endl;
            const TMVA::Types::ETreeType treetype = istraining ? TMVA::Types::kTraining : TMVA::Types::kTesting;
            for(const auto& m_entry : data_entry.second) {
                const SampleId m = m_entry.first;
                const std::string samplename = m.IsSignal() ? "Signal" : "Background";
                std::cout<<ToString(m)<<samplename<<std::endl;
                for(const auto& vars : m_entry.second) {
                    loader->AddEvent(samplename, treetype, vars, weights.at(m));
                }
            }
        }
    }

};

class MVATraining {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    MVATraining(const Arguments& _args): args(_args),
        outfile(root_ext::CreateRootFile(args.output_file()+"_"+std::to_string(args.seed())+".root")),
        factory(new TMVA::Factory ("myFactory", outfile.get(),"!V:!Silent:Color:DrawProgressBar:Transformations=I,D:AnalysisType=Classification")),
        dataloader(new TMVA::DataLoader ("mydataloader")), gen(args.seed()), test_vs_training(100000, std::numeric_limits<uint_fast32_t>::max())

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

//        std::unordered_set<std::string> enabled();
        enabled_vars.insert(mva_setup.variables.begin(), mva_setup.variables.end());
        uint_fast32_t seed = test_vs_training(gen);
        std::cout<<"SEED "<<seed<<std::endl;
        vars = std::make_shared<MvaVariablesTMVA>(dataloader, args.number_sets(), seed, enabled_vars);
    }

    void Run()
    {

        std::cout<<"Variabili iniziali: "<<enabled_vars.size()<<std::endl;

        std::cout<<"arg = "<<args.seed()<<std::endl;
        std::cout<<"std::numeric_limits<uint_fast32_t>::max() = "<<std::numeric_limits<uint_fast32_t>::max()<<std::endl;
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
                            vars->variables.push_back(0);
                            vars->names["channel"] = vars->variables.size()-1;
                            vars->loader->AddVariable("channel");
                        }
                        vars->variables.at(vars->names.at("channel")) = j;


                        if (!vars->names.count("mass")){
                            vars->variables.push_back(0);
                            vars->names["mass"] = vars->variables.size()-1;
                            vars->loader->AddVariable("mass");
                        }
                        if (!entry.id.IsBackground()) {
                            vars->variables.at(vars->names.at("mass")) = entry.id.mass;
                        }
                        else {
                            vars->variables.at(vars->names.at("mass")) = mass_background;
                        }
                    }

                    vars->AddEvent(event, entry.id, entry.weight);
                }
                std::cout << entry << " number of events: " << tot_entries << std::endl;
            }
        }


        std::cout<<vars->variables.size()<<std::endl;
        std::cout << "Upload" << std::endl;
        vars->UploadEvents();
        std::cout << "Test e training" << std::endl;

        dataloader->PrepareTrainingAndTestTree( "","", "SplitMode=Random" );

        MvaOptionCollection options = mva_setup.CreateOptionCollection(true);
        const Grid_ND grid(options.GetPositionLimits());

        std::map<std::string, std::string> methods;

        for(const auto& point : grid) {
            methods[options.GetName(point)+"_"+std::to_string(args.seed())] = options.GetConfigString(point);
        }

        std::map<std::string, TMVA::MethodBDT*> methodbase;
        for(const auto& m : methods){
            methodbase[m.first] = dynamic_cast<TMVA::MethodBDT*>(factory->BookMethod(dataloader.get(), TMVA::Types::kBDT, m.first, m.second));
            if (!methodbase.at(m.first))
                throw exception("No method found.");
        }

        std::cout<<methods.size()<<" METODI!!!"<<std::endl;

        factory->TrainAllMethods();
        factory->TestAllMethods();
        factory->EvaluateAllMethods();

        std::map<std::string, std::map<std::string, int>> var_ranking;

        ntuple::MvaTuple mva_tuple(outfile.get(), false);
        auto n_trees_option = options.at("NTrees");
        auto shrinkage_option = options.at("shrinkage");
//        auto BaggedSampleFraction_option = options.at("BaggedSampleFraction");
//        auto MaxDepth_option = options.at("MaxDepth");
//        auto MinNodeSize_option = options.at("MinNodeSize");

        for(const auto& point : grid) {
            auto name = options.GetName(point)+"_"+std::to_string(args.seed());
            auto config = options.GetConfigString(point);
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
//            mva_tuple().ROCIntegral = factory->GetROCIntegral(dataloader.get(), name);
            std::cout<<methodbase.at(name)->GetKSTrainingVsTest('S')<<std::endl;
//            mva_tuple().KS_Sgn = methodbase.at(name)->GetKSTrainingVsTest('S');
//            mva_tuple().KS_Bkg = methodbase.at(name)->GetKSTrainingVsTest('B');
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
    std::mt19937 gen;
    MvaVariables::VarNameSet enabled_vars;
    std::shared_ptr<MvaVariablesTMVA> vars;
    std::uniform_int_distribution<uint_fast32_t> test_vs_training;
};

}
}

PROGRAM_MAIN(analysis::mva_study::MVATraining, Arguments) // definition of the main program function
//./run.sh MvaTraining --input_path ~/Desktop/tuples --output_file myfile.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name muTau --number_events 10000000

