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
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include <fstream>
#include <random>


struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(std::string, variables_file);
    REQ_ARG(Long64_t, range);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(Long64_t, number_events);
    OPT_ARG(int, number_sets, 2);
};

namespace analysis {
namespace mva_study{

namespace {
Range<int> low_mass(250, 280),  medium_mass1(300, 350), medium_mass2(400, 550), high_mass(600,900);
std::vector<Range<int>> ranges{low_mass, medium_mass1, medium_mass2, high_mass};
}

static std::unordered_set<std::string> ReadVars(const std::string& var_file){
    std::ifstream f(var_file);
    std::unordered_set<std::string> enabled_vars;
    while(f.good()){
        std::string var;
        std::getline ( f, var );
        if (var.size()==0)
            continue;
        enabled_vars.insert(var);
    }
    return enabled_vars;
}

class MvaVariablesTMVA : public MvaVariables{
private:
    using DataVector = std::vector<double>;
    using MassData = std::map<SampleId, std::deque<DataVector>>;

    DataVector variables;
    std::map<std::string, size_t> names;
    std::map<size_t, MassData> data;
    std::shared_ptr<TMVA::DataLoader> loader;
    VarNameSet enabled;

public:
    MvaVariablesTMVA(std::shared_ptr<TMVA::DataLoader> dloader, uint_fast32_t seed, const VarNameSet& _enabled_vars) :
        MvaVariables(), loader(dloader), enabled(_enabled_vars) {}

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
            if(m.mass == -1) weights[m] = 1;
            else weights[m] = 1. / (data[0].at(m).size() + data[1][m].size() );
        }

        for(const auto& data_entry : data) {
            const bool istraining = data_entry.first;
            const TMVA::Types::ETreeType treetype = istraining ? TMVA::Types::kTraining : TMVA::Types::kTesting;
            for(const auto& m_entry : data_entry.second) {
                const SampleId m = m_entry.first;
                const std::string samplename = m.mass != -1 ? "Signal" : "Background";
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

//    static const std::set<std::string>& GetDisabledBranches()
//    {
//        static const std::set<std::string> DisabledBranches_read = {
//            "dphi_mumet", "dphi_metsv", "dR_taumu", "mT1", "mT2", "dphi_bbmet", "dphi_bbsv", "dR_bb", "m_bb", "n_jets",
//            "btag_weight", "ttbar_weight",  "PU_weight", "shape_denominator_weight", "trigger_accepts", "trigger_matches",
//            "event.tauId_keys_1","event.tauId_keys_2","event.tauId_values_1","event.tauId_values_2"
//        };
//        return DisabledBranches_read;
//    }


    MVATraining(const Arguments& _args): args(_args), samples(SampleEntry::ReadConfig(args.cfg_file())),
        outfile(root_ext::CreateRootFile(args.output_file())),
        factory(new TMVA::Factory ("myFactory", outfile.get(),"!V:!Silent:Color:DrawProgressBar:Transformations=I,D:AnalysisType=Classification")),
        dataloader(new TMVA::DataLoader ("mydataloader")), enabled_vars(ReadVars(args.variables_file())),
        vars(dataloader, UINT_FAST32_MAX, enabled_vars)

    {
    }

    void Run()
    {

        std::cout<<"Variabili iniziali: "<<enabled_vars.size()<<std::endl;
        for(const auto& range : ranges){
            if (!range.Contains(args.range())) continue;
            for(const SampleEntry& entry : samples)
            {
                if (entry.id.mass != -1 && !range.Contains(entry.id.mass)) continue;
                if ( entry.channel != "" && args.tree_name() != entry.channel ) continue;
                auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
                EventTuple tuple(args.tree_name(), input_file.get(), true, {} ,GetMvaBranches());
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
                    vars.AddEvent(event, entry.id, entry.weight);
                }
                std::cout << entry << " number of events: " << tot_entries << std::endl;
            }
        }

        vars.UploadEvents();
        dataloader->PrepareTrainingAndTestTree( "","", "SplitMode=Random" );

        factory->BookMethod( dataloader.get(), TMVA::Types::kBDT, "BDTG",
                             "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=6");


        factory->BookMethod( dataloader.get(), TMVA::Types::kBDT, "BDTD_2","!H:!V:NTrees=200:MaxDepth=3:MinNodeSize=2.5%:"
                                        "nCuts=-1:BoostType=Bagging:UseBaggedBoost=True:BaggedSampleFraction=0.6:UseRandomisedTrees=True:"
                                        "UseNvars=2:SeparationType=GiniIndex:UseYesNoLeaf:VarTransform=Decorrelate" );


        factory->TrainAllMethods();
        factory->TestAllMethods();
        factory->EvaluateAllMethods();

    }
private:
    Arguments args;
    SampleEntryCollection samples;
    std::shared_ptr<TFile> outfile;
    std::shared_ptr<TMVA::Factory> factory;
    std::shared_ptr<TMVA::DataLoader> dataloader;
    MvaVariables::VarNameSet enabled_vars;
    MvaVariablesTMVA vars;
};

}
}


PROGRAM_MAIN(analysis::mva_study::MVATraining, Arguments) // definition of the main program function
//./run.sh MvaTraining --input_path ~/Desktop/tuples --output_file myfile.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name muTau --number_events 10000000

