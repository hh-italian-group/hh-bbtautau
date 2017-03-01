/*! Study for Mva Training
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
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
    REQ_ARG(std::string, tree_name);
    REQ_ARG(Long64_t, number_events);
};

namespace analysis {

struct SampleEntry{
  std::string filename;

  double weight;
  bool issignal;
  SampleEntry() : weight(-1), issignal(false){}
};

std::ostream& operator<<(std::ostream& os, const SampleEntry& entry)
{
    os << entry.filename << " " << entry.weight << " " << std::boolalpha << entry.issignal;
    return os;
}

std::istream& operator>>(std::istream& is, SampleEntry& entry)
{
    is >> entry.filename >> entry.weight >> std::boolalpha >> entry.issignal;
    return is;
}

class MvaVariables{
private:
    std::vector<double> variables;
    std::map<std::string, size_t> names;
    std::shared_ptr<TMVA::DataLoader> loader;
public:
    MvaVariables(std::shared_ptr<TMVA::DataLoader> dloader) : loader(dloader){}

    double& operator[](const std::string& name)
    {
      if (!names.count(name)){
          variables.push_back(0);
          names[name]=variables.size()-1;
          loader->AddVariable(name);
      }
      return variables.at(names.at(name));
    }

    void AddEvent(bool issignal, bool istraining, double weight)
    {
        const std::string samplename = issignal ? "Signal" : "Background";
        const TMVA::Types::ETreeType treetype= istraining ? TMVA::Types::kTraining : TMVA::Types::kTesting;
        loader->AddEvent(samplename, treetype, variables, weight);
    }


};

using SampleEntryCollection = std::vector<SampleEntry>;

class MvaClassification {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    static const std::set<std::string>& GetDisabledBranches()
    {
        static const std::set<std::string> DisabledBranches_read = {
            "dphi_mumet", "dphi_metsv", "dR_taumu", "mT1", "mT2", "dphi_bbmet", "dphi_bbsv", "dR_bb", "m_bb", "n_jets",
            "btag_weight", "ttbar_weight",  "PU_weight", "shape_denominator_weight", "trigger_accepts", "trigger_matches",
            "event.tauId_keys_1","event.tauId_keys_2","event.tauId_values_1","event.tauId_values_2"
        };
        return DisabledBranches_read;
    }


    static SampleEntryCollection ReadConfig(const std::string& cfg_file){
        std::ifstream f(cfg_file);
        SampleEntryCollection collection;
        while(f.good()){
            std::string line;
            std::getline(f, line);
            if (line.size()==0 || line.at(0)=='#') continue;
            std::istringstream s(line);
            SampleEntry entry;
            s>>entry;
            collection.push_back(entry);
        }
        return collection;
    }

    MvaClassification(const Arguments& _args): args(_args), samples(ReadConfig(args.cfg_file())),
        outfile(root_ext::CreateRootFile(args.output_file())),
        factory(new TMVA::Factory ("myFactory", outfile.get(),"!V:!Silent:Color:DrawProgressBar:Transformations=I,D:AnalysisType=Classification")),
        dataloader(new TMVA::DataLoader ("mydataloader")), vars(dataloader),
        gen(rd()), testvstraining(0,1)

    {
    }

    void Run()
    {
        double c=0;
        for(const SampleEntry& entry:samples)
        {
            auto input_file=root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
            EventTuple tuple(args.tree_name(), input_file.get(), true, GetDisabledBranches());
            std::cout<<entry<<" number of events: "<<std::min(tuple.GetEntries(),args.number_events())<<std::endl;

            for(Long64_t current_entry = 0; current_entry < std::min(tuple.GetEntries(),args.number_events()); ++current_entry) {
                tuple.GetEntry(current_entry);
                const Event& event = tuple.data();
                if (event.eventEnergyScale!=0 || (event.q_1+event.q_2)!=0 || event.jets_p4.size() < 2
                    || event.extraelec_veto==true || event.extramuon_veto==true) continue;

                const bool istraining= testvstraining(gen) == 0;

                LorentzVectorE_Float bb= event.jets_p4[0] + event.jets_p4[1];
//                LorentzVectorM_Float leptons= event.p4_1 + event.p4_2;

                double circular_cut=std::sqrt(pow(event.SVfit_p4.mass()-116.,2)+pow(bb.M()-111,2));
                if (circular_cut>40) continue;

                vars["dphi_l1met"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.pfMET_p4));
                vars["dphi_svmet"]=std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4));
                vars["dR_jetsPt_bb"]=std::abs(ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]))*bb.Pt();
                vars["dR_leptonsPt_sv"] = std::abs(ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2))*event.SVfit_p4.Pt();
                vars["MT_l1"]=Calculate_MT(event.p4_1,event.pfMET_p4);
                vars["MT_l2"]=Calculate_MT(event.p4_2,event.pfMET_p4);
                vars["dphi_bbmet"]=std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.pfMET_p4));
                vars["dphi_bbsv"]=std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.SVfit_p4));

//                vars["mass_sv"]=event.SVfit_p4.M();
//                vars["mass_jets"]=bb.M();
//                vars["mt_l1"]=event.p4_1.mt(); // Calculate_MT(event.p4_1,event.pfMET_p4)??
//                vars["mt_l2"]=event.p4_2.mt();
//                vars["dphi_l2met"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_2, event.pfMET_p4));
//                vars["pt_met"]=event.pfMET_p4.pt();
//                vars["dphi_leptons"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.p4_2));
//                vars["dphi_jets"]=std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.jets_p4[0], event.jets_p4[1]));
//                vars["deta_leptons"] =std::abs(event.p4_1.eta() - event.p4_2.eta());
//                vars["deta_jets"]=std::abs(event.jets_p4[0].eta() - event.jets_p4[1].eta());
//                vars["dR_leptons"] = std::abs(ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2));
//                vars["dR_jets"]=std::abs(ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]));
//                vars["transversemass_hl"]=std::sqrt(2*leptons.pt()*event.pfMET_p4.Et()*(1-std::cos(std::abs(ROOT::Math::VectorUtil::DeltaPhi(leptons, event.pfMET_p4)))));
//                vars["transversemass_hl"]=Calculate_MT(leptons, event.pfMET_p4);
//                vars["mass_H"]=ROOT::Math::VectorUtil::InvariantMass(bb,event.SVfit_p4);

                vars.AddEvent(entry.issignal, istraining, entry.weight);

            }
        }

        dataloader->PrepareTrainingAndTestTree( "","", "SplitMode=Random" );


        factory->BookMethod( dataloader.get(), TMVA::Types::kBDT, "BDTG",
                                    "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );


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
    MvaVariables vars;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<> testvstraining;
};


}


PROGRAM_MAIN(analysis::MvaClassification, Arguments) // definition of the main program function
/*./run.sh MvaTraining --input_path ~/Desktop/tuples --output_file myfile.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name muTau --number_events 10000000

