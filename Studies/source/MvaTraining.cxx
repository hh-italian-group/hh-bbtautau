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
         const double mass_top = 173.21;
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
                LorentzVectorM_Float leptons= event.p4_1 + event.p4_2;
                LorentzVectorM_Float leptonsMET= event.p4_1 + event.p4_2 + event.pfMET_p4;

                double circular_cut=std::sqrt(pow(event.SVfit_p4.mass()-116.,2)+pow(bb.M()-111,2));
                if (circular_cut>40) continue;
//                if ((args.tree_name=="eTau" && circular_cut>40)||(args.tree_name=="muTau" && circular_cut>30)) continue;

//                vars["pt_l1"] = event.p4_1.pt();
//                vars["pt_l2"] = event.p4_2.pt();
//                vars["pt_b1"] = event.jets_p4[0].pt();
//                vars["pt_b2"] = event.jets_p4[1].pt();
//                vars["pt_l1l2"] = leptons.pt();
//                vars["pt_htautau"] = event.SVfit_p4.pt();
//                vars["pt_l1l2met"] = leptonsMET.pt();
//                vars["pt_hbb"] = bb.pt();
//                vars["pt_met"] = event.pfMET_p4.pt();

//                vars["p_zeta"] = Calculate_Pzeta(event.p4_1, event.p4_2, event.pfMET_p4);
//                vars["p_zeta_visible"] = Calculate_visiblePzeta(event.p4_1,event.p4_2);

//                vars["abs(dphi_l1l2)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.p4_2));
                vars["dphi_l1l2"] = (ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.p4_2));
//                vars["abs(dphi_b1b2)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.jets_p4[0], event.jets_p4[1]));
//                vars["dphi_b1b2"] = (ROOT::Math::VectorUtil::DeltaPhi(event.jets_p4[0], event.jets_p4[1]));
//                vars["abs(dphi_l1met)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.pfMET_p4));
//                vars["dphi_l1met"] = (ROOT::Math::VectorUtil::DeltaPhi(event.p4_1, event.pfMET_p4));
//                vars["abs(dphi_l2met)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.p4_2, event.pfMET_p4));
//                vars["dphi_l2met"] = (ROOT::Math::VectorUtil::DeltaPhi(event.p4_2, event.pfMET_p4));
//                vars["abs(dphi_l1l2met)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(leptons, event.pfMET_p4));
//                vars["dphi_l1l2met"] = (ROOT::Math::VectorUtil::DeltaPhi(leptons, event.pfMET_p4));
//                vars["abs_dphi_htautaumet"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4));
                vars["dphi_htautaumet"] = (ROOT::Math::VectorUtil::DeltaPhi(event.SVfit_p4, event.pfMET_p4));
//                vars["abs(dphi_hbbmet)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.pfMET_p4));
                vars["dphi_hbbmet"] = (ROOT::Math::VectorUtil::DeltaPhi(bb, event.pfMET_p4));
//                vars["abs(dphi_hbbhatutau)"] = std::abs(ROOT::Math::VectorUtil::DeltaPhi(bb, event.SVfit_p4));
                vars["dphi_hbbhtautau"] = (ROOT::Math::VectorUtil::DeltaPhi(bb, event.SVfit_p4));

//                vars["abs_deta_l1l2"] = std::abs(event.p4_1.eta() - event.p4_2.eta());
                vars["deta_l1l2"] = (event.p4_1.eta() - event.p4_2.eta());
//                vars["abs_deta_b1b2"] = std::abs(event.jets_p4[0].eta() - event.jets_p4[1].eta());
                vars["deta_b1b2"] = (event.jets_p4[0].eta() - event.jets_p4[1].eta());
//                vars["abs_deta_l1met"] = std::abs(event.p4_1.eta()-event.pfMET_p4.eta());
//                vars["deta_l1met"] = (event.p4_1.eta()-event.pfMET_p4.eta());
//                vars["abs(deta_l2met)"] = std::abs(event.p4_2.eta()-event.pfMET_p4.eta());
//                vars["deta_l2met"] = (event.p4_2.eta()-event.pfMET_p4.eta());
//                vars["abs(deta_l1l2met)"] = std::abs(leptons.eta()-event.pfMET_p4.eta());
//                vars["deta_l1l2met"] = (leptons.eta()-event.pfMET_p4.eta());
//                vars["abs_deta_htautaumet"] = std::abs(event.SVfit_p4.eta()-event.pfMET_p4.eta());
//                vars["deta_hatutaumet"] = (event.SVfit_p4.eta()-event.pfMET_p4.eta());
//                vars["abs(deta_hbbmet)"] = std::abs(bb.eta()-event.pfMET_p4.eta());
//                vars["deta_hbbmet"] = (bb.eta()-event.pfMET_p4.eta());
//                vars["abs_deta_hbbhtautau"] = std::abs(bb.eta()-event.SVfit_p4.eta());
                vars["deta_hbbhatutau"] = bb.eta()-event.SVfit_p4.eta();

//                vars["dR_l1l2"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2));
//                vars["dR_b1b2"] = (ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]));
//                vars["dR_l1met"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.pfMET_p4));
//                vars["dR_l2met"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_2, event.pfMET_p4));
//                vars["dR_l1l2met"] = (ROOT::Math::VectorUtil::DeltaR(leptons, event.pfMET_p4));
                vars["dR_htautaumet"] = (ROOT::Math::VectorUtil::DeltaR(event.SVfit_p4, event.pfMET_p4));
                vars["dR_hbbmet"] = (ROOT::Math::VectorUtil::DeltaR(bb, event.pfMET_p4));
//                vars["dR_hbbhtautau"] = (ROOT::Math::VectorUtil::DeltaR(bb, event.SVfit_p4));

//                vars["dR_b1b2Pt_hbb"] = (ROOT::Math::VectorUtil::DeltaR(event.jets_p4[0], event.jets_p4[1]))*bb.Pt();
//                vars["dR_l1l2Pt_htautau"] = (ROOT::Math::VectorUtil::DeltaR(event.p4_1, event.p4_2))*event.SVfit_p4.Pt();

//                vars["mass_l1l2met"] = ROOT::Math::VectorUtil::InvariantMass(leptons,event.pfMET_p4); //
                vars["mass_htautau"] = event.SVfit_p4.M();
//                vars["mass_l1l2"] = std::sqrt(pow(event.p4_1.Et()+event.p4_2.Et(),2)-pow(event.p4_1.px()+event.p4_2.px(),2)+pow(event.p4_1.py()+event.p4_2.py(),2));//
                vars["mass_hbb"] = bb.M();
//                vars["MT_l1"] = Calculate_MT(event.p4_1,event.pfMET_p4);
//                vars["MT_l2"] = Calculate_MT(event.p4_2,event.pfMET_p4);
                vars["MT_hatutau"] = Calculate_MT(event.SVfit_p4, event.pfMET_p4);
                vars["MT_l1l2"] = Calculate_MT(leptons, event.pfMET_p4);
                vars["MT_tot"] = Calculate_TotalMT(event.p4_1,event.p4_2,event.pfMET_p4); //Total transverse mass
//                vars["MT2"] = Calculate_MT2(event.p4_1,event.p4_2,event.jets_p4[0], event.jets_p4[1], event.pfMET_p4); //Stransverse mass
//                vars["mass_H"] = ROOT::Math::VectorUtil::InvariantMass(bb,event.SVfit_p4);

//                LorentzVectorM_Float a1 = event.p4_1 + event.jets_p4[0] + event.pfMET_p4;
//                LorentzVectorM_Float b1 = event.p4_2 + event.jets_p4[1];
//                LorentzVectorM_Float a2 = event.p4_1 + event.jets_p4[0];
//                LorentzVectorM_Float b2 = event.p4_2 + event.jets_p4[1] + event.pfMET_p4;
//                LorentzVectorM_Float a3 = event.p4_1 + event.jets_p4[1] + event.pfMET_p4;
//                LorentzVectorM_Float b3 = event.p4_2 + event.jets_p4[0];
//                LorentzVectorM_Float a4 = event.p4_1 + event.jets_p4[1];
//                LorentzVectorM_Float b4 = event.p4_2 + event.jets_p4[0] + event.pfMET_p4;

//                double d1 = pow(std::abs(a1.mass() - mass_top),2) + pow (std::abs(b1.mass() - mass_top),2);
//                double d2 = pow(std::abs(a2.mass() - mass_top),2) + pow (std::abs(b2.mass() - mass_top),2);
//                double d3 = pow(std::abs(a3.mass() - mass_top),2) + pow (std::abs(b3.mass() - mass_top),2);
//                double d4 = pow(std::abs(a4.mass() - mass_top),2) + pow (std::abs(b4.mass() - mass_top),2);

//                if (d1<d2 && d1<d3 && d1<d4) {
//                    vars["Mass_top1"] = a1.mass();
//                    vars["Mass_top2"] = b1.mass();
//                }
//                if (d2<d1 && d2<d3 && d2<d4) {
//                    vars["Mass_top1"] = a2.mass();
//                    vars["Mass_top2"] = b2.mass();
//                }
//                if (d3<d1 && d3<d2 && d3<d4) {
//                    vars["Mass_top1"] = a3.mass();
//                    vars["Mass_top2"] = b3.mass();
//                }
//                if (d4<d1 && d4<d3 && d4<d2) {
//                    vars["Mass_top1"] = a4.mass();
//                    vars["Mass_top2"] = b4.mass();
//                }

                const analysis::LorentzVectorXYZ sv(event.SVfit_p4.px(),event.SVfit_p4.py(),event.SVfit_p4.pz(),event.SVfit_p4.e());
                const auto boosted_tau1 = ROOT::Math::VectorUtil::boost(event.p4_1, sv.BoostToCM());
//                vars["theta_l1(h)"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_tau1, sv)); //theta angle between the first final state lepton and the direction of flight of h_tautau in the h_tautau rest frame
//                vars["phi_l1(h)"] = std::atan(boosted_tau1.py()/boosted_tau1.px()); //phi angle between the first final state lepton and the direction of flight of h_tautau in the h_tautau rest frame
                const auto boosted_tau2 = ROOT::Math::VectorUtil::boost(event.p4_2, sv.BoostToCM());
//                vars["theta_l2(h)"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_tau2, sv)); //angle between the second final state lepton and the direction of flight of h_tautau in the h_tautau rest frame
//                vars["phi_l2(h)"] = std::atan(boosted_tau2.py()/boosted_tau2.px()); //phi angle between the second final state lepton and the direction of flight of h_tautau in the h_tautau rest frame
//                vars["R_l1l2"] = ROOT::Math::VectorUtil::DeltaR(boosted_tau1, boosted_tau2); // R between the two final state leptons in the h_tautau rest frame

                const analysis::LorentzVectorXYZ hbb(bb.px(),bb.py(),bb.pz(),bb.e());
                const auto boosted_b1 = ROOT::Math::VectorUtil::boost(event.jets_p4[0], hbb.BoostToCM());
//                vars["theta_b1(h)"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_b1, hbb)); //angle between the first final state bjet and the direction of flight of h_bb in the h_bb rest frame
//                vars["phi_b1(h)"] = std::atan(boosted_b1.py()/boosted_b1.px()); //phi angle between the first final state bjet and the direction of flight of h_bb in the h_bb rest frame
                const auto boosted_b2 = ROOT::Math::VectorUtil::boost(event.jets_p4[2], hbb.BoostToCM());
//                vars["theta_b2(h)"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_b2, hbb)); //angle between the second final state bjet and the direction of flight of h_bb in the h_bb rest frame
//                if (boosted_b2.px()!=0) vars["phi_b2(h)"] = std::atan(boosted_b2.py()/boosted_b2.px()); //phi angle between the second final state bjet and the direction of flight of h_bb in the h_bb rest frame
//                vars["R_b1b2"] = ROOT::Math::VectorUtil::DeltaR(boosted_b1, boosted_b2); // R between the two final state b-jetsin the h_bb rest frame

                LorentzVectorE_Float H = bb + event.SVfit_p4;
                const analysis::LorentzVectorXYZ vec_H(H.px(),H.py(),H.pz(),H.e());
                const auto boosted_l1 = ROOT::Math::VectorUtil::boost(event.p4_1, vec_H.BoostToCM());
                const auto boosted_l2 = ROOT::Math::VectorUtil::boost(event.p4_2, vec_H.BoostToCM());
                const auto boosted_j1 = ROOT::Math::VectorUtil::boost(event.jets_p4[0], vec_H.BoostToCM());
                const auto boosted_j2 = ROOT::Math::VectorUtil::boost(event.jets_p4[1], vec_H.BoostToCM());
                const TVector3 vec_l1(boosted_l1.px(),boosted_l1.py(),boosted_l1.pz());
                const TVector3 vec_l2(boosted_l2.px(),boosted_l2.py(),boosted_l2.pz());
                const TVector3 vec_j1(boosted_j1.px(),boosted_j1.py(),boosted_j1.pz());
                const TVector3 vec_j2(boosted_j2.px(),boosted_j2.py(),boosted_j2.pz());
                const auto n1 = vec_l1.Cross(vec_l2);
                const auto n2 = vec_j1.Cross(vec_j2);
                vars["phi"] = ROOT::Math::VectorUtil::Angle(n1, n2); //angle between the decay planes of the four final state elements expressed in the H rest frame

                const auto boosted_htautau = ROOT::Math::VectorUtil::boost(event.SVfit_p4, vec_H.BoostToCM());
                vars["theta_star1"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_htautau, ROOT::Math::Cartesian3D<>(0, 0, 1))); // Is the production angle of the h_tautau defined in the H rest frame

                const auto boosted_hbb = ROOT::Math::VectorUtil::boost(bb, vec_H.BoostToCM());
                vars["theta_star2"] = std::acos(ROOT::Math::VectorUtil::CosTheta(boosted_hbb, ROOT::Math::Cartesian3D<>(0, 0, 1)));// Is the production angle of the h_bb defined in the H rest frame

                const TVector3 vec_htautau(boosted_htautau.px(),boosted_htautau.py(),boosted_htautau.pz());
                TVector3 z_axis(0,0,1);
                const auto n3 = vec_htautau.Cross(z_axis);
                vars["phi_1"] = ROOT::Math::VectorUtil::Angle(n1,n3); //Angle between the decay plane of the lepton pair and a plane defined by the vector of the h_tautau in the H rest frame and the positive direction of z axis

                const TVector3 vec_hbb(boosted_hbb.px(),boosted_hbb.py(),boosted_hbb.pz());
                const auto n4 = vec_hbb.Cross(z_axis);
                vars["phi_2"] = ROOT::Math::VectorUtil::Angle(n2,n4); //Angle between the decay plane of the b-jets pair and a plane defined by the vector of the h_bb in the H rest frame and the positive direction of z axis


                vars.AddEvent(entry.issignal, istraining, entry.weight);

            }
        }

        dataloader->PrepareTrainingAndTestTree( "","", "SplitMode=Random" );


        factory->BookMethod( dataloader.get(), TMVA::Types::kBDT, "BDTG",
                             "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=6");


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

//./run.sh MvaTraining --input_path ~/Desktop/tuples --output_file BDT_lm_muTau.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name muTau --number_events 10000000

