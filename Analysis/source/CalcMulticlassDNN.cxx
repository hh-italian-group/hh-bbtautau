/*! Writes multiclass outputs to a new tree that can be used as a friend. */

#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventInfo.h"

#include "MulticlassInference/MulticlassInference/interface/hmc.h"

#define CHECK_EMPTY(COND, EXPR) ((COND) ? (EXPR) : (hmc::features::EMPTY))

namespace analysis {

struct CalcMulticlassDNNArguments {
  REQ_ARG(Channel, channel);
  REQ_ARG(Period, period);
  REQ_ARG(std::string, input);
  REQ_ARG(std::string, output);
  OPT_ARG(int, end, -1);
  OPT_ARG(int, progress, 100);
};

class FeatureProvider {
  public:
  FeatureProvider(Period period, Channel channel, TTree* inTree);

  FeatureProvider(const FeatureProvider&) = delete;

  void calculate();

  inline hmc::EventId getEventId() const {
    return hmc::EventId(ulong64Inputs_.at("evt"));
  }

  inline void add(const std::string& featureName) {
    features_.emplace(featureName, 0.);
  }

  inline bool has(const std::string& featureName) {
    return features_.count(featureName) == 1;
  }

  inline float get(const std::string& featureName) const {
    const auto& it = features_.find(featureName);
    if (it == features_.end()) {
      throw exception("FeatureProvider: unknown feature '" + featureName + "'");
    }
    return it->second;
  }

  private:
  Period period_;
  Channel channel_;
  std::map<std::string, bool> boolInputs_;
  std::map<std::string, int> intInputs_;
  std::map<std::string, ULong64_t> ulong64Inputs_;
  std::map<std::string, float> floatInputs_;
  std::map<std::string, float> features_;

  // currently not needed
  // bool passBaseline_() const;
  // bool passEllipseMassCut_(const TLorentzVector& bH, const TLorentzVector& tauH) const;
  // bool passRectMassCut_(const TLorentzVector& bH, const TLorentzVector& tauH) const;
};

class CalcMulticlassDNN {
  public:
  CalcMulticlassDNN(const CalcMulticlassDNNArguments& args)
      : args_(args)
      , inputFile_(root_ext::OpenRootFile(args_.input()))
      , outputFile_(root_ext::CreateRootFile(args_.output())) {
  }

  void Run() {
    std::cout << "running multiclass classification inference" << std::endl;
    std::cout << "channel : " << args_.channel() << std::endl;
    std::cout << "period  : " << args_.period() << std::endl;
    std::cout << "input   : " << args_.input() << std::endl;
    std::cout << "output  : " << args_.output() << std::endl;
    std::cout << "end     : " << args_.end() << std::endl;
    std::cout << "progress: " << args_.progress() << std::endl;

    // disable potential multi-threading to preserve the order of entries
    ROOT::DisableImplicitMT();

    // read the input tree
    std::string channelName = EnumNameMap<Channel>::GetDefault().EnumToString(args_.channel());
    TTree* inTree = (TTree*)inputFile_->Get(channelName.c_str());

    // create the input feature provider
    FeatureProvider features(args_.period(), args_.channel(), inTree);

    // create the output tree
    outputFile_->cd();
    TTree* outTree = new TTree(channelName.c_str(), channelName.c_str());

    // load models and define output branches
    std::vector<std::pair<std::string, std::string>> modelSpecs = { { "v0", "kl1_c2v1_c31" } };
    std::vector<hmc::Model*> models;
    for (const auto& modelSpec : modelSpecs) {
      const std::string& version = modelSpec.first;
      const std::string& tag = modelSpec.second;

      // load the model
      models.push_back(hmc::loadModel(int(args_.period()), version, tag));
      hmc::Model*& model = models.back();

      // register required features with the feature provider
      for (const auto& featureName : model->getFeatureNames()) {
        features.add(featureName);
      }

      // define branches per output node
      for (const auto& nodeName : model->getAllNodeNames()) {
        auto branchName = "mdnn__" + version + "__" + tag + "__" + nodeName;
        outTree->Branch(branchName.c_str(), model->output.getOutputAddress(nodeName),
            (branchName + "/F").c_str());
      }
    }

    // start iterating
    int nEntries = args_.end() > 0 ? args_.end() : inTree->GetEntries();
    for (int i = 0; i < nEntries; i++) {
      // load the entry and calculate features
      inTree->GetEntry(i);
      features.calculate();

      // fill features of all models and run them
      for (hmc::Model*& model : models) {
        model->input.clear();
        for (const auto& featureName : model->getFeatureNames()) {
          model->input.setValue(featureName, features.get(featureName));
        }

        model->run(features.getEventId());
      }

      // fill the output tree
      outTree->Fill();

      // print progress
      if (i == 0 || (i + 1) % args_.progress() == 0 || i == nEntries - 1) {
        std::cout << "processed entry " << (i + 1) << std::endl;
      }
    }

    // write the output file
    std::cout << "writing output tree" << std::endl;
    outputFile_->Write();

    // close files and finish
    outputFile_->Close();
    inputFile_->Close();
    std::cout << "done" << std::endl;
  }

  private:
  CalcMulticlassDNNArguments args_;
  std::shared_ptr<TFile> inputFile_;
  std::shared_ptr<TFile> outputFile_;
};

FeatureProvider::FeatureProvider(Period period, Channel channel, TTree* inTree)
    : period_(period)
    , channel_(channel) {
  // define names of variables to read
  std::vector<std::string> boolNames = { "pass_VBF_trigger" };
  std::vector<std::string> intNames = { "b1_valid", "b2_valid", "VBF1_valid", "VBF2_valid" };
  std::vector<std::string> ulong64Names = { "evt" };
  std::vector<std::string> floatNames = { "b1_pt", "b1_eta", "b1_phi", "b1_m", "b1_DeepFlavour",
    "b1_DeepFlavour_CvsB", "b1_DeepFlavour_CvsL", "b1_HHbtag", "b2_pt", "b2_eta", "b2_phi", "b2_m",
    "b2_DeepFlavour", "b2_DeepFlavour_CvsB", "b2_DeepFlavour_CvsL", "b2_HHbtag", "VBF1_pt",
    "VBF1_eta", "VBF1_phi", "VBF1_m", "VBF1_DeepFlavour", "VBF1_DeepFlavour_CvsB",
    "VBF1_DeepFlavour_CvsL", "VBF1_HHbtag", "VBF2_pt", "VBF2_eta", "VBF2_phi", "VBF2_m",
    "VBF2_DeepFlavour", "VBF2_DeepFlavour_CvsB", "VBF2_DeepFlavour_CvsL", "VBF2_HHbtag", "tau1_pt",
    "tau1_eta", "tau1_phi", "tau1_m", "tau2_pt", "tau2_eta", "tau2_phi", "tau2_m", "MET_pt",
    "MET_phi" };

  // register them in input maps and set branch addresses
  for (const auto& name : boolNames) {
    boolInputs_.emplace(name, 0.);
    inTree->SetBranchAddress(name.c_str(), &boolInputs_.at(name));
  }
  for (const auto& name : intNames) {
    intInputs_.emplace(name, 0.);
    inTree->SetBranchAddress(name.c_str(), &intInputs_.at(name));
  }
  for (const auto& name : ulong64Names) {
    ulong64Inputs_.emplace(name, 0);
    inTree->SetBranchAddress(name.c_str(), &ulong64Inputs_.at(name));
  }
  for (const auto& name : floatNames) {
    floatInputs_.emplace(name, 0.);
    inTree->SetBranchAddress(name.c_str(), &floatInputs_.at(name));
  }
}

void FeatureProvider::calculate() {
  // check if objects are set
  bool b1Set = intInputs_.at("b1_valid") == 1;
  bool b2Set = intInputs_.at("b2_valid") == 1;
  bool vbfj1Set = intInputs_.at("VBF1_valid") == 1;
  bool vbfj2Set = intInputs_.at("VBF2_valid") == 1;
  bool lep1Set = floatInputs_.at("tau1_pt") > 0;
  bool lep2Set = floatInputs_.at("tau2_pt") > 0;
  bool bHSet = b1Set && b2Set;
  bool tauHSet = lep1Set && lep2Set;
  bool vbfjjSet = vbfj1Set && vbfj2Set;

  // define vectors for objects
  TLorentzVector b1, b2, vbfj1, vbfj2, lep1, lep2, bH, tauH, vbfjj;
  b1.SetPtEtaPhiM(floatInputs_.at("b1_pt"), floatInputs_.at("b1_eta"),
      floatInputs_.at("b1_phi"), floatInputs_.at("b1_m"));
  b2.SetPtEtaPhiM(floatInputs_.at("b2_pt"), floatInputs_.at("b2_eta"),
      floatInputs_.at("b2_phi"), floatInputs_.at("b2_m"));
  vbfj1.SetPtEtaPhiM(floatInputs_.at("VBF1_pt"), floatInputs_.at("VBF1_eta"),
      floatInputs_.at("VBF1_phi"), floatInputs_.at("VBF1_m"));
  vbfj2.SetPtEtaPhiM(floatInputs_.at("VBF2_pt"), floatInputs_.at("VBF2_eta"),
      floatInputs_.at("VBF2_phi"), floatInputs_.at("VBF2_m"));
  lep1.SetPtEtaPhiM(floatInputs_.at("tau1_pt"), floatInputs_.at("tau1_eta"),
      floatInputs_.at("tau1_phi"), floatInputs_.at("tau1_m"));
  lep2.SetPtEtaPhiM(floatInputs_.at("tau2_pt"), floatInputs_.at("tau2_eta"),
      floatInputs_.at("tau2_phi"), floatInputs_.at("tau2_m"));
  bH = b1 + b2;
  tauH = lep1 + lep2;
  vbfjj = vbfj1 + vbfj2;

  // loop through features and set values
  for (auto& it : features_) {
    if (it.first == "is_mutau") {
      it.second = float(channel_ == Channel::MuTau);
    } else if (it.first == "is_etau") {
      it.second = float(channel_ == Channel::ETau);
    } else if (it.first == "is_tautau") {
      it.second = float(channel_ == Channel::TauTau);
    } else if (it.first == "bjet1_pt") {
      it.second = CHECK_EMPTY(b1Set, b1.Pt());
    } else if (it.first == "bjet1_eta") {
      it.second = CHECK_EMPTY(b1Set, b1.Eta());
    } else if (it.first == "bjet1_phi") {
      it.second = CHECK_EMPTY(b1Set, b1.Phi());
    } else if (it.first == "bjet1_e") {
      it.second = CHECK_EMPTY(b1Set, b1.E());
    } else if (it.first == "bjet1_deepflavor_b") {
      it.second = CHECK_EMPTY(b1Set, floatInputs_.at("b1_DeepFlavour"));
    } else if (it.first == "bjet1_deepflavor_cvsb") {
      it.second = CHECK_EMPTY(b1Set, floatInputs_.at("b1_DeepFlavour_CvsB"));
    } else if (it.first == "bjet1_deepflavor_cvsl") {
      it.second = CHECK_EMPTY(b1Set, floatInputs_.at("b1_DeepFlavour_CvsL"));
    } else if (it.first == "bjet1_hhbtag") {
      it.second = CHECK_EMPTY(b1Set, floatInputs_.at("b1_HHbtag"));
    } else if (it.first == "bjet2_pt") {
      it.second = CHECK_EMPTY(b2Set, b2.Pt());
    } else if (it.first == "bjet2_eta") {
      it.second = CHECK_EMPTY(b2Set, b2.Eta());
    } else if (it.first == "bjet2_phi") {
      it.second = CHECK_EMPTY(b2Set, b2.Phi());
    } else if (it.first == "bjet2_e") {
      it.second = CHECK_EMPTY(b2Set, b2.E());
    } else if (it.first == "bjet2_deepflavor_b") {
      it.second = CHECK_EMPTY(b2Set, floatInputs_.at("b2_DeepFlavour"));
    } else if (it.first == "bjet2_deepflavor_cvsb") {
      it.second = CHECK_EMPTY(b2Set, floatInputs_.at("b2_DeepFlavour_CvsB"));
    } else if (it.first == "bjet2_deepflavor_cvsl") {
      it.second = CHECK_EMPTY(b2Set, floatInputs_.at("b2_DeepFlavour_CvsL"));
    } else if (it.first == "bjet2_hhbtag") {
      it.second = CHECK_EMPTY(b2Set, floatInputs_.at("b2_HHbtag"));
    } else if (it.first == "vbfjet1_pt") {
      it.second = CHECK_EMPTY(vbfj1Set, vbfj1.Pt());
    } else if (it.first == "vbfjet1_eta") {
      it.second = CHECK_EMPTY(vbfj1Set, vbfj1.Eta());
    } else if (it.first == "vbfjet1_phi") {
      it.second = CHECK_EMPTY(vbfj1Set, vbfj1.Phi());
    } else if (it.first == "vbfjet1_e") {
      it.second = CHECK_EMPTY(vbfj1Set, vbfj1.E());
    } else if (it.first == "vbfjet1_deepflavor_b") {
      it.second = CHECK_EMPTY(vbfj1Set, floatInputs_.at("VBF1_DeepFlavour"));
    } else if (it.first == "vbfjet1_deepflavor_cvsb") {
      it.second = CHECK_EMPTY(vbfj1Set, floatInputs_.at("VBF1_DeepFlavour_CvsB"));
    } else if (it.first == "vbfjet1_deepflavor_cvsl") {
      it.second = CHECK_EMPTY(vbfj1Set, floatInputs_.at("VBF1_DeepFlavour_CvsL"));
    } else if (it.first == "vbfjet1_hhbtag") {
      it.second = CHECK_EMPTY(vbfj1Set, floatInputs_.at("VBF1_HHbtag"));
    } else if (it.first == "vbfjet2_pt") {
      it.second = CHECK_EMPTY(vbfj2Set, vbfj2.Pt());
    } else if (it.first == "vbfjet2_eta") {
      it.second = CHECK_EMPTY(vbfj2Set, vbfj2.Eta());
    } else if (it.first == "vbfjet2_phi") {
      it.second = CHECK_EMPTY(vbfj2Set, vbfj2.Phi());
    } else if (it.first == "vbfjet2_e") {
      it.second = CHECK_EMPTY(vbfj2Set, vbfj2.E());
    } else if (it.first == "vbfjet2_deepflavor_b") {
      it.second = CHECK_EMPTY(vbfj2Set, floatInputs_.at("VBF2_DeepFlavour"));
    } else if (it.first == "vbfjet2_deepflavor_cvsb") {
      it.second = CHECK_EMPTY(vbfj2Set, floatInputs_.at("VBF2_DeepFlavour_CvsB"));
    } else if (it.first == "vbfjet2_deepflavor_cvsl") {
      it.second = CHECK_EMPTY(vbfj2Set, floatInputs_.at("VBF2_DeepFlavour_CvsL"));
    } else if (it.first == "vbfjet2_hhbtag") {
      it.second = CHECK_EMPTY(vbfj1Set, floatInputs_.at("VBF2_HHbtag"));
    } else if (it.first == "lep1_pt") {
      it.second = CHECK_EMPTY(lep1Set, lep1.Pt());
    } else if (it.first == "lep1_eta") {
      it.second = CHECK_EMPTY(lep1Set, lep1.Eta());
    } else if (it.first == "lep1_phi") {
      it.second = CHECK_EMPTY(lep1Set, lep1.Phi());
    } else if (it.first == "lep1_e") {
      it.second = CHECK_EMPTY(lep1Set, lep1.E());
    } else if (it.first == "lep2_pt") {
      it.second = CHECK_EMPTY(lep2Set, lep2.Pt());
    } else if (it.first == "lep2_eta") {
      it.second = CHECK_EMPTY(lep2Set, lep2.Eta());
    } else if (it.first == "lep2_phi") {
      it.second = CHECK_EMPTY(lep2Set, lep2.Phi());
    } else if (it.first == "lep2_e") {
      it.second = CHECK_EMPTY(lep2Set, lep2.E());
    } else if (it.first == "met_pt") {
      it.second = floatInputs_.at("MET_pt");
    } else if (it.first == "met_phi") {
      it.second = floatInputs_.at("MET_phi");
    } else if (it.first == "bh_pt") {
      it.second = CHECK_EMPTY(bHSet, bH.Pt());
    } else if (it.first == "bh_eta") {
      it.second = CHECK_EMPTY(bHSet, bH.Eta());
    } else if (it.first == "bh_phi") {
      it.second = CHECK_EMPTY(bHSet, bH.Phi());
    } else if (it.first == "bh_e") {
      it.second = CHECK_EMPTY(bHSet, bH.E());
    } else if (it.first == "tauh_sv_pt") {
      it.second = CHECK_EMPTY(tauHSet, tauH.Pt());
    } else if (it.first == "tauh_sv_eta") {
      it.second = CHECK_EMPTY(tauHSet, tauH.Eta());
    } else if (it.first == "tauh_sv_phi") {
      it.second = CHECK_EMPTY(tauHSet, tauH.Phi());
    } else if (it.first == "tauh_sv_e") {
      it.second = CHECK_EMPTY(tauHSet, tauH.E());
    } else {
      throw exception("MulticlassInference: unhandled feature '" + it.first + "'");
    }
  }
}

// bool FeatureProvider::passBaseline_() const {
//   // currently no distinction between years
//   if (channel_ == Channel::MuTau) {
//     return floatInputs_.at("tau1_pt") > 20. && fabs(floatInputs_.at("tau1_eta")) < 2.1
//         && floatInputs_.at("tau2_pt") > 20. && fabs(floatInputs_.at("tau1_eta")) < 2.3
//         && boolInputs_.at("pass_VBF_trigger");
//   } else if (channel_ == Channel::ETau) {
//     return floatInputs_.at("tau1_pt") > 20. && fabs(floatInputs_.at("tau1_eta")) < 2.1
//         && floatInputs_.at("tau2_pt") > 20. && fabs(floatInputs_.at("tau1_eta")) < 2.3
//         && boolInputs_.at("pass_VBF_trigger");
//   } else if (channel_ == Channel::TauTau) {
//     return floatInputs_.at("tau1_pt") > 40. && fabs(floatInputs_.at("tau1_eta")) < 2.1
//         && floatInputs_.at("tau2_pt") > 40. && fabs(floatInputs_.at("tau1_eta")) < 2.1
//         && boolInputs_.at("pass_VBF_trigger");
//   } else {
//     throw exception("unknown channel " + std::to_string(channel_));
//   }
// }

// bool FeatureProvider::passEllipseMassCut_(
//     const TLorentzVector& bH, const TLorentzVector& tauH) const {
//   auto bHPull = (bH.M() - 111.) / 45.;
//   auto tauHPull = (tauH.M() - 116.) / 35.;
//   return (bHPull * bHPull + tauHPull * tauHPull) < 1.;
// }

// bool FeatureProvider::passRectMassCut_(const TLorentzVector& bH, const TLorentzVector& tauH)
// const {
//   // TODO: use fatjet softdrop mass instead of bH mass?
//   return bH.M() > 90. && bH.M() < 160. && tauH.M() > 79.5 && tauH.M() < 152.5;
// }

} // namespace analysis

PROGRAM_MAIN(analysis::CalcMulticlassDNN, analysis::CalcMulticlassDNNArguments)
