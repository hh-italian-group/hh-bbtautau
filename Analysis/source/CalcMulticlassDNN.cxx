/*! Writes multiclass outputs to a new tree that can be used as a friend. */

#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventInfo.h"

#include "MulticlassInference/MulticlassInference/interface/hmc.h"

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
  std::map<std::string, float> floatInputs_;
  std::map<std::string, bool> boolInputs_;
  std::map<std::string, ULong64_t> ulong64Inputs_;
  std::map<std::string, float> features_;
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
      hmc::Model* model = models.back();

      // register required features with the feature provider
      for (const auto& featureName : model->getFeatureNames()) {
        features.add(featureName);
      }

      // define branches per output node
      for (const auto& nodeName : model->getNodeNames()) {
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
  // TODO
  std::vector<std::string> ulong64Names = { "evt" };
  std::vector<std::string> floatNames = { "tau1_pt" };
  std::vector<std::string> boolNames = {};

  // register them in inputs and set branch addresses
  for (const auto& name : ulong64Names) {
    ulong64Inputs_.emplace(name, 0);
    inTree->SetBranchAddress(name.c_str(), &ulong64Inputs_.at(name));
  }
  for (const auto& name : floatNames) {
    floatInputs_.emplace(name, 0.);
    inTree->SetBranchAddress(name.c_str(), &floatInputs_.at(name));
  }
  for (const auto& name : boolNames) {
    boolInputs_.emplace(name, 0.);
    inTree->SetBranchAddress(name.c_str(), &boolInputs_.at(name));
  }
}

void FeatureProvider::calculate() {
  for (auto& it : features_) {
    if (it.first == "is_mutau") {
      it.second = float(channel_ == Channel::MuTau);
    } else if (it.first == "is_etau") {
      it.second = float(channel_ == Channel::ETau);
    } else if (it.first == "is_tautau") {
      it.second = float(channel_ == Channel::TauTau);
    } else if (it.first == "n_jets_20") {
      it.second = 0.5;
    } else if (it.first == "n_bjets_20") {
      it.second = 0.5;
    } else if (it.first == "is_vbf_loose_cat") {
      it.second = 0.5;
    } else if (it.first == "is_vbf_tight_cat") {
      it.second = 0.5;
    } else if (it.first == "is_resolved_1b_cat") {
      it.second = 0.5;
    } else if (it.first == "is_resolved_2b_cat") {
      it.second = 0.5;
    } else if (it.first == "is_boosted_cat") {
      it.second = 0.5;
    } else if (it.first == "bjet1_pt") {
      it.second = 0.5;
    } else if (it.first == "bjet1_eta") {
      it.second = 0.5;
    } else if (it.first == "bjet1_phi") {
      it.second = 0.5;
    } else if (it.first == "bjet1_e") {
      it.second = 0.5;
    } else if (it.first == "bjet1_deepflavor_b") {
      it.second = 0.5;
    } else if (it.first == "bjet1_deepflavor_c") {
      it.second = 0.5;
    } else if (it.first == "bjet2_pt") {
      it.second = 0.5;
    } else if (it.first == "bjet2_eta") {
      it.second = 0.5;
    } else if (it.first == "bjet2_phi") {
      it.second = 0.5;
    } else if (it.first == "bjet2_e") {
      it.second = 0.5;
    } else if (it.first == "bjet2_deepflavor_b") {
      it.second = 0.5;
    } else if (it.first == "bjet2_deepflavor_c") {
      it.second = 0.5;
    } else if (it.first == "jet3_pt") {
      it.second = 0.5;
    } else if (it.first == "jet3_eta") {
      it.second = 0.5;
    } else if (it.first == "jet3_phi") {
      it.second = 0.5;
    } else if (it.first == "jet3_e") {
      it.second = 0.5;
    } else if (it.first == "jet3_deepflavor_b") {
      it.second = 0.5;
    } else if (it.first == "jet3_deepflavor_c") {
      it.second = 0.5;
    } else if (it.first == "jet4_pt") {
      it.second = 0.5;
    } else if (it.first == "jet4_eta") {
      it.second = 0.5;
    } else if (it.first == "jet4_phi") {
      it.second = 0.5;
    } else if (it.first == "jet4_e") {
      it.second = 0.5;
    } else if (it.first == "jet4_deepflavor_b") {
      it.second = 0.5;
    } else if (it.first == "jet4_deepflavor_c") {
      it.second = 0.5;
    } else if (it.first == "jet5_pt") {
      it.second = 0.5;
    } else if (it.first == "jet5_eta") {
      it.second = 0.5;
    } else if (it.first == "jet5_phi") {
      it.second = 0.5;
    } else if (it.first == "jet5_e") {
      it.second = 0.5;
    } else if (it.first == "jet5_deepflavor_b") {
      it.second = 0.5;
    } else if (it.first == "jet5_deepflavor_c") {
      it.second = 0.5;
    } else if (it.first == "vbfjet1_pt") {
      it.second = 0.5;
    } else if (it.first == "vbfjet1_eta") {
      it.second = 0.5;
    } else if (it.first == "vbfjet1_phi") {
      it.second = 0.5;
    } else if (it.first == "vbfjet1_e") {
      it.second = 0.5;
    } else if (it.first == "vbfjet1_deepflavor_b") {
      it.second = 0.5;
    } else if (it.first == "vbfjet1_deepflavor_c") {
      it.second = 0.5;
    } else if (it.first == "vbfjet2_pt") {
      it.second = 0.5;
    } else if (it.first == "vbfjet2_eta") {
      it.second = 0.5;
    } else if (it.first == "vbfjet2_phi") {
      it.second = 0.5;
    } else if (it.first == "vbfjet2_e") {
      it.second = 0.5;
    } else if (it.first == "vbfjet2_deepflavor_b") {
      it.second = 0.5;
    } else if (it.first == "vbfjet2_deepflavor_c") {
      it.second = 0.5;
    } else if (it.first == "lep1_pt") {
      it.second = 0.5;
    } else if (it.first == "lep1_eta") {
      it.second = 0.5;
    } else if (it.first == "lep1_phi") {
      it.second = 0.5;
    } else if (it.first == "lep1_e") {
      it.second = 0.5;
    } else if (it.first == "lep2_pt") {
      it.second = 0.5;
    } else if (it.first == "lep2_eta") {
      it.second = 0.5;
    } else if (it.first == "lep2_phi") {
      it.second = 0.5;
    } else if (it.first == "lep2_e") {
      it.second = 0.5;
    } else if (it.first == "met_pt") {
      it.second = 0.5;
    } else if (it.first == "met_phi") {
      it.second = 0.5;
    } else if (it.first == "bh_pt") {
      it.second = 0.5;
    } else if (it.first == "bh_eta") {
      it.second = 0.5;
    } else if (it.first == "bh_phi") {
      it.second = 0.5;
    } else if (it.first == "bh_e") {
      it.second = 0.5;
    } else if (it.first == "tauh_sv_pt") {
      it.second = 0.5;
    } else if (it.first == "tauh_sv_eta") {
      it.second = 0.5;
    } else if (it.first == "tauh_sv_phi") {
      it.second = 0.5;
    } else if (it.first == "tauh_sv_e") {
      it.second = 0.5;
    } else {
      throw exception("MulticlassInference: undhandled feature '" + it.first + "'");
    }
  }
}

} // namespace analysis

PROGRAM_MAIN(analysis::CalcMulticlassDNN, analysis::CalcMulticlassDNNArguments)
