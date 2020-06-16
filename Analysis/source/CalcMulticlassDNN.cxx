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
  FeatureProvider(TTree* inTree)
      : inTree_(inTree) {
    // define variables to read
    ulong64Inputs_.emplace("evt", 0);
    floatInputs_.emplace("tau1_pt", 0.);

    // set branch adresses
    for (auto& it : floatInputs_) {
      inTree->SetBranchAddress(it.first.c_str(), &it.second);
    }
    for (auto& it : boolInputs_) {
      inTree->SetBranchAddress(it.first.c_str(), &it.second);
    }
    for (auto& it : ulong64Inputs_) {
      inTree->SetBranchAddress(it.first.c_str(), &it.second);
    }
  }

  FeatureProvider(const FeatureProvider&) = delete;

  inline hmc::EventId getEventId() const {
    return hmc::EventId(ulong64Inputs_.at("evt"));
  }

  inline void add(const std::string& featureName) {
    features_.emplace(featureName, 0.);
  }

  void calculate(int i) {
    inTree_->GetEntry(i);

    for (auto& it : features_) {
      // TODO
      it.second = 0.5;
    }
  }

  float get(const std::string& featureName) const {
    const auto& it = features_.find(featureName);
    if (it == features_.end()) {
      throw exception("FeatureProvider: unknown feature '" + featureName + "'");
    }
    return it->second;
  }

  private:
  TTree* inTree_;
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
    std::cout << "channel: " << args_.channel() << std::endl;
    std::cout << "period : " << args_.period() << std::endl;
    std::cout << "input  : " << args_.input() << std::endl;
    std::cout << "output : " << args_.output() << std::endl;

    // disable potential multi-threading to preserve the order of entries
    ROOT::DisableImplicitMT();

    // read the input tree
    std::string channelName = EnumNameMap<Channel>::GetDefault().EnumToString(args_.channel());
    TTree* inTree = (TTree*)inputFile_->Get(channelName.c_str());

    // create the input feature provider
    FeatureProvider features(inTree);

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

      // propagate required features to the feature provider
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
      // calculate features for this entry
      features.calculate(i);

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

} // namespace analysis

PROGRAM_MAIN(analysis::CalcMulticlassDNN, analysis::CalcMulticlassDNNArguments)
