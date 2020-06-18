/*! Writes multiclass outputs to a new tree that can be used as a friend. */

#include "AnalysisTools/Core/include/ProgressReporter.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "hh-bbtautau/Analysis/include/AnaTuple.h"

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
  FeatureProvider(Period period, Channel channel, bbtautau::AnaTuple& anaTuple)
      : period_(period)
      , channel_(channel)
      , anaTuple_(anaTuple) {
  }

  FeatureProvider(const FeatureProvider&) = delete;

  void calculate(Long64_t i);

  inline hmc::EventId getEventId() const {
    return hmc::EventId(anaTuple_().evt);
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
  bbtautau::AnaTuple& anaTuple_;
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

    // read the input tree as an AnaTuple
    std::string channelName = EnumNameMap<Channel>::GetDefault().EnumToString(args_.channel());
    bbtautau::AnaTuple anaTuple(channelName, inputFile_.get(), true);

    // create the input feature provider
    FeatureProvider features(args_.period(), args_.channel(), anaTuple);

    // create the output tree
    outputFile_->cd();
    auto outTree = std::make_unique<TTree>(channelName.c_str(), channelName.c_str());

    // load models and define output branches
    std::vector<std::pair<std::string, std::string>> modelSpecs
        = { { "v0", "kl1_c2v1_c31" }, { "v0", "kl1_c2v1_c31_vbfbsm" } };
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
    const Long64_t nEntries = args_.end() > 0 ? args_.end() : anaTuple.GetEntries();
    tools::ProgressReporter progressReporter(10, std::cout);
    progressReporter.SetTotalNumberOfEvents(nEntries);
    for (Long64_t i = 0; i < nEntries; i++) {
      // calculate features
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
      if (i % args_.progress() == 0) {
        progressReporter.Report(i + 1, false);
      }
    }
    progressReporter.Report(nEntries, true);

    // write the output file
    std::cout << "writing output tree" << std::endl;
    outputFile_->Write();

    // close files and finish
    outputFile_->Close();
    inputFile_->Close();
    outTree.release();
    std::cout << "done" << std::endl;
  }

  private:
  CalcMulticlassDNNArguments args_;
  std::shared_ptr<TFile> inputFile_;
  std::shared_ptr<TFile> outputFile_;
};

void FeatureProvider::calculate(Long64_t i) {
  anaTuple_.GetEntry(i);
  const bbtautau::AnaEvent& event = anaTuple_.data();

  // check if objects are set
  bool b1Set = event.b1_valid == 1;
  bool b2Set = event.b2_valid == 1;
  bool vbfj1Set = event.VBF1_valid == 1;
  bool vbfj2Set = event.VBF2_valid == 1;
  bool bHSet = b1Set && b2Set;
  bool vbfjjSet = vbfj1Set && vbfj2Set;

  // define vectors for objects
  TLorentzVector b1, b2, vbfj1, vbfj2, lep1, lep2, bH, tauH, vbfjj;
  b1.SetPtEtaPhiM(event.b1_pt, event.b1_eta, event.b1_phi, event.b1_m);
  b2.SetPtEtaPhiM(event.b2_pt, event.b2_eta, event.b2_phi, event.b2_m);
  vbfj1.SetPtEtaPhiM(event.VBF1_pt, event.VBF1_eta, event.VBF1_phi, event.VBF1_m);
  vbfj2.SetPtEtaPhiM(event.VBF2_pt, event.VBF2_eta, event.VBF2_phi, event.VBF2_m);
  lep1.SetPtEtaPhiM(event.tau1_pt, event.tau1_eta, event.tau1_phi, event.tau1_m);
  lep2.SetPtEtaPhiM(event.tau2_pt, event.tau2_eta, event.tau2_phi, event.tau2_m);
  bH = b1 + b2;
  tauH = lep1 + lep2;
  vbfjj = vbfj1 + vbfj2;

  // set features
  size_t setCounter = 0;
  const auto setFeature = [&](const std::string& name, float value, bool condition) {
    auto iter = features_.find(name);
    if (iter != features_.end()) {
      iter->second = condition ? value : hmc::features::EMPTY;
      ++setCounter;
    }
  };

  setFeature("is_mutau", float(channel_ == Channel::MuTau), true);
  setFeature("is_etau", float(channel_ == Channel::ETau), true);
  setFeature("is_tautau", float(channel_ == Channel::TauTau), true);
  setFeature("bjet1_pt", b1.Pt(), b1Set);
  setFeature("bjet1_eta", b1.Eta(), b1Set);
  setFeature("bjet1_phi", b1.Phi(), b1Set);
  setFeature("bjet1_e", b1.E(), b1Set);
  setFeature("bjet1_deepflavor_b", event.b1_DeepFlavour, b1Set);
  setFeature("bjet1_deepflavor_cvsb", event.b1_DeepFlavour_CvsB, b1Set);
  setFeature("bjet1_deepflavor_cvsl", event.b1_DeepFlavour_CvsL, b1Set);
  setFeature("bjet1_hhbtag", event.b1_HHbtag, b1Set);
  setFeature("bjet2_pt", b2.Pt(), b2Set);
  setFeature("bjet2_eta", b2.Eta(), b2Set);
  setFeature("bjet2_phi", b2.Phi(), b2Set);
  setFeature("bjet2_e", b2.E(), b2Set);
  setFeature("bjet2_deepflavor_b", event.b2_DeepFlavour, b2Set);
  setFeature("bjet2_deepflavor_cvsb", event.b2_DeepFlavour_CvsB, b2Set);
  setFeature("bjet2_deepflavor_cvsl", event.b2_DeepFlavour_CvsL, b2Set);
  setFeature("bjet2_hhbtag", event.b2_HHbtag, b2Set);
  setFeature("vbfjet1_pt", vbfj1.Pt(), vbfj1Set);
  setFeature("vbfjet1_eta", vbfj1.Eta(), vbfj1Set);
  setFeature("vbfjet1_phi", vbfj1.Phi(), vbfj1Set);
  setFeature("vbfjet1_e", vbfj1.E(), vbfj1Set);
  setFeature("vbfjet1_deepflavor_b", event.VBF1_DeepFlavour, vbfj1Set);
  setFeature("vbfjet1_deepflavor_cvsb", event.VBF1_DeepFlavour_CvsB, vbfj1Set);
  setFeature("vbfjet1_deepflavor_cvsl", event.VBF1_DeepFlavour_CvsL, vbfj1Set);
  setFeature("vbfjet1_hhbtag", event.VBF1_HHbtag, vbfj1Set);
  setFeature("vbfjet2_pt", vbfj2.Pt(), vbfj2Set);
  setFeature("vbfjet2_eta", vbfj2.Eta(), vbfj2Set);
  setFeature("vbfjet2_phi", vbfj2.Phi(), vbfj2Set);
  setFeature("vbfjet2_e", vbfj2.E(), vbfj2Set);
  setFeature("vbfjet2_deepflavor_b", event.VBF2_DeepFlavour, vbfj2Set);
  setFeature("vbfjet2_deepflavor_cvsb", event.VBF2_DeepFlavour_CvsB, vbfj2Set);
  setFeature("vbfjet2_deepflavor_cvsl", event.VBF2_DeepFlavour_CvsL, vbfj2Set);
  setFeature("vbfjet2_hhbtag", event.VBF2_HHbtag, vbfj1Set);
  setFeature("lep1_pt", lep1.Pt(), true);
  setFeature("lep1_eta", lep1.Eta(), true);
  setFeature("lep1_phi", lep1.Phi(), true);
  setFeature("lep1_e", lep1.E(), true);
  setFeature("lep2_pt", lep2.Pt(), true);
  setFeature("lep2_eta", lep2.Eta(), true);
  setFeature("lep2_phi", lep2.Phi(), true);
  setFeature("lep2_e", lep2.E(), true);
  setFeature("met_pt", event.MET_pt, true);
  setFeature("met_phi", event.MET_phi, true);
  setFeature("bh_pt", bH.Pt(), bHSet);
  setFeature("bh_eta", bH.Eta(), bHSet);
  setFeature("bh_phi", bH.Phi(), bHSet);
  setFeature("bh_e", bH.E(), bHSet);
  setFeature("tauh_sv_pt", tauH.Pt(), true);
  setFeature("tauh_sv_eta", tauH.Eta(), true);
  setFeature("tauh_sv_phi", tauH.Phi(), true);
  setFeature("tauh_sv_e", tauH.E(), true);

  if (setCounter != features_.size()) {
    throw exception("only calculated " + std::to_string(setCounter) + " out of "
        + std::to_string(features_.size()) + " features");
  }
}

// bool FeatureProvider::passBaseline_() const {
//   // currently no distinction between years
//   if (channel_ == Channel::MuTau) {
//     return event.tau1_pt > 20. && fabs(event.tau1_eta) < 2.1
//         && event.tau2_pt > 20. && fabs(event.tau1_eta) < 2.3
//         && event.pass_VBF_trigger;
//   } else if (channel_ == Channel::ETau) {
//     return event.tau1_pt > 20. && fabs(event.tau1_eta) < 2.1
//         && event.tau2_pt > 20. && fabs(event.tau1_eta) < 2.3
//         && event.pass_VBF_trigger;
//   } else if (channel_ == Channel::TauTau) {
//     return event.tau1_pt > 40. && fabs(event.tau1_eta) < 2.1
//         && event.tau2_pt > 40. && fabs(event.tau1_eta) < 2.1
//         && event.pass_VBF_trigger;
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
