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
  REQ_ARG(std::string, model);
  OPT_ARG(Long64_t, max_events, std::numeric_limits<Long64_t>::max());
  OPT_ARG(Long64_t, begin_entry_index, 0);
  OPT_ARG(Long64_t, end_entry_index, std::numeric_limits<Long64_t>::max());
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
};

class CalcMulticlassDNN {
  public:
  CalcMulticlassDNN(const CalcMulticlassDNNArguments& args)
      : args_(args)
      , inputFile_(root_ext::OpenRootFile(args_.input()))
      , outputFile_(root_ext::CreateRootFile(args_.output(), ROOT::kLZMA, 9)) {
  }

  void Run() {
    std::cout << "running multiclass classification inference\n"
              << "channel           : " << args_.channel() << '\n'
              << "period            : " << args_.period() << '\n'
              << "input             : " << args_.input() << '\n'
              << "output            : " << args_.output() << '\n'
              << "max_events        : " << args_.max_events() << '\n'
              << "begin_entry_index : " << args_.begin_entry_index() << '\n'
              << "end_entry_index   : " << args_.end_entry_index() << '\n'
              << "progress          : " << args_.progress() << std::endl;

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
    static const std::vector<std::pair<std::string, std::string>> modelSpecs = {
      { "v3", "kl1_c2v1_c31_vbf" },
      { "v3", "kl1_c2v1_c31_vr" },
      { "v4", "kl1_c2v1_c31_vbf" },
      { "v4", "kl1_c2v1_c31_vr" },
      { "v5", "kl1_c2v1_c31_vbf" },
    };
    const auto activeModels = SplitValueListT<std::string, std::set<std::string>>(args_.model(), false, ",", true);

    std::vector<hmc::Model*> models;
    for (const auto& [version, tag] : modelSpecs) {
      if(!activeModels.count(version)) continue;

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
    const Long64_t endEntry = std::min(anaTuple.GetEntries(), args_.end_entry_index());
    const Long64_t nEntries = std::min(endEntry - args_.begin_entry_index(), args_.max_events());
    tools::ProgressReporter progressReporter(10, std::cout);
    progressReporter.SetTotalNumberOfEvents(static_cast<unsigned long>(nEntries));

    for (Long64_t i = 0; i < nEntries; ++i) {
      const Long64_t entry_index = args_.begin_entry_index() + i;

      // calculate features
      features.calculate(entry_index);

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
  bool ct1Set = event.central_jet1_valid;
  bool ct2Set = event.central_jet2_valid;
  bool ct3Set = event.central_jet3_valid;
  bool fw1Set = event.forward_jet1_valid;
  bool fw2Set = event.forward_jet2_valid;

  // define vectors for objects
  TLorentzVector b1, b2, vbfj1, vbfj2, lep1, lep2, bH, tauH, vbfjj, ct1, ct2, ct3, fw1, fw2;
  b1.SetPtEtaPhiM(event.b1_pt, event.b1_eta, event.b1_phi, event.b1_m);
  b2.SetPtEtaPhiM(event.b2_pt, event.b2_eta, event.b2_phi, event.b2_m);
  vbfj1.SetPtEtaPhiM(event.VBF1_pt, event.VBF1_eta, event.VBF1_phi, event.VBF1_m);
  vbfj2.SetPtEtaPhiM(event.VBF2_pt, event.VBF2_eta, event.VBF2_phi, event.VBF2_m);
  lep1.SetPtEtaPhiM(event.tau1_pt, event.tau1_eta, event.tau1_phi, event.tau1_m);
  lep2.SetPtEtaPhiM(event.tau2_pt, event.tau2_eta, event.tau2_phi, event.tau2_m);
  ct1.SetPtEtaPhiM(event.central_jet1_pt, event.central_jet1_eta, event.central_jet1_phi, event.central_jet1_m);
  ct2.SetPtEtaPhiM(event.central_jet2_pt, event.central_jet2_eta, event.central_jet2_phi, event.central_jet2_m);
  ct3.SetPtEtaPhiM(event.central_jet3_pt, event.central_jet3_eta, event.central_jet3_phi, event.central_jet3_m);
  fw1.SetPtEtaPhiM(event.forward_jet1_pt, event.forward_jet1_eta, event.forward_jet1_phi, event.forward_jet1_m);
  fw2.SetPtEtaPhiM(event.forward_jet2_pt, event.forward_jet2_eta, event.forward_jet2_phi, event.forward_jet2_m);
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
  setFeature("is_2016", float(period_ == Period::Run2016), true);
  setFeature("is_2017", float(period_ == Period::Run2017), true);
  setFeature("is_2018", float(period_ == Period::Run2018), true);
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
  setFeature("ctjet1_pt", ct1.Pt(), ct1Set);
  setFeature("ctjet1_eta", ct1.Eta(), ct1Set);
  setFeature("ctjet1_phi", ct1.Phi(), ct1Set);
  setFeature("ctjet1_e", ct1.E(), ct1Set);
  setFeature("ctjet1_deepflavor_b", event.central_jet1_DeepFlavour, ct1Set);
  setFeature("ctjet1_hhbtag", event.central_jet1_HHbtag, ct1Set);
  setFeature("ctjet2_pt", ct2.Pt(), ct2Set);
  setFeature("ctjet2_eta", ct2.Eta(), ct2Set);
  setFeature("ctjet2_phi", ct2.Phi(), ct2Set);
  setFeature("ctjet2_e", ct2.E(), ct2Set);
  setFeature("ctjet2_deepflavor_b", event.central_jet2_DeepFlavour, ct2Set);
  setFeature("ctjet2_hhbtag", event.central_jet2_HHbtag, ct2Set);
  setFeature("ctjet3_pt", ct3.Pt(), ct3Set);
  setFeature("ctjet3_eta", ct3.Eta(), ct3Set);
  setFeature("ctjet3_phi", ct3.Phi(), ct3Set);
  setFeature("ctjet3_e", ct3.E(), ct3Set);
  setFeature("ctjet3_deepflavor_b", event.central_jet3_DeepFlavour, ct3Set);
  setFeature("ctjet3_hhbtag", event.central_jet3_HHbtag, ct3Set);
  setFeature("fwjet1_pt", fw1.Pt(), fw1Set);
  setFeature("fwjet1_eta", fw1.Eta(), fw1Set);
  setFeature("fwjet1_phi", fw1.Phi(), fw1Set);
  setFeature("fwjet1_e", fw1.E(), fw1Set);
  setFeature("fwjet2_pt", fw2.Pt(), fw2Set);
  setFeature("fwjet2_eta", fw2.Eta(), fw2Set);
  setFeature("fwjet2_phi", fw2.Phi(), fw2Set);
  setFeature("fwjet2_e", fw2.E(), fw2Set);
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
  setFeature("tauh_sv_ez", pow(pow(tauH.Pz(), 2) + pow(tauH.M(), 2), 0.5), true);

  if (setCounter != features_.size()) {
    throw exception("only calculated " + std::to_string(setCounter) + " out of "
        + std::to_string(features_.size()) + " features");
  }
}

} // namespace analysis

PROGRAM_MAIN(analysis::CalcMulticlassDNN, analysis::CalcMulticlassDNNArguments)
