/*! Adds multiclass outputs into a new tree that can be used as a friend */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"

#include "MulticlassInference/MulticlassInference/interface/hmc.h"

namespace analysis {

struct CalcMulticlassDNNArguments {
  REQ_ARG(Channel, channel);
  REQ_ARG(Period, period);
  REQ_ARG(std::string, input);
  REQ_ARG(std::string, output);
};

class CalcMulticlassDNN {
public:
  CalcMulticlassDNN(const CalcMulticlassDNNArguments &args)
      : args_(args),
        inputFile_(root_ext::OpenRootFile(args_.input())),
        outputFile_(root_ext::CreateRootFile(args_.output())) {}

  void Run() {
    std::cout << "running multiclass classification inference" << std::endl;
    std::cout << "channel: " << args_.channel() << std::endl;
    std::cout << "period : " << args_.period() << std::endl;
    std::cout << "input  : " << args_.input() << std::endl;
    std::cout << "output : " << args_.output() << std::endl;

    // open the input tree
    std::string channelName = EnumNameMap<Channel>::GetDefault().EnumToString(args_.channel());
    TTree* inTree = (TTree*)inputFile_->Get(channelName.c_str());

    // create the output tree
    outputFile_->cd();
    TTree* outTree = new TTree(channelName.c_str(), channelName.c_str());

    // TODO: setup and run inference here
    // test
    hmc::Model* model = hmc::loadModel(2018, "v0", "kl1_c2v1_c31");
    model->input.clear();
    for (int f = 0; f < (int)model->getNumberOfFeatures(); f++) {
        std::string featureName = model->getFeatureName(f);
        model->input.setValue(featureName, (float)f);
    }
    model->run(12);
    for (const auto& it : model->output) {
        std::cout << "Node " << it.first << ", value " << it.second << std::endl;
    }
    // test end

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
