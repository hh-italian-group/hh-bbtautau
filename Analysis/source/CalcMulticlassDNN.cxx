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
    // test
    std::cout << args_.channel() << std::endl;
    std::cout << args_.period() << std::endl;
    std::cout << args_.input() << std::endl;
    std::cout << args_.output() << std::endl;

    std::string channelName = EnumNameMap<Channel>::GetDefault().EnumToString(args_.channel());

    // open the input tree
    TTree* inTree = (TTree*)inputFile_->Get(channelName.c_str());

    // create the output tree
    outputFile_->cd();
    TTree* outTree = new TTree(channelName.c_str(), channelName.c_str());

    // write the output file
    outputFile_->Write();
    outputFile_->Close();
  }

private:
  CalcMulticlassDNNArguments args_;
  std::shared_ptr<TFile> inputFile_;
  std::shared_ptr<TFile> outputFile_;
};

} // namespace analysis

PROGRAM_MAIN(analysis::CalcMulticlassDNN, analysis::CalcMulticlassDNNArguments)
