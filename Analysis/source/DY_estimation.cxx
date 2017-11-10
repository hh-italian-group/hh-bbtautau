/*! This code is to estimate scale factor for DY normalization using the muMu channel in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"

//Root and Roofit Headers
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_file); // required argument "input_file"
    REQ_ARG(std::string, output_file); // required argument "output_file"
    OPT_ARG(bool, flag, false); // optional argument "flag" with the default value = false
};
using namespace RooFit;
namespace analysis {
class DY_estimation { // simple analyzer definition
public:
    DY_estimation(const Arguments& _args) : args(_args),
    input_file(root_ext::OpenRootFile(args.input_file())),
    output_file(root_ext::CreateRootFile(args.output_file()))
    {
        // Analyzer initialization (e.g. open input/output files, parse configs...)
    }
    void Run()
    {
        // analyzer code
        std::cout << boost::format("Processing input file '%1%' into output file '%2%' wiht flag = %3%.\n")
                     % args.input_file() % args.output_file() % args.flag();
        RooRealVar mass("mass","mass",0,500);
        TH1D* dataHisto0b = getTH1D("2jets0btagR/mh/OS_Isolated/Central/Data_SingleMuon");
        //RooDataHist data0b("data0b","data for 0b category",mass,Import(*dataHisto0b));
    }
private:
    Arguments args;
    std::shared_ptr<TFile> input_file, output_file;

    TH1D* getTH1D(std::string path){
        TH1D* h1;
        input_file->GetObject(path.c_str(),h1);
        return h1;
    }


};

} // namesapce analysis
PROGRAM_MAIN(analysis::DY_estimation, Arguments) // definition of the main program function

