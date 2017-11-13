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
#include "RooHistPdf.h"
#include "RooExtendPdf.h"
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
        //Initialize the x axis
        RooRealVar mass("mass","x axis for the fits",0,500);
        //mass.setRange(0,500);
        std::cout<<"x axis is initialized"<<std::endl;

        //Initialize the scale factors
        RooRealVar sf_0b("sf_0b","Scale factor for 0b contibution",1,0.1,10);
        RooRealVar sf_1b("sf_1b","Scale factor for 1b contibution",1,0.1,10);
        RooRealVar sf_2b("sf_2b","Scale factor for 2b contibution",1,0.1,10);
        RooRealVar sf_ob("sf_ob","Scale factor for other bkg contibution",1,0.1,10);
        std::cout<<"Variables are initialized"<<std::endl;


        //Data Histograms
        // For 0b category
        TH1D* dataHisto0b = dynamic_cast<TH1D*>(input_file->Get("2jets0btagR/mh/OS_Isolated/Central/Data_SingleMuon"));
        RooDataHist data0b("data0b","data for 0b category",mass,Import(*dataHisto0b));
        // For 1b Category
        TH1D* dataHisto1b = dynamic_cast<TH1D*>(input_file->Get("2jets1btagR/mh/OS_Isolated/Central/Data_SingleMuon"));
        RooDataHist data1b("data1b","data for 1b category",mass,Import(*dataHisto1b));
        // For 2b Category
        TH1D* dataHisto2b = dynamic_cast<TH1D*>(input_file->Get("2jets2btagR/mh/OS_Isolated/Central/Data_SingleMuon"));
        RooDataHist data2b("data1b","data for 2b category",mass,Import(*dataHisto2b));


        //Create the Monte-Carlo Pdfs

        //<<<<<<<<<<<<<<<<<<< For 0b Category >>>>>>>>>>>>>>>>>>>>>>>>>>>>
        // 0b contribution
        TH1D* mchist_0b_0b = dynamic_cast<TH1D*>(input_file->Get("2jets0btagR/mh/OS_Isolated/Central/DY_0b/m_tt_vis"));
        //RooRealVar norm_0b_0b("norm_0b_0b","norm_0b_0b",mchist_0b_0b->Integral());
        //RooFormulaVar frac_0b_0b("frac_0b_0b","frac_0b_0b","@0*@1",RooArgList(sf_0b,norm_0b_0b));
        //mchist_0b_0b->Scale(1./mchist_0b_0b->Integral());
        //RooDataHist DY_0b_0b("DY_0b_0b","DY_0b_0b",mass,Import(*mchist_0b_0b)) ;
        //RooHistPdf DY_0b_0b_pdf("DY_0b_0b_pdf","DY_0b_0b_pdf",mass,DY_0b_0b);
        RooHistPdf DY_0b_0b_pdf = getPdf("DY_0b_0b","for 0b contribution in 0b category",
                                                    mchist_0b_0b,mass,sf_0b);
        /*// 1b contribution
        TH1D* mchist_1b_0b = dynamic_cast<TH1D*>(input_file->Get("2jets0btagR/mh/OS_Isolated/Central/DY_1b/m_tt_vis"));
        RooHistPdf DY_1b_0b_pdf = getPdf("DY_1b_0b","for 1b contribution in 0b category",
                                                     mchist_1b_0b,mass,sf_1b);
        // 2b contribution
        TH1D* mchist_2b_0b = dynamic_cast<TH1D*>(input_file->Get("2jets0btagR/mh/OS_Isolated/Central/DY_2b/m_tt_vis"));
        RooHistPdf DY_2b_0b_pdf = getPdf("DY_2b_0b","for 2b contribution in 0b category",
                                                     mchist_2b_0b,mass,sf_2b);
        // Other Background contribution
        TH1D* mchist_ob_0b = dynamic_cast<TH1D*>(input_file->Get("2jets0btagR/mh/OS_Isolated/Central/other_bkg/m_tt_vis"));
        RooHistPdf DY_ob_0b_pdf = getPdf("DY_ob_0b","for other bkg contibution in 0b category",
                                                     mchist_ob_0b,mass,sf_ob);

        // Adding all contributions and create the total pdf
        RooAddPdf mcpdf_0b("mcpdf_0b","0b+1b+2b+other_bkg for 0 btag category",
                           RooArgList(DY_0b_0b_pdf,DY_1b_0b_pdf,DY_2b_0b_pdf,DY_ob_0b_pdf));

        //<<<<<<<<<<<<<<<< For 1b Category >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        // 0b contribution
        TH1D* mchist_0b_1b = dynamic_cast<TH1D*>(input_file->Get("2jets1btagR/mh/OS_Isolated/Central/DY_0b/m_tt_vis"));
        RooHistPdf DY_0b_1b_pdf = getPdf("DY_0b_1b","for 0b contribution in 1b category",
                                                     mchist_0b_1b,mass,sf_0b);
        // 1b contribution
        TH1D* mchist_1b_1b = dynamic_cast<TH1D*>(input_file->Get("2jets1btagR/mh/OS_Isolated/Central/DY_1b/m_tt_vis"));
        RooHistPdf DY_1b_1b_pdf = getPdf("DY_1b_1b","for 1b contribution in 1b category",
                                                     mchist_1b_1b,mass,sf_1b);
        // 2b contribution
        TH1D* mchist_2b_1b = dynamic_cast<TH1D*>(input_file->Get("2jets1btagR/mh/OS_Isolated/Central/DY_2b/m_tt_vis"));
        RooHistPdf DY_2b_1b_pdf = getPdf("DY_2b_1b","for 2b contribution in 1b category",
                                                     mchist_2b_1b,mass,sf_2b);
        // Other Background contribution
        TH1D* mchist_ob_1b = dynamic_cast<TH1D*>(input_file->Get("2jets1btagR/mh/OS_Isolated/Central/other_bkg/m_tt_vis"));
        RooHistPdf DY_ob_1b_pdf = getPdf("DY_ob_1b","for other bkg contribution in 1b category",
                                                     mchist_ob_1b,mass,sf_ob);
        // Adding all contributions and create the total pdf
        RooAddPdf mcpdf_1b("mcpdf_1b","0b+1b+2b+other_bkg for 1 btag category",
                           RooArgList(DY_0b_1b_pdf,DY_1b_1b_pdf,DY_2b_1b_pdf,DY_ob_1b_pdf));

        //<<<<<<<<<<<<<<<< For 2b Category >>>>>>>>>>>>>>>>>>>
        // 0b contribution
        TH1D* mchist_0b_2b = dynamic_cast<TH1D*>(input_file->Get("2jets2btagR/mh/OS_Isolated/Central/DY_0b/m_tt_vis"));
        RooHistPdf DY_0b_2b_pdf = getPdf("DY_0b_2b","for 0b contribution in 2b category",
                                                     mchist_0b_2b,mass,sf_0b);
        // 1b contribution
        TH1D* mchist_1b_2b = dynamic_cast<TH1D*>(input_file->Get("2jets2btagR/mh/OS_Isolated/Central/DY_1b/m_tt_vis"));
        RooHistPdf DY_1b_2b_pdf = getPdf("DY_1b_2b","for 1b contribution in 2b category",
                                                     mchist_1b_2b,mass,sf_1b);
        // 2b contribution
        TH1D* mchist_2b_2b = dynamic_cast<TH1D*>(input_file->Get("2jets2btagR/mh/OS_Isolated/Central/DY_2b/m_tt_vis"));
        RooHistPdf DY_2b_2b_pdf = getPdf("DY_2b_2b","for 2b contribution in 2b category",
                                                     mchist_2b_2b,mass,sf_2b);
        // Other Background contribution
        TH1D* mchist_ob_2b = dynamic_cast<TH1D*>(input_file->Get("2jets2btagR/mh/OS_Isolated/Central/other_bkg/m_tt_vis"));
        RooHistPdf DY_ob_2b_pdf = getPdf("DY_ob_2b","for other bkg contribuion in 2b category",
                                                     mchist_ob_2b,mass,sf_ob);
        // Adding all contributions and create the total pdf
        RooAddPdf mcpdf_2b("mcpdf_2b","0b+1b+2b+other_bkg for 2 btag category",
                           RooArgList(DY_0b_2b_pdf,DY_1b_2b_pdf,DY_2b_2b_pdf,DY_ob_2b_pdf));
        std::cout<<"Monte carlo histograms are created"<<std::endl;

        //Define Category
        RooCategory categories("categories","categories") ;
        categories.defineType("0b_tag") ;
        categories.defineType("1b_tag") ;
        categories.defineType("2b_tag") ;
        std::cout<<"Categories are created"<<std::endl;

        // Construct combined dataset in (mass,category)
        RooDataHist combData("combData","combined data",mass,Index(categories),Import("0b_tag",data0b),
                             Import("1b_tag",data1b),Import("2b_tag",data2b));
        std::cout<<"combdata is created"<<std::endl;
        // Construct a simultaneous pdf in (mass,categories)

        // Construct a simultaneous pdf using category "categories" as index
        RooSimultaneous simPdf("simPdf","simultaneous pdf",categories) ;

        // Associate pdfs with the categories
        simPdf.addPdf(mcpdf_0b,"0b_tag") ;
        simPdf.addPdf(mcpdf_1b,"1b_tag") ;
        simPdf.addPdf(mcpdf_2b,"2b_tag") ;

        // Perform simultaneous fit
        //simPdf.fitTo(combData,Extended(kTRUE)) ;

        //Plotting
        RooPlot* frame0 = mass.frame() ;
        combData.plotOn(frame0,Cut("categories==categories::0b_tag")) ;
        simPdf.plotOn(frame0,Slice(categories,"0b_tag"),ProjWData(categories,combData)) ;
        simPdf.plotOn(frame0,Slice(categories,"0b_tag"),Components("DY_ob_0b_pdf"),ProjWData(categories,combData),LineStyle(kDashed)) ;


        RooPlot* frame1 = mass.frame() ;
        combData.plotOn(frame1,Cut("categories==categories::1b_tag")) ;
        simPdf.plotOn(frame1,Slice(categories,"1b_tag"),ProjWData(categories,combData)) ;
        simPdf.plotOn(frame1,Slice(categories,"1b_tag"),Components("DY_ob_1b_pdf"),ProjWData(categories,combData),LineStyle(kDashed)) ;

        RooPlot* frame2 = mass.frame() ;
        combData.plotOn(frame2,Cut("categories==categories::2b_tag")) ;
        simPdf.plotOn(frame2,Slice(categories,"2b_tag"),ProjWData(categories,combData)) ;
        simPdf.plotOn(frame2,Slice(categories,"2b_tag"),Components("DY_ob_2b_pdf"),ProjWData(categories,combData),LineStyle(kDashed)) ;

        TCanvas* c = new TCanvas("rf501_simultaneouspdf","rf403_simultaneouspdf",800,400) ;
        c->Divide(4) ;
        c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame0->GetYaxis()->SetTitleOffset(1.4) ; frame0->Draw() ;
        c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
        c->cd(3) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;

        output_file->cd();
        c->Write();*/

        TCanvas* c = new TCanvas("rf501_simultaneouspdf","rf403_simultaneouspdf",800,400) ;
        RooPlot* fram0 = mass.frame();
        DY_0b_0b_pdf.plotOn(fram0);
        fram0->Draw();
        output_file->cd();
        c->Write();

    }
private:
    Arguments args;
    std::shared_ptr<TFile> input_file, output_file;

    RooHistPdf getPdf(std::string name, std::string title, TH1D* h, RooAbsReal& mass,
                                RooAbsReal& scale_factor){
        RooRealVar norm("norm","norm",h->Integral());
        RooFormulaVar frac("frac","frac","@0*@1",RooArgList(scale_factor,norm));
        h->Scale(1./h->Integral());
        RooDataHist dataHist(name.c_str(),("Histogram for "+title).c_str(),mass,Import(*h)) ;
        RooHistPdf histpdf((name+"_pdf").c_str(),("Pdf for "+title).c_str(),mass,dataHist);
        //RooExtendPdf expdf((name+"_expdf").c_str(),("Extended pdf "+title).c_str(),histpdf,frac);
        return histpdf;
    }

};

} // namesapce analysis
PROGRAM_MAIN(analysis::DY_estimation, Arguments) // definition of the main program function

