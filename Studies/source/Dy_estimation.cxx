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
    OPT_ARG(analysis::Range<double>, fit_range, analysis::Range<double>(0, 500));
    OPT_ARG(std::string, var_name, "mass");
    OPT_ARG(std::string, histo_name, "m_tt_vis");

};

using namespace RooFit;
namespace analysis {
struct Contribution{ // list of Variables to create an extended pdf for a contribution
    std::shared_ptr<TH1D> histogram;
    std::shared_ptr<RooDataHist> rooHistogram;
    std::shared_ptr<RooHistPdf> pdf;
    std::shared_ptr<RooExtendPdf> expdf;
    std::shared_ptr<RooRealVar> norm;
    std::shared_ptr<RooFormulaVar> fraction;

    Contribution() :
    histogram(nullptr),
    rooHistogram(nullptr),
    pdf(nullptr),
    expdf(nullptr),
    norm(nullptr),
    fraction(nullptr)
    {

    }
};

struct CategoryModel{
    std::map<std::string,Contribution> mc_contributions,
    RooAddPdf* sum_pdf;
};

class Dy_estimation { // simple analyzer definition
public:
    Dy_estimation(const Arguments& _args) : args(_args),
    x(args.var_name().c_str(), args.var_name().c_str(), args.fit_range().min(), args.fit_range().max()),
    _histo_name(args.histo_name()),
    input_file(root_ext::OpenRootFile(args.input_file())),
    output_file(root_ext::CreateRootFile(args.output_file()))
    {
        // Analyzer initialization (e.g. open input/output files, parse configs...)
    }
    void Run()
    {
        //Data Histograms
        // For 0b category
        std::shared_ptr<TH1D> dataHisto0b(input_file->Get(("2jets0btagR/mh/OS_Isolated/Central/Data_SingleMuon/"
                                                                + _histo_name).c_str()));
        RooDataHist data0b("data0b","data for 0b category",x,Import(*dataHisto0b));
        // For 1b Category
        std::shared_ptr<TH1D> dataHisto1b(input_file->Get(("2jets1btagR/mh/OS_Isolated/Central/Data_SingleMuon/"
                                                                + _histo_name).c_str()));
        RooDataHist data1b("data1b","data for 1b category",x,Import(*dataHisto1b));
        // For 2b Category
        std::shared_ptr<TH1D> dataHisto2b(input_file->Get(("2jets2btagR/mh/OS_Isolated/Central/Data_SingleMuon/"
                                                                + _histo_name).c_str()));
        RooDataHist data2b("data1b","data for 2b category",x,Import(*dataHisto2b));

        //Create the Categories
        //For 0b Category
        std::shared_ptr<TH1D> mchist_0b_0b(input_file->Get(("2jets0btagR/mh/OS_Isolated/Central/DY_0b/"
                                                                 + _histo_name).c_str()));
        std::shared_ptr<TH1D> mchist_1b_0b(input_file->Get(("2jets0btagR/mh/OS_Isolated/Central/DY_1b/"
                                                                  + _histo_name).c_str()));
        std::shared_ptr<TH1D> mchist_2b_0b(input_file->Get(("2jets0btagR/mh/OS_Isolated/Central/DY_2b/"
                                                                  + _histo_name).c_str()));
        std::shared_ptr<TH1D> mchist_ob_0b(input_file->Get(("2jets0btagR/mh/OS_Isolated/Central/other_bkg/"
                                                                 + _histo_name).c_str()));
        createCategory("0b",mchist_0b_0b,mchist_1b_0b,mchist_2b_0b,mchist_ob_0b);

        //For 1b Category
        std::shared_ptr<TH1D> mchist_0b_1b(input_file->Get(("2jets1btagR/mh/OS_Isolated/Central/DY_0b/"
                                                                 + _histo_name).c_str()));
        std::shared_ptr<TH1D> mchist_1b_1b(input_file->Get(("2jets1btagR/mh/OS_Isolated/Central/DY_1b/"
                                                                 + _histo_name).c_str()));
        std::shared_ptr<TH1D> mchist_2b_1b(input_file->Get(("2jets1btagR/mh/OS_Isolated/Central/DY_2b/"
                                                                  + _histo_name).c_str()));
        std::shared_ptr<TH1D> mchist_ob_1b(input_file->Get(("2jets1btagR/mh/OS_Isolated/Central/other_bkg/"
                                                                  +_histo_name).c_str()));
        createCategory("1b",mchist_0b_1b,mchist_1b_1b,mchist_2b_1b,mchist_ob_1b);

        //For 2b Category
        std::shared_ptr<TH1D> mchist_0b_2b(input_file->Get(("2jets2btagR/mh/OS_Isolated/Central/DY_0b/"
                                                                  + _histo_name).c_str()));
        std::shared_ptr<TH1D> mchist_1b_2b(input_file->Get(("2jets2btagR/mh/OS_Isolated/Central/DY_1b/"
                                                                  + _histo_name).c_str()));
        std::shared_ptr<TH1D> mchist_2b_2b(input_file->Get(("2jets2btagR/mh/OS_Isolated/Central/DY_2b/"
                                                                  + _histo_name).c_str()));
        std::shared_ptr<TH1D> mchist_ob_2b(input_file->Get(("2jets2btagR/mh/OS_Isolated/Central/other_bkg/"
                                                                  + _histo_name).c_str()));
        createCategory("2b",mchist_0b_2b,mchist_1b_2b,mchist_2b_2b,mchist_ob_2b);

        //Define Category
        RooCategory categories("categories","categories") ;
        categories.defineType("0b_tag") ;
        categories.defineType("1b_tag") ;
        categories.defineType("2b_tag") ;
        std::cout<<"Categories are created"<<std::endl;

        // Construct combined dataset in (mass,category)
        RooDataHist combData("combData","combined data",x,Index(categories),Import("0b_tag",data0b),
                             Import("1b_tag",data1b),Import("2b_tag",data2b));
        std::cout<<"combdata is created"<<std::endl;
        // Construct a simultaneous pdf in (mass,categories)

        // Construct a simultaneous pdf using category "categories" as index
        RooSimultaneous simPdf("simPdf","simultaneous pdf",categories) ;

        // Associate pdfs with the categories
        CategoryModel cat_0b = Category_Map["0b"];
        CategoryModel cat_1b = Category_Map["1b"];
        CategoryModel cat_2b = Category_Map["2b"];

        simPdf.addPdf(*(cat_0b.sum_pdf),"0b_tag") ;
        simPdf.addPdf(*(cat_1b.sum_pdf),"1b_tag") ;
        simPdf.addPdf(*(cat_2b.sum_pdf),"2b_tag") ;

        // Perform simultaneous fit
        simPdf.fitTo(combData,Extended(kTRUE)) ;

        //Plotting
        RooPlot* frame0 = x.frame() ;
        combData.plotOn(frame0,Cut("categories==categories::0b_tag")) ;
        simPdf.plotOn(frame0,Slice(categories,"0b_tag"),ProjWData(categories,combData)) ;
        simPdf.plotOn(frame0,Slice(categories,"0b_tag"),Components("DY_ob_0b_pdf"),ProjWData(categories,combData),LineStyle(kDashed)) ;


        RooPlot* frame1 = x.frame() ;
        combData.plotOn(frame1,Cut("categories==categories::1b_tag")) ;
        simPdf.plotOn(frame1,Slice(categories,"1b_tag"),ProjWData(categories,combData)) ;
        simPdf.plotOn(frame1,Slice(categories,"1b_tag"),Components("DY_ob_1b_pdf"),ProjWData(categories,combData),LineStyle(kDashed)) ;

        RooPlot* frame2 = x.frame() ;
        combData.plotOn(frame2,Cut("categories==categories::2b_tag")) ;
        simPdf.plotOn(frame2,Slice(categories,"2b_tag"),ProjWData(categories,combData)) ;
        simPdf.plotOn(frame2,Slice(categories,"2b_tag"),Components("DY_ob_2b_pdf"),ProjWData(categories,combData),LineStyle(kDashed)) ;

        TCanvas* c = new TCanvas("rf501_simultaneouspdf","rf403_simultaneouspdf",800,400) ;
        c->Divide(4) ;
        c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame0->GetYaxis()->SetTitleOffset(1.4) ; frame0->Draw() ;
        c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
        c->cd(3) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;

        output_file->cd();
        c->Write();
  }

private:
    Arguments args;
    std::shared_ptr<TFile> input_file, output_file;
    std::string _histo_name;
    enum scale_factors {_0b, _1b, _2b, _ob};
    //X axis
    RooRealVar x;

    //Scale Factors
    RooRealVar sf_0b{"sf_0b","Scale factor for 0b contibution",1,0.1,10};
    RooRealVar sf_1b{"sf_1b","Scale factor for 1b contibution",1,0.1,10};
    RooRealVar sf_2b{"sf_2b","Scale factor for 2b contibution",1,0.1,10};
    RooRealVar sf_ob{"sf_ob","Scale factor for other bkg contibution",1,0.1,10};

    std::map<std::string,CategoryModel> Category_Map;

    Contribution setPdf(std::string name, std::string title,std::shared_ptr<TH1D> h, scale_factors sf){
        Contribution temp ;
        RooRealVar* scale_factor;
        switch(sf) {
            case _0b :
                scale_factor = &sf_0b;
                break;
            case  _1b:
                scale_factor = &sf_1b;
                break;
            case _2b:
                scale_factor = &sf_2b;
                break;
            case _ob:
                scale_factor = &sf_ob;
                break;
        }
        temp.histogram = h;
        std::shared_ptr<RooRealVar> norm(("norm_"+name).c_str(),("norm for "+title).c_str(),h->Integral());
        temp.norm = norm;
        std::shared_ptr<RooFormulaVar> frac(("frac_"+name).c_str(),("fraction of "+title).c_str(),
                                                "@0*@1",RooArgList(*scale_factor,*norm));
        temp.fraction = frac;
        h->Scale(1./h->Integral());
        std::shared_ptr<RooDataHist> dataHist(("dataHist_"+name).c_str(),("RooHistogram for "+title).c_str(),
                                               x,Import(*h)) ;
        temp.rooHistogram = dataHist;
        std::shared_ptr<RooHistPdf> histpdf(("DY_"+name+"_pdf").c_str(),("Pdf for "+title).c_str(),
                                             x,*dataHist);
        temp.pdf = histpdf;
        std::shared_ptr<RooExtendPdf> expdf(("DY_"+name+"_expdf").c_str(),("Extended pdf "+title).c_str(),
                                               *histpdf,*frac);
        temp.expdf = expdf;

        return temp;

    }

    void createCategory(std::string cat_name, std::shared_ptr<TH1D> histo_0b, std::shared_ptr<TH1D> histo_1b,
                        std::shared_ptr<TH1D> histo_2b, std::shared_ptr<TH1D> histo_ob){
        // 0b contribution
        scale_factors sf = _0b;
        Contribution DY_0b = setPdf("0b_"+cat_name,"for 0b contribution in "+ cat_name +" category",histo_0b,sf);

        // 1b contribution
        sf = _1b;
        Contribution DY_1b = setPdf("1b_"+cat_name,"for 1b contribution in "+ cat_name +" category",histo_1b,sf);

        // 2b contribution
        sf = _2b;
        Contribution DY_2b = setPdf("2b_"+cat_name,"for 2b contribution in "+ cat_name + " category",histo_2b,sf);

        // Other Background contribution
        sf = _ob;
        Contribution DY_ob = setPdf("ob_"+cat_name,"for other bbkg contribution in " + cat_name + " category",
                                    histo_ob,sf);

        // Adding all contributions and create the total pdf
        RooAddPdf* sumpdf = new RooAddPdf(("sumpdf_"+cat_name).c_str(),
                                          ("0b+1b+2b+other_bkg for " + cat_name + " category").c_str(),
                                            RooArgList(*(DY_0b.expdf),*(DY_1b.expdf),*(DY_2b.expdf),*(DY_ob.expdf)));

        CategoryModel category;
        category.cont_0b = DY_0b;
        category.cont_1b = DY_1b;
        category.cont_2b = DY_2b;
        category.cont_ob = DY_ob;
        category.sum_pdf = sumpdf;

        Category_Map[cat_name] = category;

    }

};

} // namesapce analysis
PROGRAM_MAIN(analysis::Dy_estimation, Arguments) // definition of the main program function

