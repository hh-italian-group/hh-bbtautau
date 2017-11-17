
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "Analysis/include/EventAnalyzerDataId.h"
#include "Analysis/include/AnalysisCategories.h"

//Root and Roofit Headers
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

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_file); // required argument "input_file"
    REQ_ARG(std::string, output_file); // required argument "output_file"
    OPT_ARG(analysis::Range<double>, fit_range, analysis::Range<double>(0, 500));
    OPT_ARG(std::string, var_name, "mass");
    OPT_ARG(std::string, histo_name, "m_tt_vis");
    OPT_ARG(analysis::Range<double>, scale_factor_range, analysis::Range<double>(0.1, 2));
};

using namespace RooFit;
namespace analysis {
struct Contribution{ // list of Variables to create an extended pdf for a contribution
    std::string name;
    std::shared_ptr<TH1D> histogram;
    RooRealVar norm;
    RooFormulaVar fraction;
    RooDataHist rooHistogram;
    RooHistPdf pdf;
    RooExtendPdf expdf;

    Contribution(std::shared_ptr<TFile> input_file, const EventAnalyzerDataId& dataId, const std::string& hist_name,
                 const RooRealVar& x, const RooRealVar& scale_factor) :
    name(boost::str(boost::format("%1%_%2%") % dataId.Get<EventCategory>() % dataId.Get<std::string>())),
    histogram(root_ext::ReadObject<TH1D>(*input_file, dataId.GetName()+ "/" + hist_name)),
    norm(("norm_"+name).c_str(),("norm for "+name).c_str(),histogram->Integral()),
    fraction(("frac_"+name).c_str(),("fraction of "+name).c_str(),"@0*@1",RooArgList(scale_factor,norm)),
    rooHistogram(("rooHistogram_"+name).c_str(),("RooHistogram for "+name).c_str(),x,Import(*histogram)),
    pdf(("pdf_"+name).c_str(),("Pdf for "+name).c_str(),x,rooHistogram),
    expdf(("expdf_"+name).c_str(),("Extended pdf for "+name).c_str(),pdf,fraction)
    {
    }
};

struct CategoryModel{
    std::map<std::string,std::shared_ptr<Contribution>> mc_contributions;
    std::string name;
    std::shared_ptr<RooAddPdf> sum_pdf;
    CategoryModel(std::shared_ptr<TFile> input_file, const EventAnalyzerDataId& catId,
                  const std::set<std::string>& contribution_names,const std::string& hist_name, const RooRealVar& x,
                  const std::map<std::string,std::shared_ptr<RooRealVar>>& scale_factor_map) :
    name(ToString(catId.Get<EventCategory>()))
    {
        RooArgList pdf_list;
        for(const std::string& contrib_name : contribution_names){
            EventAnalyzerDataId dataId = catId.Set(contrib_name);
             mc_contributions[contrib_name] = std::make_shared<Contribution>(input_file,dataId,hist_name,x,
                                                                             *(scale_factor_map.at(contrib_name)));
             pdf_list.add(mc_contributions[contrib_name]->expdf);
        }
        sum_pdf = std::make_shared<RooAddPdf>(("sumpdf_"+name).c_str(),("Total Pdf for "+name).c_str(),pdf_list);
    }
};
class Dy_estimation { // simple analyzer definition
public:
    Dy_estimation(const Arguments& _args) : args(_args),
    x(args.var_name().c_str(), args.var_name().c_str(), args.fit_range().min(), args.fit_range().max()),
    input_file(root_ext::OpenRootFile(args.input_file())),
    output_file(root_ext::CreateRootFile(args.output_file()))
    {
        // Analyzer initialization (e.g. open input/output files, parse configs...)
    }
    void Run()
    {
        std::map<std::string, std::shared_ptr<RooRealVar>> scale_factor_map;
        for (const std::string& contrib_name: contribution_names){
            scale_factor_map[contrib_name] = std::make_shared<RooRealVar>
                    (("sf_"+contrib_name).c_str(),("Scale Factor for contribution "+contrib_name).c_str(),1,
                     args.scale_factor_range().min(),args.scale_factor_range().max());
        }
        std::string data_folder = "Data_SingleMuon";
        //std::map<std::string,std::shared_ptr<RooDataHist>> dataCategories;
        std::map<std::string,TH1*> dataCategories;
        std::map<std::string,std::shared_ptr<CategoryModel>> categories;
        for(const EventCategory& cat : eventCategories ){
            EventAnalyzerDataId catId = metaId.Set(cat);
            std::string category = ToString(catId.Get<EventCategory>());
            categories[category] = std::make_shared<CategoryModel>(input_file,catId,contribution_names,args.histo_name()
                                                                   ,x,scale_factor_map);
            EventAnalyzerDataId dataId = catId.Set(data_folder);
            //std::shared_ptr<TH1D> dataHisto = std::make_shared<TH1D>(*(root_ext::ReadObject<TH1D>(*input_file,
            //                                                            dataId.GetName()+ "/" + args.histo_name())));
            //dataCategories[category] = std::make_shared<RooDataHist>(("data_"+category).c_str(),
            //                                    ("RooHistogram for data for "+category).c_str(),x,Import(*dataHisto));
            dataCategories[category] = root_ext::ReadObject<TH1>(*input_file,dataId.GetName()+ "/" + args.histo_name());
        }

        RooCategory rooCategories("rooCategories","rooCategories") ;
        for(const EventCategory& cat : eventCategories ){
            rooCategories.defineType(ToString(cat).c_str()) ;
        }

        // Construct combined dataset in (mass,category)
        RooDataHist combData("combData","combined data",x,rooCategories,dataCategories);
        //RooDataHist combData("combData","combined data",x,Index(categories),Import(,data0b),
        //                     Import("1b_tag",data1b),Import("2b_tag",data2b));
/*        std::cout<<"combdata is created"<<std::endl;
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
*/
        //Plotting
        RooPlot* frame0 = x.frame() ;
        combData.plotOn(frame0,Cut(((std::string)("rooCategories==rooCategories::")+
                                    ToString(EventCategory::TwoJets_ZeroBtag_Resolved())).c_str())) ;
        //simPdf.plotOn(frame0,Slice(categories,"0b_tag"),ProjWData(categories,combData)) ;
        //simPdf.plotOn(frame0,Slice(categories,"0b_tag"),Components("DY_ob_0b_pdf"),ProjWData(categories,combData),LineStyle(kDashed)) ;

/*
        RooPlot* frame1 = x.frame() ;
        combData.plotOn(frame1,Cut("categories==categories::1b_tag")) ;
        simPdf.plotOn(frame1,Slice(categories,"1b_tag"),ProjWData(categories,combData)) ;
        simPdf.plotOn(frame1,Slice(categories,"1b_tag"),Components("DY_ob_1b_pdf"),ProjWData(categories,combData),LineStyle(kDashed)) ;

        RooPlot* frame2 = x.frame() ;
        combData.plotOn(frame2,Cut("categories==categories::2b_tag")) ;
        simPdf.plotOn(frame2,Slice(categories,"2b_tag"),ProjWData(categories,combData)) ;
        simPdf.plotOn(frame2,Slice(categories,"2b_tag"),Components("DY_ob_2b_pdf"),ProjWData(categories,combData),LineStyle(kDashed)) ;
*/
        TCanvas* c = new TCanvas("rf501_simultaneouspdf","rf403_simultaneouspdf",800,400) ;
        c->Divide(4) ;
        c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame0->GetYaxis()->SetTitleOffset(1.4) ; frame0->Draw() ;
        //c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
        //c->cd(3) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;

        output_file->cd();
        c->Write();
  }

private:
    Arguments args;
    std::shared_ptr<TFile> input_file, output_file;
    enum scale_factors {_0b, _1b, _2b, _ob};
    //X axis
    RooRealVar x;

    EventSubCategory subCategory = EventSubCategory().SetCutResult(SelectionCut::mh, true);
    EventAnalyzerDataId metaId{subCategory,EventRegion::SignalRegion(),EventEnergyScale::Central};

    EventCategorySet eventCategories{EventCategory::TwoJets_ZeroBtag_Resolved(),
                EventCategory::TwoJets_OneBtag_Resolved(),EventCategory::TwoJets_TwoBtag_Resolved()};

    std::set<std::string> contribution_names{"DY_0b","DY_1b","DY_2b","other_bkg"};


/*
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
*/
};

} // namesapce analysis
PROGRAM_MAIN(analysis::Dy_estimation, Arguments) // definition of the main program function

