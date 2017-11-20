
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
        std::map<std::string,TH1*> dataCategories;
        std::map<std::string,std::shared_ptr<CategoryModel>> categories;
        for(const EventCategory& cat : eventCategories ){
            EventAnalyzerDataId catId = metaId.Set(cat);
            std::string category = ToString(catId.Get<EventCategory>());
            categories[category] = std::make_shared<CategoryModel>(input_file,catId,contribution_names,args.histo_name()
                                                                   ,x,scale_factor_map);
            EventAnalyzerDataId dataId = catId.Set(data_folder);
            dataCategories[category] = root_ext::ReadObject<TH1>(*input_file,dataId.GetName()+ "/" + args.histo_name());
        }

        //Define a RooCategory
        RooCategory rooCategories("rooCategories","rooCategories") ;
        for(const EventCategory& cat : eventCategories ){
            rooCategories.defineType(ToString(cat).c_str()) ;
        }

        // Construct combined dataset in (x,category)
        RooDataHist combData("combData","combined data",x,rooCategories,dataCategories);

        // Construct a simultaneous pdf in (mass,categories)
        RooSimultaneous simPdf("simPdf","simultaneous pdf",rooCategories) ;
        for(const EventCategory& cat : eventCategories ){
            std::string category = ToString(cat);
            simPdf.addPdf(*(categories[category]->sum_pdf),category.c_str());
        }

        // Perform simultaneous fit
        simPdf.fitTo(combData,Extended(kTRUE)) ;

        //Plotting
        TCanvas* c = new TCanvas("rf501_simultaneouspdf","rf403_simultaneouspdf",800,400) ;
        c->Divide(4) ;
        size_t i=0;
        for(const EventCategory& cat : eventCategories ){
            RooPlot* frame = x.frame() ;
            combData.plotOn(frame,Cut(((std::string)("rooCategories==rooCategories::")+
                                    ToString(cat)).c_str())) ;
            simPdf.plotOn(frame,Slice(rooCategories,ToString(cat).c_str()),ProjWData(rooCategories,combData)) ;
            simPdf.plotOn(frame,Slice(rooCategories,ToString(cat).c_str()),
                      Components(((std::string)("expdf_")+ToString(cat)+(std::string)("_other_bkg")).c_str() ),
                      ProjWData(rooCategories,combData),LineStyle(kDashed)) ;
             c->cd(i+1) ; gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.4) ; frame->Draw() ;
             i++;
        }
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
};

} // namesapce analysis
PROGRAM_MAIN(analysis::Dy_estimation, Arguments) // definition of the main program function

