/*! Final analysis step to estimate scale-factors for DY normalization using the muMu channel
 in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "Analysis/include/EventAnalyzerDataId.h"
#include "Analysis/include/AnalysisCategories.h"

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
#include "TLatex.h"
#include "RooPlot.h"
#include "RooFitResult.h"

struct Arguments {
    REQ_ARG(std::string, input_file);
    REQ_ARG(std::string, output_file);
    OPT_ARG(analysis::Range<double>, fit_range, analysis::Range<double>(0, 500));
    OPT_ARG(std::string, var_name, "mass");
    OPT_ARG(std::string, histo_name, "m_tt_vis");
    OPT_ARG(analysis::Range<double>, scale_factor_range, analysis::Range<double>(0.1, 2));
    OPT_ARG(analysis::DYFitModel, fit_method,analysis::DYFitModel::NbjetBins);
    OPT_ARG(int, med_ht, 80);
    OPT_ARG(int, high_ht, 150);
    OPT_ARG(int, low_nJet, 0);
    OPT_ARG(int, high_nJet, 4);
    OPT_ARG(int, vlowPt, 0);
    OPT_ARG(int, lowPt, 20);
    OPT_ARG(int, medPt, 40);
    OPT_ARG(int, medPt1, 30);
    OPT_ARG(int, medPt2, 50);
    OPT_ARG(int, highPt, 100);
    OPT_ARG(int, vhighPt,200);
    OPT_ARG(std::string, sub_category_option, "mh");
    OPT_ARG(std::string, sample_order, "LO");
};


namespace analysis {
using namespace RooFit;

static const std::string& GetCategoryName(const EventCategory& category){

    static std::map<EventCategory,std::string> eventCategories{ {  EventCategory::Parse("2j0b"), "2j0b"} ,
        {EventCategory::Parse("2j1b"), "2j1b"}, { EventCategory::Parse("2j2b+"), "2j2b" }};

    return eventCategories.at(category);
}

struct Contribution{
    std::string name;
    std::shared_ptr<TH1D> histogram;
    RooRealVar norm;
    RooFormulaVar fraction;
    RooDataHist rooHistogram;
    RooHistPdf pdf;
    RooExtendPdf expdf;

    Contribution(std::shared_ptr<TFile> input_file, const EventAnalyzerDataId& dataId, const std::string& hist_name,
                 const RooRealVar& x, const RooRealVar& scale_factor) :
    name(boost::str(boost::format("%1%_%2%_%3%") % GetCategoryName(dataId.Get<EventCategory>()) % dataId.Get<EventSubCategory>()
                    % dataId.Get<std::string>())),
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
                  const std::vector<std::string>& contribution_names,const std::string& hist_name, const RooRealVar& x,
                  const std::map<std::string,std::shared_ptr<RooRealVar>>& scale_factor_map) :
    name(GetCategoryName(catId.Get<EventCategory>())+"_"+ToString(catId.Get<EventSubCategory>()))
    {
        RooArgList pdf_list;
        for(const std::string& contrib_name : contribution_names){
            EventAnalyzerDataId dataId = catId.Set(contrib_name);
             mc_contributions[contrib_name] = std::make_shared<Contribution>(input_file,dataId,hist_name,x,
                                                                             *(scale_factor_map.at(contrib_name)));
             if(mc_contributions[contrib_name]->histogram->Integral() > 0)
                pdf_list.add(mc_contributions[contrib_name]->expdf);
             mc_contributions[contrib_name]->expdf.Print();
        }
        sum_pdf = std::make_shared<RooAddPdf>(("sumpdf_"+name).c_str(),("Total Pdf for "+name).c_str(),pdf_list);
    }
};
class DY_norm_baseline {
public:
    DY_norm_baseline(const Arguments& _args) :
    args(_args), input_file(root_ext::OpenRootFile(args.input_file())),
    output_file(root_ext::CreateRootFile(args.output_file())),
    x(args.var_name().c_str(), args.var_name().c_str(), args.fit_range().min(), args.fit_range().max()),
    fit_model(args.fit_method())
    {
    }

    void Run()
    {
        static std::string dy_contrib_prefix;
        if(args.sample_order() == "LO") dy_contrib_prefix= "DY_lo";
        else if(args.sample_order() == "NLO") dy_contrib_prefix = "DY_nlo";
        EventSubCategory base_sub_category;
        if(args.sub_category_option() == "mh") base_sub_category =
                EventSubCategory().SetCutResult(SelectionCut::mh, true).SetCutResult(SelectionCut::lowMET,true);
        else if(args.sub_category_option() == "mtt") base_sub_category =
                EventSubCategory().SetCutResult(SelectionCut::mtt, true).SetCutResult(SelectionCut::lowMET,true);
        const std::vector<int> ht_points = { 0, args.med_ht(), args.high_ht() };
        const std::vector<int> nJet_points = {args.low_nJet(), args.high_nJet() };
        std::vector<int> pt_points;
        if(args.sample_order() == "NLO") pt_points = {args.vlowPt(), args.lowPt(), args.medPt(), args.highPt()};
        else if(args.sample_order() == "LO") pt_points = {args.vlowPt(), args.lowPt(), args.medPt1(), args.medPt2(), args.highPt(),
                                                          args.vhighPt()};
        static const size_t max_n_b = 2;
        if(fit_model == DYFitModel::NbjetBins) {
            subCategories = { base_sub_category };
            for(size_t nb = 0; nb <= max_n_b; ++nb) {
                //for(int nJet : nJet_points){
                    const std::string name = boost::str(boost::format("%1%_%2%b") % dy_contrib_prefix % nb);
                    contribution_names.push_back(name);
                //}
            }

        } else {
            /*subCategories = { EventSubCategory(base_sub_category).SetCutResult(SelectionCut::lowHT, true),
                              EventSubCategory(base_sub_category).SetCutResult(SelectionCut::medHT, true),
                              EventSubCategory(base_sub_category).SetCutResult(SelectionCut::highHT, true)};*/
            subCategories = { base_sub_category };
            for(size_t nb = 0; nb <= max_n_b; ++nb) {
                if(fit_model==DYFitModel::NbjetBins_htBins){
                    for(int ht : ht_points) {
                        const std::string name = boost::str(boost::format("%1%_%2%b_%3%ht") % dy_contrib_prefix % nb % ht);
                        contribution_names.push_back(name);
                    }
                }
                else if(fit_model==DYFitModel::NbjetBins_NjetBins){
                    for(int nJet : nJet_points) {
                        const std::string name = boost::str(boost::format("%1%_%2%b_%3%Jet") % dy_contrib_prefix % nb % nJet);
                        contribution_names.push_back(name);
                    }
                }
                else if(fit_model==DYFitModel::NbjetBins_ptBins){
                    if(args.sample_order() == "NLO")
                        subCategories = { EventSubCategory(base_sub_category).SetCutResult(SelectionCut::vlowPtNLO, true),
                                        EventSubCategory(base_sub_category).SetCutResult(SelectionCut::lowPtNLO, true),
                                        EventSubCategory(base_sub_category).SetCutResult(SelectionCut::medPtNLO, true),
                                        EventSubCategory(base_sub_category).SetCutResult(SelectionCut::highPtNLO, true)};
                    else if(args.sample_order() == "LO")
                        subCategories = { EventSubCategory(base_sub_category).SetCutResult(SelectionCut::vlowPtLO, true),
                                        EventSubCategory(base_sub_category).SetCutResult(SelectionCut::lowPtLO, true),
                                        EventSubCategory(base_sub_category).SetCutResult(SelectionCut::medPt1LO, true),
                                        EventSubCategory(base_sub_category).SetCutResult(SelectionCut::medPt2LO, true),
                                        EventSubCategory(base_sub_category).SetCutResult(SelectionCut::highPtLO, true),
                                        EventSubCategory(base_sub_category).SetCutResult(SelectionCut::vhighPtLO, true)};
                    for(int pt : pt_points){
                        const std::string name = boost::str(boost::format("%1%_%2%b_%3%Pt") % dy_contrib_prefix % nb % pt);
                        contribution_names.push_back(name);
                    }
                }
            }
        }
        contribution_names.push_back("other_bkg_muMu");
        contribution_names.push_back("QCD");
        /*contribution_names.push_back("WW");
        contribution_names.push_back("WZ");
        contribution_names.push_back("Wjets");
        contribution_names.push_back("tW");
        contribution_names.push_back("EWK");
        contribution_names.push_back("ZH");
        contribution_names.push_back("ZZ");
        contribution_names.push_back("TT");*/

        /*std::map<std::string,std::pair<double,double>>scale_factor_values =
                                                          {{"DY_MC_0b_0Jet",{1.09725,2.37998e-03}} ,
                                                           {"DY_MC_0b_4Jet",{0.806768,8.59454e-03}} ,
                                                           {"DY_MC_1b_0Jet",{1.15362,2.37680e-02}} ,
                                                           {"DY_MC_1b_4Jet",{1.06926,7.12087e-02}} ,
                                                           {"DY_MC_2b_0Jet",{1.2323,5.77241e-02}} ,
                                                           {"DY_MC_2b_4Jet",{0.962531,6.87868e-02}} ,
                                                           {"WW",{1.0, 0 }},
                                                           {"WZ",{1.0, 0 }},
                                                           {"Wjets",{1.0, 0 }},
                                                           {"tW",{1.0, 0 }},
                                                           {"EWK",{1.0, 0 }},
                                                           {"ZH",{1.0, 0 }},
                                                           {"ZZ",{1.0, 0 }},
                                                           {"TT",{1.0, 0 }},};*/


        std::map<std::string, std::shared_ptr<RooRealVar>> scale_factor_map;
        for (const std::string& contrib_name: contribution_names){
            //if(contrib_name.find("DY") != std::string::npos){
            //    double value = scale_factor_values[contrib_name].first;
            //    double error = scale_factor_values[contrib_name].second;
                scale_factor_map[contrib_name] = std::make_shared<RooRealVar>
                    (("sf_"+contrib_name).c_str(),("Scale Factor for contribution "+contrib_name).c_str(),
                    1,args.scale_factor_range().min(),args.scale_factor_range().max());
            //}
           /*else{
                scale_factor_map[contrib_name] = std::make_shared<RooRealVar>
                        (("sf_"+contrib_name).c_str(),("Scale Factor for contribution "+contrib_name).c_str(),
                         1);
           }*/
        }
        std::string data_folder = "Data_SingleMuon";
        std::map<std::string,TH1*> dataCategories;
        std::map<std::string,std::shared_ptr<CategoryModel>> categories;

        RooCategory rooCategories("rooCategories","rooCategories") ;
        RooSimultaneous simPdf("simPdf","simultaneous pdf",rooCategories) ;

        for(const EventCategory& cat : eventCategories_vec){
            for (const EventSubCategory& subCategory : subCategories){
                EventAnalyzerDataId catId{cat, subCategory, EventRegion::SignalRegion(), UncertaintySource::None,
                                          UncertaintyScale::Central};
                std::string category = GetCategoryName(cat);
                std::string subcategory= ToString(catId.Get<EventSubCategory>());
                categories[category+subcategory] = std::make_shared<CategoryModel>(input_file,catId,contribution_names,
                                                                                   args.histo_name(),x,scale_factor_map);
                EventAnalyzerDataId dataId = catId.Set(data_folder);
                dataCategories[category+subcategory] = root_ext::ReadObject<TH1>(*input_file,dataId.GetName()+ "/" +
                                                                                 args.histo_name());
                rooCategories.defineType((category+subcategory).c_str()) ;
                categories[category+subcategory]->sum_pdf->Print();
                simPdf.addPdf(*(categories[category+subcategory]->sum_pdf),(category+subcategory).c_str());
            }
        }

        // Construct combined dataset in (x,rooCategories)
        RooDataHist combData("combData","combined data",x,rooCategories,dataCategories);

        // Perform simultaneous fit
        RooFitResult* result = simPdf.fitTo(combData,Extended(kTRUE),Save(),Verbose(false));//,PrintLevel(3)) ;

        //Saving Results
        output_file->mkdir(ToString(fit_model).c_str());
        output_file->cd(ToString(fit_model).c_str());

        const TMatrixDSym& correaltion_matrix = result->correlationMatrix();
        const TMatrixDSym& covariance_matrix = result->covarianceMatrix();
        int nRows = covariance_matrix.GetNrows();
        int nColumns = covariance_matrix.GetNcols();
        auto cov_hist = std::make_shared<TH2D>("covariance_matrix","covariance matrix",
                                                                        nRows,0.5,0.5+nRows,nColumns,0.5,0.5+nColumns);
        auto cor_hist = std::make_shared<TH2D>("correlation_matrix","correlation matrix",
                                                                        nRows,0.5,0.5+nRows,nColumns,0.5,0.5+nColumns);
        for(int i = 0; i<nRows;i++){
            for(int j = 0; j<nColumns;j++){
                double cov = covariance_matrix[i][j];
                cov_hist->SetBinContent(i+1,j+1,cov);
                double cor = correaltion_matrix[i][j];
                cor_hist->SetBinContent(i+1,j+1,cor);
            }
        }

        auto scale_factors_hist = std::make_shared<TH1D>("scale_factors","Scale factors afte the fit",
                                                                     nRows,0.5,0.5+nRows);
        int i=1;
        for (const std::string& contrib_name: contribution_names){
            //if(contrib_name != "other_bkg_muMu"){
                cov_hist->GetXaxis()->SetBinLabel(i,contrib_name.c_str());
                cov_hist->GetYaxis()->SetBinLabel(i,contrib_name.c_str());
                cor_hist->GetXaxis()->SetBinLabel(i,contrib_name.c_str());
                cor_hist->GetYaxis()->SetBinLabel(i,contrib_name.c_str());
            //}

            scale_factors_hist->GetXaxis()->SetBinLabel(i,contrib_name.c_str());
            scale_factors_hist->SetBinContent(i,scale_factor_map[contrib_name]->getValV());
            scale_factors_hist->SetBinError(i,scale_factor_map[contrib_name]->getError());
            i++;
            //}
        }
        cov_hist->Write();
        cor_hist->Write();
        scale_factors_hist->Write();

        //Plotting
        for(auto& cat : eventCategories_vec){
            for(const EventSubCategory& sub_cat: subCategories){
                TCanvas* c = new TCanvas(("fit_"+GetCategoryName(cat)+"_"+ToString(sub_cat)).c_str(),
                                     ("fit in eventCategory " + GetCategoryName(cat) + ToString(sub_cat)).c_str(),800,400) ;
                RooPlot* frame = x.frame() ;
                //simPdf.getPdf((it->second+ToString(sub_cat)).c_str())->Print();
                //std::cout<<" RooCategory cat = "<<ToString(cat)<<std::endl;
                combData.plotOn(frame,Cut((static_cast<std::string>("rooCategories==rooCategories::")+
                                    GetCategoryName(cat)+ToString(sub_cat)).c_str())) ;
                simPdf.plotOn(frame,Slice(rooCategories,(GetCategoryName(cat)+ToString(sub_cat)).c_str()),
                                  ProjWData(rooCategories,combData)) ;
                simPdf.plotOn(frame,Slice(rooCategories,(GetCategoryName(cat)+ToString(sub_cat)).c_str()),
                              Components((static_cast<std::string>("expdf_")+GetCategoryName(cat)+
                                                   static_cast<std::string>("_")+ToString(sub_cat) +
                                                   static_cast<std::string>("_other_bkg_muMu")).c_str()),
                                ProjWData(rooCategories,combData),LineStyle(kDashed)) ;
                /*simPdf.plotOn(frame,Slice(rooCategories,(GetCategoryName(cat)+ToString(sub_cat)).c_str()),
                              Components((static_cast<std::string>("expdf_")+GetCategoryName(cat)+
                                                   static_cast<std::string>("_")+ToString(sub_cat) +
                                                   static_cast<std::string>("_QCD")).c_str()),
                                ProjWData(rooCategories,combData),LineColor(kBlack),LineStyle(kDashed)) ;*/

                gPad->SetLeftMargin(0.15f) ; frame->GetYaxis()->SetTitleOffset(1.4f) ; frame->Draw() ;

                TLatex *tex = new TLatex(0.87,   0.83,"CMS");
                tex->SetNDC(true);
                tex->SetTextFont(61);
                tex->SetTextAlign(31);
                tex->SetTextSize(0.06f);
                tex->Draw();

                tex = new TLatex(0.87,   0.76,GetCategoryName(cat).c_str());
                tex->SetNDC(true);
                tex->SetTextFont(52);
                tex->SetTextAlign(31);
                tex->SetTextSize(0.0456f);
                tex->Draw();

                tex = new TLatex(0.87,   0.78,"Unpublished");
                tex->SetNDC(true);
                tex->SetTextFont(52);
                tex->SetTextAlign(31);
                tex->SetTextSize(0.04f);
                tex->Draw();

                tex = new TLatex(0.9,   0.91,"28.6 fb^{-1} (13 TeV)");
                tex->SetNDC(true);
                tex->SetTextFont(42);
                tex->SetTextAlign(31);
                tex->SetTextSize(0.048f);
                tex->Draw();

                c->Write();
            }
        }
        /*TCanvas* c = new TCanvas("fit_","fit in eventCategory ",800,400) ;
        RooPlot* frame = x.frame();
        simPdf.getPdf("2j2bmh_lowMET")->Print();
        simPdf.getPdf("2j2bmh_lowMET")->plotOn(frame,Components("expdf_2j2b+_mh_lowMET_other_bkg_muMu"));
        simPdf.getPdf("2j1bmh_lowMET")->Print();
        simPdf.getPdf("2j1bmh_lowMET")->plotOn(frame,Components("expdf_2j1b_mh_lowMET_other_bkg_muMu"),LineColor(kBlack));
        frame->Draw();
        c->Write();
        output_file->Close();*/
  }

private:
    Arguments args;
    std::shared_ptr<TFile> input_file, output_file;
    //X axis
    RooRealVar x;


    std::vector<EventCategory> eventCategories_vec{ EventCategory::Parse("2j0b"), EventCategory::Parse("2j1b"), EventCategory::Parse("2j2b+") };

    std::vector<EventSubCategory> subCategories;

    std::vector<std::string> contribution_names;
    DYFitModel fit_model;
};

} // namesapce analysis
PROGRAM_MAIN(analysis::DY_norm_baseline, Arguments) // definition of the main program function
