/*! Final analysis step to estimate scale-factors for DY normalization using the muMu channel
 in the HH->bbtautau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
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
};
namespace analysis {
using namespace RooFit;

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
    name(boost::str(boost::format("%1%_%2%_%3%") % dataId.Get<EventCategory>() % dataId.Get<EventSubCategory>()
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
    name(ToString(catId.Get<EventCategory>())+"_"+ToString(catId.Get<EventSubCategory>()))
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
class DY_estimation {
public:
    DY_estimation(const Arguments& _args) :
    args(_args), input_file(root_ext::OpenRootFile(args.input_file())),
    output_file(root_ext::CreateRootFile(args.output_file())),
    x(args.var_name().c_str(), args.var_name().c_str(), args.fit_range().min(), args.fit_range().max()),
    fit_model(args.fit_method())
    {
    }

    void Run()
    {
        static const std::string dy_contrib_prefix = "DY_MC";
        auto base_sub_category = EventSubCategory().SetCutResult(SelectionCut::mh, true)
                                                         .SetCutResult(SelectionCut::lowMET,true);
        const std::vector<int> ht_points = { 0, args.med_ht(), args.high_ht() };
        const std::vector<int> nJet_points = {args.low_nJet(), args.high_nJet() };
        static const size_t max_n_b = 2;
        if(fit_model == DYFitModel::NbjetBins) {
            subCategories = { base_sub_category };
            for(size_t nb = 0; nb <= max_n_b; ++nb) {
                const std::string name = boost::str(boost::format("%1%_%2%b") % dy_contrib_prefix % nb);
                contribution_names.push_back(name);
            }

        } else {
            subCategories = { EventSubCategory(base_sub_category).SetCutResult(SelectionCut::lowPt, true),
                              EventSubCategory(base_sub_category).SetCutResult(SelectionCut::medPt, true),
                              EventSubCategory(base_sub_category).SetCutResult(SelectionCut::highPt, true)};
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
            }
        }
        contribution_names.push_back("other_bkg_muMu");

        std::map<std::string, std::shared_ptr<RooRealVar>> scale_factor_map;
        for (const std::string& contrib_name: contribution_names){
            scale_factor_map[contrib_name] = std::make_shared<RooRealVar>
                        (("sf_"+contrib_name).c_str(),("Scale Factor for contribution "+contrib_name).c_str(),1,
                        args.scale_factor_range().min(),args.scale_factor_range().max());
        }
        std::string data_folder = "Data_SingleMuon";
        std::map<std::string,TH1*> dataCategories;
        std::map<std::string,std::shared_ptr<CategoryModel>> categories;

        RooCategory rooCategories("rooCategories","rooCategories") ;
        RooSimultaneous simPdf("simPdf","simultaneous pdf",rooCategories) ;

        for(const EventCategory& cat : eventCategories ){
            for (const EventSubCategory& subCategory : subCategories){
                EventAnalyzerDataId catId{cat,subCategory,EventRegion::SignalRegion(),EventEnergyScale::Central};
                std::string category = ToString(catId.Get<EventCategory>());
                std::string subcategory= ToString(catId.Get<EventSubCategory>());
                categories[category+subcategory] = std::make_shared<CategoryModel>(input_file,catId,contribution_names,
                                                                                   args.histo_name(),x,scale_factor_map);
                EventAnalyzerDataId dataId = catId.Set(data_folder);
                dataCategories[category+subcategory] = root_ext::ReadObject<TH1>(*input_file,dataId.GetName()+ "/" +
                                                                                 args.histo_name());

                rooCategories.defineType((category+subcategory).c_str()) ;
                simPdf.addPdf(*(categories[category+subcategory]->sum_pdf),(category+subcategory).c_str());
            }
        }

        // Construct combined dataset in (x,rooCategories)
        RooDataHist combData("combData","combined data",x,rooCategories,dataCategories);

        // Perform simultaneous fit
        RooFitResult* result = simPdf.fitTo(combData,Extended(kTRUE),Save()) ;

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
            cov_hist->GetXaxis()->SetBinLabel(i,contrib_name.c_str());
            cov_hist->GetYaxis()->SetBinLabel(i,contrib_name.c_str());
            cor_hist->GetXaxis()->SetBinLabel(i,contrib_name.c_str());
            cor_hist->GetYaxis()->SetBinLabel(i,contrib_name.c_str());

            scale_factors_hist->GetXaxis()->SetBinLabel(i,contrib_name.c_str());
            scale_factors_hist->SetBinContent(i,scale_factor_map[contrib_name]->getValV());
            scale_factors_hist->SetBinError(i,scale_factor_map[contrib_name]->getError());
            i++;
        }
        cov_hist->Write();
        cor_hist->Write();
        scale_factors_hist->Write();

        //Plotting
        for(const EventCategory& cat : eventCategories ){
            for(const EventSubCategory& sub_cat: subCategories){
                TCanvas* c = new TCanvas(("fit_"+ToString(cat)+"_"+ToString(sub_cat)).c_str(),
                                     ("fit in eventCategory " + ToString(cat) + ToString(sub_cat)).c_str(),800,400) ;
                RooPlot* frame = x.frame() ;
                combData.plotOn(frame,Cut((static_cast<std::string>("rooCategories==rooCategories::")+
                                    ToString(cat)+ToString(sub_cat)).c_str())) ;
                simPdf.plotOn(frame,Slice(rooCategories,(ToString(cat)+ToString(sub_cat)).c_str()),
                                  ProjWData(rooCategories,combData)) ;
                simPdf.plotOn(frame,Slice(rooCategories,(ToString(cat)+ToString(sub_cat)).c_str()),
                            Components((static_cast<std::string>("expdf_")+ToString(cat)+"_"+ToString(sub_cat)+
                                        static_cast<std::string>("_other_bkg_muMu")).c_str() ),
                            ProjWData(rooCategories,combData),LineStyle(kDashed)) ;
                gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(static_cast<float>(1.4)) ; frame->Draw() ;
                c->Write();
            }
        }
  }

private:
    Arguments args;
    std::shared_ptr<TFile> input_file, output_file;
    //X axis
    RooRealVar x;

    std::vector<EventCategory> eventCategories{EventCategory::TwoJets_ZeroBtag(),
                EventCategory::TwoJets_OneBtag(),EventCategory::TwoJets_TwoBtag()};
    std::vector<EventSubCategory> subCategories;

    std::vector<std::string> contribution_names;
    DYFitModel fit_model;
};

} // namesapce analysis
PROGRAM_MAIN(analysis::DY_estimation, Arguments) // definition of the main program function

