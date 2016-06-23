/*! Code to randomly re-tag b-jets to correct MC.
Original code: https://github.com/ajgilbert/ICHiggsTauTau/blob/master/Analysis/Utilities/src/BTagWeight.cc
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <cmath>

#include <TRandom3.h>
#include <TH1.h>

#include "AnalyzerData.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisMath.h"
#include "exception.h"
#include "h-tautau/Analysis/include/Particles.h"
#include "h-tautau/Analysis/include/Htautau_2015.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "BTagCalibrationStandalone.h"

namespace analysis {
namespace btag {

// namespace Run2{

  struct jetEffInfo{
    float pt, eta, CSV;
    int   partonFlavour;
    float eff, SF;

    jetEffInfo(const ntuple::Event& event, int jetIndex):eff(0.),SF(0.){
        pt  = event.jets_p4.at(jetIndex).pt();
        eta = event.jets_p4.at(jetIndex).eta();
        CSV = event.jets_csv.at(jetIndex);
        partonFlavour = 0;/*event.partonFlavour_jets.at(jetIndex);*/
    }

  };

  typedef std::vector<jetEffInfo>    jetEffInfoVector;
  typedef std::pair<bool,jetEffInfo> bTagOutcome;
  typedef std::vector<bTagOutcome>   bTagMap;

  class BjetEffWeight {

  public:
    BjetEffWeight(const std::string& bTagEffName = "hh-bbtautau/Analysis/data/bTagEff_Loose.root",
                  const std::string& bjetSFName = "hh-bbtautau/Analysis/data/CSVv2.csv")
    :bTagEffFile(root_ext::OpenRootFile(bTagEffName))
    {
       using namespace btag_calibration;
       calib        = new BTagCalibration("CSVv2", bjetSFName);
       reader       = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "mujets", "central");
       reader_light = new BTagCalibrationReader(calib, BTagEntry::OP_MEDIUM, "incl", "central");

       eff_b = (TH2F*) bTagEffFile->Get("eff_b");
       eff_c = (TH2F*) bTagEffFile->Get("eff_c");
       eff_l = (TH2F*) bTagEffFile->Get("eff_l");
     }

     double Compute(const ntuple::Event& event){
           //std::cout<<" SyncTree Entry :   "<<current_entry<<"\n";
           using namespace cuts::Htautau_2015;

           bTagMap bTagMedium;
           for (int jetIndex = 0; jetIndex < event.jets_p4.size(); jetIndex++){
             jetEffInfo jetInfo(event, jetIndex);
             if      (abs (jetInfo.partonFlavour) == 5 ) {
               jetInfo.SF  = reader->eval(btag_calibration::BTagEntry::FLAV_B, jetInfo.eta, jetInfo.pt);
               jetInfo.eff = GetEfficiency(eff_b,jetInfo.pt,TMath::Abs(jetInfo.eta));
             }
             else if (abs (jetInfo.partonFlavour) == 4 ) {
               jetInfo.SF = reader->eval(btag_calibration::BTagEntry::FLAV_C, jetInfo.eta, jetInfo.pt);
               jetInfo.eff = GetEfficiency(eff_c,jetInfo.pt,TMath::Abs(jetInfo.eta));
             }
             else {
               jetInfo.SF = reader_light->eval(btag_calibration::BTagEntry::FLAV_UDSG, jetInfo.eta, jetInfo.pt);
               jetInfo.eff = GetEfficiency(eff_l,jetInfo.pt,TMath::Abs(jetInfo.eta));
             }

             bTagOutcome jetOutcome(false,jetInfo);
             if ( jetOutcome.second.CSV > cuts::Htautau_2015::btag::CSVM ) jetOutcome.first = true;
             bTagMedium.push_back(jetOutcome);

             // std::cout<<"\t  jet number "<<jetIndex<<" Flavour : "<<jetOutcome.second.partonFlavour
             //                                       <<" ---->   pt = "<<jetOutcome.second.pt
             //                                              <<" eta = "<<jetOutcome.second.eta
             //                                              <<" CSV = "<<jetOutcome.second.CSV
             //                                              <<" eff = "<<jetOutcome.second.eff
             //                                              <<"  SF = "<<jetOutcome.second.SF<<std::endl;
           }
           return GetBtagWeight(bTagMedium);
   }

    ~BjetEffWeight() {}

  private:
    double GetEfficiency (const TH2F* h, const double pt, const double eta){
      int xBin = h->GetXaxis()->FindFixBin(pt);
      int yBin = h->GetYaxis()->FindFixBin(eta);

      return h->GetBinContent(xBin,yBin);
    }

    double GetBtagWeight (const bTagMap map){
      double MC=1;
      double Data=1;
//      std::cout<<"  -------------   Weight Computation ------------ "<<std::endl;
      for (auto element : map){
        // std::cout<<"\t\t\t Flavour : "<<element.second.partonFlavour
        //                                       <<" ---->   pt = "<<element.second.pt
        //                                              <<" eta = "<<element.second.eta
        //                                              <<" CSV = "<<element.second.CSV
        //                                              <<" eff = "<<element.second.eff
        //                                              <<"  SF = "<<element.second.SF<<std::endl;

        if(element.first){
          MC=MC*element.second.eff;
          Data=Data*(element.second.eff*element.second.SF);
        }
        if(!(element.first)){
          MC=MC*(1-element.second.eff);
          Data=Data*(1-(element.second.eff*element.second.SF));
        }
      }
//      std::cout<<"   Data Weight = "<< Data <<"   MC Weight  = "<<MC<<std::endl;
      if (!MC) return 0;
      return Data/MC;
    }


  private:
    std::shared_ptr<TFile> bTagEffFile;
    btag_calibration::BTagCalibration *calib;
    btag_calibration::BTagCalibrationReader *reader, *reader_light;

    TH2F* eff_b = 0;
    TH2F* eff_c = 0;
    TH2F* eff_l = 0;
  };
//}


namespace Run1{
enum class payload { ALL2011, MORIOND2013, EPS13 };
enum class tagger { SSVHEM, SSVHPT, CSVM };

inline TH1F* SFb_error_2011()
{
    static const float SFb_error_2011_array[] = {
          0.0295675,
          0.0295095,
          0.0210867,
          0.0219349,
          0.0227033,
          0.0204062,
          0.0185857,
          0.0256242,
          0.0383341,
          0.0409675,
          0.0420284,
          0.0541299,
          0.0578761,
          0.0655432 };
    static const float ptbins_2011[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 670};
    static TH1F* SFb_error_2011_ = nullptr;
    if(!SFb_error_2011_) {
        const size_t nbins_2011 = sizeof(ptbins_2011) / sizeof(ptbins_2011[0]) - 1;
        SFb_error_2011_ = new TH1F("SFb_error_2011", "SFb_error_2011",nbins_2011, ptbins_2011);
        for(unsigned i=1; i<nbins_2011; i++)
            SFb_error_2011_->SetBinContent(i+1, SFb_error_2011_array[i-1]) ;
        //Adding number for 20-30 bin
        SFb_error_2011_->SetBinContent(1, 0.12);
    }
    return SFb_error_2011_;
}

inline TH1F* SFb_error_2012()
{
    static const float SFb_error_2012_array[] = {
          0.0415707,
          0.0204209,
          0.0223227,
          0.0206655,
          0.0199325,
          0.0174121,
          0.0202332,
          0.0182446,
          0.0159777,
          0.0218531,
          0.0204688,
          0.0265191,
          0.0313175,
          0.0415417,
          0.0740446,
          0.0596716 };
    static const float ptbins_2012[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
    static TH1F* SFb_error_2012_ = nullptr;
    if(!SFb_error_2012_) {
        const size_t nbins_2012 = sizeof(ptbins_2012) / sizeof(ptbins_2012[0]) - 1;
        SFb_error_2012_ = new TH1F("SFb_error_2012", "SFb_error_2012", nbins_2012, ptbins_2012);
        for(unsigned i = 0; i < nbins_2012; i++)
        SFb_error_2012_->SetBinContent(i+1, SFb_error_2012_array[i]) ;
    }
    return SFb_error_2012_;
}

inline double BEff(unsigned flavour, tagger algo, double pt, double eta)
{
    eta = std::abs(eta);
    if (eta >= 2.4)
        throw std::runtime_error("Jet eta is too large");

    const double x = std::min(pt, 670.0);
    if (flavour == 5) {
        if (algo == tagger::SSVHEM) return 0.64939;
        if (algo == tagger::SSVHPT) return 0.48337;
        if (algo == tagger::CSVM)   return 0.72973;
    } else if (flavour == 4) {
        if (algo == tagger::SSVHEM) return 0.17134;
        if (algo == tagger::SSVHPT) return 0.06522;
        if (algo == tagger::CSVM)   return 0.19249;
    } else {
        if (algo == tagger::SSVHEM) return (((-0.000420178+(0.00029105*x))+(-8.9398e-07*(x*x)))+(1.35401e-09*(x*(x*x))))+(-7.93826e-13*(x*(x*(x*x))));
        if (algo == tagger::SSVHPT) return (-2.9605e-05+(2.35624e-05*x))+(-1.77552e-08*(x*x));
        if (algo == tagger::CSVM)   {
            if (eta < 0.8)                return (0.00967751+(2.54564e-05*x))+(-6.92256e-10*(x*x));
            if (eta >= 0.8 && eta < 1.6)  return (0.00974141+(5.09503e-05*x))+(2.0641e-08*(x*x));
            if (eta >= 1.6 && eta < 2.4)  return (0.013595+(0.000104538*x))+(-1.36087e-08*(x*x));
        }
    }

    return 1.0;
}

inline double SF(payload set, unsigned flavour, tagger algo, double pt, double eta, int Btag_mode, int Bfake_mode)
{
    eta = std::abs(eta);
    if (eta >= 2.4)
        throw std::runtime_error("Jet eta is too large");

    double sf = 1.0;
    if (flavour == 5 || flavour == 4) {
        if (set == payload::ALL2011) {
            double x = std::max(30., std::min(pt, 670.0));
            if (algo == tagger::SSVHEM) {
                sf = 0.896462*((1.+(0.00957275*x))/(1.+(0.00837582*x)));
                if (Btag_mode > 0) {
                    double sf_err = 0.0;
                    static std::vector<double> unc_bins = {
                        40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 120.0,
                        160.0, 210.0, 260.0, 320.0, 400.0, 500.0, 670.0
                    };
                    static std::vector<double> unc_vals = {
                        0.0316234, 0.0310149, 0.02381, 0.0223228, 0.023461, 0.0202517, 0.0156249,
                        0.0214799, 0.0399369, 0.0416666, 0.0431031, 0.0663209, 0.0687731, 0.0793305
                    };
                    if (pt < 30.0) {
                        sf_err = 0.12;
                    } else if (pt >= 670.0) {
                        sf_err = 2.0 * unc_vals[13];
                    } else {
                        auto it = std::upper_bound(unc_bins.begin(), unc_bins.end(), pt);
                        if (it != unc_bins.end()) {
                            sf_err = unc_vals[unsigned(it - unc_bins.begin())];
                        }
                    }
                    if (flavour == 4) sf_err *= 2.0;
                    if (Btag_mode == 1) sf -= sf_err;
                    if (Btag_mode == 2) sf += sf_err;
                }
            }
            if (algo == tagger::SSVHPT) sf = 0.422556*((1.+(0.437396*x))/(1.+(0.193806*x)));
            if (algo == tagger::CSVM) {
            sf = 0.6981*((1.+(0.414063*x))/(1.+(0.300155*x)));
            //Adding uncertainty shifts for HiggsTauTau systematic study. Currently only implemented for CSVM.
            double err_SFb=SFb_error_2011()->GetBinContent(SFb_error_2011()->FindBin(pt));
            if(Btag_mode==1) sf=sf-err_SFb;
            if(Btag_mode==2) sf=sf+err_SFb;
            }
        } else if (set == payload::MORIOND2013) {
            double x = std::max(20., std::min(pt, 800.0));
            if (algo == tagger::CSVM)   sf = 0.726981*((1.+(0.253238*x))/(1.+(0.188389*x)));
        } else if (set == payload::EPS13) {
            double x = std::max(20., std::min(pt, 800.0));
            if (algo == tagger::CSVM) {
                sf = (0.938887+(0.00017124*x))+(-2.76366e-07*(x*x));
                //Adding uncertainty shifts for HiggsTauTau systematic study. Currently only implemented for CSVM.
                double err_SFb=SFb_error_2012()->GetBinContent(SFb_error_2012()->FindBin(pt));
                if(Btag_mode==1) sf=sf-err_SFb;
                if(Btag_mode==2) sf=sf+err_SFb;
            }
        }
    } else {
        if (set == payload::ALL2011) {
            double x = std::max(20., std::min(pt, 670.0));
            if (algo == tagger::SSVHEM)     sf = ((0.890254+(0.000553319*x))+(-1.29993e-06*(x*x)))+(4.19294e-10*(x*(x*x)));
            if (algo == tagger::SSVHPT)     sf = ((0.97409+(0.000646241*x))+(-2.86294e-06*(x*x)))+(2.79484e-09*(x*(x*x)));
            if (algo == tagger::CSVM) {
                if(Bfake_mode==0) {
                    if (eta < 0.8)                sf = ((1.06182+(0.000617034*x))+(-1.5732e-06*(x*x)))+(3.02909e-10*(x*(x*x)));
                    if (eta >= 0.8 && eta < 1.6)  sf = ((1.111+(-9.64191e-06*x))+(1.80811e-07*(x*x)))+(-5.44868e-10*(x*(x*x)));
                    if (eta >= 1.6 && eta < 2.4)  sf = ((1.08498+(-0.000701422*x))+(3.43612e-06*(x*x)))+(-4.11794e-09*(x*(x*x)));
                }
                //Adding uncertainty shifts for HiggsTauTau systematic study. Currently only implemented for CSVM.
                if(Bfake_mode==1) {
                    if (eta < 0.8)                sf = ((0.972455+(7.51396e-06*x))+(4.91857e-07*(x*x)))+(-1.47661e-09*(x*(x*x)));
                    if (eta >= 0.8 && eta < 1.6)  sf = ((1.02055+(-0.000378856*x))+(1.49029e-06*(x*x)))+(-1.74966e-09*(x*(x*x)));
                    if (eta >= 1.6 && eta < 2.4)  sf = ((0.983476+(-0.000607242*x))+(3.17997e-06*(x*x)))+(-4.01242e-09*(x*(x*x)));
                }
                if(Bfake_mode==2) {
                    if (eta < 0.8)                sf = ((1.15116+(0.00122657*x))+(-3.63826e-06*(x*x)))+(2.08242e-09*(x*(x*x)));
                    if (eta >= 0.8 && eta < 1.6)  sf = ((1.20146+(0.000359543*x))+(-1.12866e-06*(x*x)))+(6.59918e-10*(x*(x*x)));
                    if (eta >= 1.6 && eta < 2.4)  sf = ((1.18654+(-0.000795808*x))+(3.69226e-06*(x*x)))+(-4.22347e-09*(x*(x*x)));
                }
            }
        } else if (set == payload::MORIOND2013) {
            double x = std::max(20., std::min(pt, 800.0));
            if (algo == tagger::CSVM) {
                if (eta < 0.8)                sf = ((1.06238+(0.00198635*x))+(-4.89082e-06*(x*x)))+(3.29312e-09*(x*(x*x)));
                if (eta >= 0.8 && eta < 1.6)  sf = ((1.08048+(0.00110831*x))+(-2.96189e-06*(x*x)))+(2.16266e-09*(x*(x*x)));
                if (eta >= 1.6 && eta < 2.4)  sf = ((1.09145+(0.000687171*x))+(-2.45054e-06*(x*x)))+(1.7844e-09*(x*(x*x)));
            }
        } else if (set == payload::EPS13) {
            double x = std::max(20., std::min(pt, 1000.0));
            if (algo == tagger::CSVM) {
                if(Bfake_mode==0) {
                    if (eta < 0.8)                sf = ((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)));
                    if (eta >= 0.8 && eta < 1.6)  sf = ((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)));
                    if (eta >= 1.6 && eta < 2.4)  sf = ((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)));
                }
                //Adding uncertainty shifts for HiggsTauTau systematic study. Currently only implemented for CSVM.
                if(Bfake_mode==1) {
                    if (eta < 0.8)                sf = ((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)));
                    if (eta >= 0.8 && eta < 1.6)  sf = ((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)));
                    if (eta >= 1.6 && eta < 2.4)  sf = ((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)));
                }
                if(Bfake_mode==2) {
                    if (eta < 0.8)                sf = ((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)));
                    if (eta >= 0.8 && eta < 1.6)  sf = ((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)));
                    if (eta >= 1.6 && eta < 2.4)  sf = ((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)));
                }
            }
        }
    }
    return sf;
}

}

} // btag
} // analysis
