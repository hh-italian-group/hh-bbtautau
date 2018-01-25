/*
** class  : HHReweight5D
** author : L. Cadamuro (LLR)
** date   : 17 Feb 2017
** brief  : Compute a reweight coefficient to transform an input distribution into a general point in kl, kt plane. 31/05/2017: Updated to 5D reweight
** link   : https://github.com/LLRCMS/KLUBAnalysis/blob/7475ab49acb704a912104ca9e134319e3e479e58/interface/HHReweight5D.h
**          https://github.com/LLRCMS/KLUBAnalysis/blob/master/src/HHReweight5D.cc
*/

#include "TH2D.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TFile.h"
#include "AnalysisTools/Core/include/exception.h"


#define NCOEFFSA 15
#define DEBUG false

struct BenchmarkParameters{
    double kl, kt, c2, cg, c2g;
    BenchmarkParameters(int_kl, double _kt, double _c2, double _cg, double _c2g):
        kl(_kl), kt(_kt), c2(_c2), cg(_cg), c2g(_c2g) {}
};

//static const  std::map<int, BenchmarkParameters> BenchmarkId ={{-1,{0.,1.,0.,0.,0.}}, {0,{1.,1.,0.,0.,0.}}, {1, {7.5, 1.0,-1.,0.,0.}},
//                                                               {2,{1.,1.,0.5,-0.8,0.6}}, {3,{1.,1.,-1.5,0.,-0.8}}, {4,{-3.5,1.5,-3.,0.,0.}},
//                                                               {5,{1.,1.,0.,0.8,-1.}}, {6,{2.4,1.0,0.,0.2,-0.2}}, {7,{5.,1.,0.,0.2,-0.2}},
//                                                               {8,{15.,1.,0.,-1.,1.}}, {9,{1.,1.,1.,-0.6,0.6}}, {10,{10,1.5,-1.0,0.,0.,}},
//                                                               {11,{2.4,1.,0.,1.,-1.}}, {12, {15.,1.,1.,0.,0.}}};

//static const  std::map<int, BenchmarkParameters> BenchmarkId ={{-1,{0.,1.,0.,0.,0.}}, {11,{1.,1.,0.,0.,0.}}, {5, {7.5, 1.0,-1.,0.,0.}},
//                                                               {8,{1.,1.,0.5,-0.8,0.6}}, {10,{1.,1.,-1.5,0.,-0.8}}, {9,{-3.5,1.5,-3.,0.,0.}},
//                                                               {7,{1.,1.,0.,0.8,-1.}}, {4,{2.4,1.0,0.,0.2,-0.2}}, {1,{5.,1.,0.,0.2,-0.2}},
//                                                               {6,{15.,1.,0.,-1.,1.}}, {12,{1.,1.,1.,-0.6,0.6}}, {0,{10,1.5,-1.0,0.,0.,}},
//                                                               {2,{2.4,1.,0.,1.,-1.}}, {3, {15.,1.,1.,0.,0.}}};

static const  std::map<int, BenchmarkParameters> BenchmarkId ={{0,{-20,1.,0.,0.,0.}}, {1,{-10.,1.,0.,0.,0.}}, {2,{-5,1.,0.,0.,0.}},
                                                               {3,{-1.,1.,0.,0.,0.}}, {4,{0,1.,0.,0.,0.}}, {5,{1.,1.,0.,0.,0.}},
                                                               {6,{2.,1.,0.,0.,0.}}, {7,{5.,1.,0.,0.,0.}}, {8,{10.,1.,0.,0.,0.}},
                                                               {9,{30.,1.,0.,0.,0.}}};

class HHReweight5D{
    public:

    HHReweight5D() {}
    HHReweight5D(std::string _coeffFile, const  TH2* _hInput, bool _useAbsEta=true): coeffFile(_coeffFile), useAbsEta(_useAbsEta)
    {
        readInputFile(coeffFile); // initialize the reweight parameters from the txt input

        TH2* cloneH = dynamic_cast<TH2*>(_hInput->Clone("h_input"));
        if (!CheckConsistency(cloneH, h_A_vec.at(0).get()))
        {
            throw analysis::exception("* Error : the input histogram to HHReweight is not compatible with the reweight file, did you use the correct binning?");
        }
        hInput.reset(cloneH);

        A_13TeV = {2.09078, 10.1517, 0.282307, 0.101205, 1.33191, -8.51168, -1.37309, 2.82636, 1.45767, -4.91761, -0.675197, 1.86189, 0.321422, -0.836276, -0.568156};

        useAbsEta = _useAbsEta;
    }

        // return the weight to be applied for the reweight
        // NOTE: normalization is arbitrary you'll have to scale by the sum of weights
        double getWeight(const BenchmarkParameters& parameters, double mhh, double cth)
        {
            double kl = parameters.kl;
            double kt = parameters.kt;
            double c2 = parameters.c2;
            double cg = parameters.cg;
            double c2g = parameters.c2g;

            if (useAbsEta) cth = TMath::Abs(cth);

            std::pair<int,int> bins = find2DBin(hInput.get(), mhh, cth);
            double denom = hInput->GetBinContent(bins.first, bins.second);
            if (denom == 0)
                return 0;

            double nEvSM = h_SM->GetBinContent(bins.first, bins.second);

            std::array<double, NCOEFFSA> Acoeffs;
            for (size_t ic = 0; ic < NCOEFFSA; ++ic)
            {
              Acoeffs[ic] = (h_A_vec.at(ic))->GetBinContent(bins.first, bins.second);
            }
            double effBSM = nEvSM * functionGF(kl,kt,c2,cg,c2g,Acoeffs)/functionGF(kl,kt,c2,cg,c2g,A_13TeV);

            if (DEBUG && effBSM/denom < 0){
              std::cout << "** HHReweight5D : warning : I am getting a negative weight "
                << "kl, kt, c2, cg, c2g, mhh, cth " << kl << " " << kt << " " << c2 << " " << cg << " " << c2g << " " << mhh << " " << cth << " | vals: "
                << nEvSM << " " << functionGF(kl,kt,c2,cg,c2g,Acoeffs) << " " << functionGF(kl,kt,c2,cg,c2g,A_13TeV) << " " << denom
                << std::endl;
            }

            if (effBSM/denom < 0) return 0; // sometimes I get negative coeffs.. this should be a temporary fix!
            return (effBSM/denom) ;
        }

        double getWeight(double kl, double kt, double mhh, double cth) // kl, kt only
        {
            BenchmarkParameters parameters(kl, kt, 0, 0, 0);
            return getWeight(parameters, mhh, cth);
        }

    private:
        std::shared_ptr<TH2> hInput;
        std::string coeffFile;
        bool useAbsEta;

        std::array<std::shared_ptr<TH2>, NCOEFFSA> h_A_vec;
        std::array<double, NCOEFFSA> A_13TeV;


        // the coefficients of the reweight - read from the input file
        std::shared_ptr<TH2> h_SM;
        std::shared_ptr<TH2> h_sumv1;

        void readInputFile(const std::string& coeffFile)
        {
            if (DEBUG) std::cout << " -- Reading input file" << coeffFile << std::endl;

            // create histograms to be filled from file
            // this is the binning of input file histogram
            double binning_mHH [56] = { 250,260,270,280,290,300,310,320,330,340,
                                        350,360,370,380,390,400,410,420,430,440,
                                        450,460,470,480,490,
                                        500,510,520,530,540,550,600,610,620,630,
                                        640,650,660,670,680,690,700,750,800,850,
                                        900,950,1000,1100,1200,1300,1400,1500.,1750,2000,50000};
            double binning_cth [5]  = {0.0, 0.4, 0.6, 0.8, 1.0} ;
            int nbins_mHH = 55; // size of arrays - 1
            int nbins_cth = 4;  // size of arrays - 1

            for (size_t ic = 0; ic < NCOEFFSA; ++ic)
            {
              std::string name = "h_A" + std::to_string(ic);
              h_A_vec.at(ic) = std::make_shared<TH2D> (name.c_str(), name.c_str(), nbins_mHH, binning_mHH, nbins_cth, binning_cth );
            }

            h_SM    = std::make_shared<TH2D> ("h_SM", "h_SM",       nbins_mHH, binning_mHH, nbins_cth, binning_cth );
            h_sumv1 = std::make_shared<TH2D> ("h_sumv1", "h_sumv1", nbins_mHH, binning_mHH, nbins_cth, binning_cth );

            if (DEBUG) std::cout << " -- Histograms done" <<std::endl;

            // read and fill from the file
            std::ifstream infile;
            infile.open(coeffFile);
            if (!infile.is_open())
                throw analysis::exception("Could not open input file");

            std::string line;
            while (std::getline(infile, line))
            {
                if (DEBUG) std::cout << " -- Reading line " << line << std::endl;
                line = line.substr(0, line.find("#", 0)); // remove comments introduced by #
                if (!line.empty())
                {
                    std::vector<std::string> tokens = tokenize(line);
                    if (tokens.size() != 35)
                    {
                        throw analysis::exception  ("Error in reading input file: cannot interpret line: '%1%'") % line;
                    }
                    //The columns are respectively: nbins GenMhh GenCostStar NenventsSM NenventsSumV1 A1 A3 A7 errorA1 errorA3 errorA7

                    double mHH = std::stod(tokens.at(1));
                    double cth = std::stod(tokens.at(2));

                    int ibin = h_A_vec.at(0)->FindBin(mHH, cth);
                    h_SM->SetBinContent(ibin, std::stod(tokens.at(3)));
                    h_sumv1->SetBinContent(ibin, std::stod(tokens.at(4)));

                    for (size_t ic = 0; ic < NCOEFFSA; ++ic)
                    {
                      (h_A_vec.at(ic))->SetBinContent(ibin, std::stod(tokens.at(ic+5)) );
                      (h_A_vec.at(ic))->SetBinError(ibin, std::stod(tokens.at(ic+5+15)) );
                    }

                    if (DEBUG)
                    {
                        std::cout << " -- I'll store a file with the histograms" << std::endl;
                        TFile* fOut = TFile::Open("HHReweight_histograms.root", "recreate");
                        fOut->cd();
                        for (size_t ic = 0; ic < NCOEFFSA; ++ic)
                        {
                          (h_A_vec.at(ic))->Write();
                        }
                        h_SM->Write();
                        h_sumv1->Write();
                        fOut->Close();
                    }
                }
            }
        }

        // double functionGF(double kl, double kt, double c2, double cg, double c2g, double A1, double A3, double A7);
        double functionGF(double kl, double kt, double c2, double cg, double c2g, std::array<double, NCOEFFSA> const &A)
        {
            // this can be extended to 5D coefficients; currently c2, cg, c2g are unused
            // return ( A1*pow(kt,4) + A3*pow(kt,2)*pow(kl,2) + A7*kl*pow(kt,3) );
            return ( A[0]*pow(kt,4) + A[1]*pow(c2,2) + (A[2]*pow(kt,2) + A[3]*pow(cg,2))*pow(kl,2) + A[4]*pow(c2g,2) + ( A[5]*c2 + A[6]*kt*kl )*pow(kt,2) + (A[7]*kt*kl + A[8]*cg*kl )*c2 + A[9]*c2*c2g + (A[10]*cg*kl + A[11]*c2g)*pow(kt,2)+ (A[12]*kl*cg + A[13]*c2g )*kt*kl + A[14]*cg*c2g*kl );
        }


        // split a string on whitespaces, return elements
        std::vector<std::string> tokenize(std::string input)
        {
            if (DEBUG) std::cout << " -- Tokenizing input " << input << std::endl;

            std::istringstream buffer(input);
            std::vector<std::string> ret((std::istream_iterator<std::string>(buffer)), std::istream_iterator<std::string>());

            if (DEBUG){
                std::cout << "I got " << ret.size() << " tokens" << std::endl;
                std::cout << "     --> " ;
                for (std::string x : ret) std::cout << ":" << x << ": ";
                std::cout << std::endl;
            }

            return ret;
        }


        // return bin in 2D isto wihtout under/over flow (e.g. if ibin > ibinmax , ibin = ibinmax)
        std::pair<int,int> find2DBin(TH2* h, double x, double y)
        {
            int ibinx = h->GetXaxis()->FindBin(x);
            int ibiny = h->GetYaxis()->FindBin(y);

            if (ibinx <= 0) ibinx = 1;
            if (ibinx > h->GetNbinsX()) ibinx = h->GetNbinsX();

            if (ibiny <= 0) ibiny = 1;
            if (ibiny > h->GetNbinsY()) ibiny = h->GetNbinsY();

            return std::make_pair(ibinx, ibiny);

        }

        // adapted from ROOT to check histogram consistency
        bool CheckConsistency(const TH2* h1, const TH2* h2)
        {
           if (h1 == h2) return true;

           if (h1->GetDimension() != h2->GetDimension() ) {
              // throw DifferentDimension();
              return false;
           }
           Int_t dim = h1->GetDimension();

           // returns kTRUE if number of bins and bin limits are identical
           Int_t nbinsx = h1->GetNbinsX();
           Int_t nbinsy = h1->GetNbinsY();
           Int_t nbinsz = h1->GetNbinsZ();

           // Check whether the histograms have the same number of bins.
           if (nbinsx != h2->GetNbinsX() ||
               (dim > 1 && nbinsy != h2->GetNbinsY())  ||
               (dim > 2 && nbinsz != h2->GetNbinsZ()) ) {
              // throw DifferentNumberOfBins();
              return false;
           }

           bool ret = true;

           // check axis limits
          ret &= CheckAxisLimits(h1->GetXaxis(), h2->GetXaxis());
          if (dim > 1) ret &= CheckAxisLimits(h1->GetYaxis(), h2->GetYaxis());
          if (dim > 2) ret &= CheckAxisLimits(h1->GetZaxis(), h2->GetZaxis());

          // check bin limits
          ret &= CheckBinLimits(h1->GetXaxis(), h2->GetXaxis());
          if (dim > 1) ret &= CheckBinLimits(h1->GetYaxis(), h2->GetYaxis());
          if (dim > 2) ret &= CheckBinLimits(h1->GetZaxis(), h2->GetZaxis());


           return ret;
        }

        bool CheckAxisLimits(const TAxis *a1, const TAxis *a2 )
        {
           if ( ! TMath::AreEqualRel(a1->GetXmin(), a2->GetXmin(),1.E-12) ||
                ! TMath::AreEqualRel(a1->GetXmax(), a2->GetXmax(),1.E-12) ) {
              // throw DifferentAxisLimits();
              return false;
           }
           return true;
        }

        bool CheckBinLimits(const TAxis* a1, const TAxis * a2)
        {
           const TArrayD * h1Array = a1->GetXbins();
           const TArrayD * h2Array = a2->GetXbins();
           Int_t fN = h1Array->fN;
           if ( fN != 0 ) {
              if ( h2Array->fN != fN ) {
                 // throw DifferentBinLimits();
                 return false;
              }
              else {
                 for ( int i = 0; i < fN; ++i ) {
                    if ( ! TMath::AreEqualRel( h1Array->GetAt(i), h2Array->GetAt(i), 1E-10 ) ) {
                       // throw DifferentBinLimits();
                       return false;
                    }
                 }
              }
           }

           return true;
        }

};


// double HHReweight5D::functionGF(double kl, double kt, double c2, double cg, double c2g, double A1, double A3, double A7)
// {
//     // this can be extended to 5D coefficients; currently c2, cg, c2g are unused
//     return ( A1*pow(kt,4) + A3*pow(kt,2)*pow(kl,2) + A7*kl*pow(kt,3) );
// }






