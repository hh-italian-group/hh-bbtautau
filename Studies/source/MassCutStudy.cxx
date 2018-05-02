/*! Study of elliptical mass cut
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <random>
#include "Math/Minimizer.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/StatEstimators.h"
#include "hh-bbtautau/Analysis/include/MvaVariables.h"
#include "hh-bbtautau/Studies/include/MvaConfiguration.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"
#include "hh-bbtautau/Studies/include/MvaMethods.h"
#include "hh-bbtautau/Analysis/include/MvaConfigurationReader.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(int, spin);
    REQ_ARG(double, efficiency);
};

namespace analysis {
namespace mva_study{

using clock = std::chrono::system_clock;

class MassCutStudy {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    MassCutStudy(const Arguments& _args): args(_args), reporter(std::make_shared<TimeReporter>())
    {
        MvaSetupCollection setups;
        SampleEntryListCollection samples_list;

        ConfigReader configReader;
        MvaConfigReader setupReader(setups);
        configReader.AddEntryReader("SETUP", setupReader, true);
        SampleConfigReader sampleReader(samples_list);
        configReader.AddEntryReader("FILES", sampleReader, false);
        configReader.ReadConfig(args.cfg_file());

        samples = samples_list.at("Samples").files;
    }

    struct HH_candidates{
        double mbb;
        double mtt;
        double weight;

        HH_candidates(double _mbb, double _mtt, double _weight)
            : mbb(_mbb), mtt(_mtt), weight(_weight) {}
    };
    std::map<SampleId, std::vector<HH_candidates>> element;
    double N_bkg{0}, N_sgn{0};
    double max_histo_mbb, max_histo_mtt;

    void TimeReport(bool tot = false) const
    {
        reporter->TimeReport(tot);
    }

    double BkgMin (double x0, double y0, const double* param) const {
        static const double k = 1000;
        const double a = param[0];
        const double b = param[1];
        double n_b = 0, n_s = 0;

        if(a < 0 || b < 0 || x0 < 0 || y0 < 0) return 2;

        for (const auto& entry: element.at(SampleId::MassTot())){
            analysis::EllipseParameters ellipse_params{x0, a, y0, b};
            if(ellipse_params.IsInside(entry.mbb, entry.mtt))
                n_s += entry.weight;
        }
        for (const auto& entry: element.at(SampleId::Bkg())){
            analysis::EllipseParameters ellipse_params{x0, a, y0, b};
            if(ellipse_params.IsInside(entry.mbb, entry.mtt))
                n_b += entry.weight;
        }
        double result = n_b / N_bkg;
        if( args.efficiency() > n_s / N_sgn)
            result += 2*std::erf(k*(args.efficiency()-n_s/N_sgn));
//            result += std::exp(k*std::pow(args.efficiency()-n_s/N_sgn, 2)) - 1;

        std::cout <<"a: "<< a << " b: "<< b << " x0: " << x0 << " y0: " << y0
                  <<"  eff_s: "<<n_s/N_sgn<<"   eff_b"<<n_b / N_bkg<<"  result:"<< result << "\n";
        return result;
    }

    double SignalMax(double a, double b, const double* param) const
    {
        const double x0 = param[0];
        const double y0 = param[1];
        double n_s = 0;

        if(a < 0 || b < 0 || x0 < 0 || y0 < 0) return std::numeric_limits<double>::infinity();
        for (const auto& entry: element.at(SampleId::MassTot())){
            analysis::EllipseParameters ellipse_params{x0, a, y0, b};
            if(ellipse_params.IsInside(entry.mbb, entry.mtt))
                n_s += entry.weight;
        }
        std::cout <<"a: "<< a << " b: "<< b << " x0: " << x0 << " y0: " << y0
                  <<"  eff_s: "<<n_s/N_sgn << "\n";

        return -n_s / N_sgn;
    }


    void FindCenter(double a, double b, double& x0, double& y0, double step) const
    {
        auto min = std::shared_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit", "Combined"));
        min->SetMaxFunctionCalls(100000);
        min->SetMaxIterations(100000);
        min->SetPrintLevel(2);
//        min->SetStrategy(2);
        min->SetTolerance(0.1);
        min->SetPrecision(0.01);

        auto fn = std::bind(&MassCutStudy::SignalMax, this, a, b, std::placeholders::_1);
        ROOT::Math::Functor f(fn, 2);
        min->SetFunction(f);
        min->SetVariable(0,"max_mbb", x0, step);
        min->SetVariable(1,"max_mtt", y0, step);


        min->Minimize();

        x0 = min->X()[0];
        y0 = min->X()[1];
    }

    void FindRadius(double x0, double y0, double& a, double& b, double step) const
    {
        auto min = std::shared_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit", "Combined"));
        min->SetMaxFunctionCalls(100000);
        min->SetMaxIterations(100000);
        min->SetPrintLevel(2);
        min->SetTolerance(0.05);
        min->SetPrecision(0.01);

        auto fn = std::bind(&MassCutStudy::BkgMin, this, x0, y0, std::placeholders::_1);
        ROOT::Math::Functor f(fn, 2);
        min->SetFunction(f);
        min->SetVariable (0,"a(bb)", a, step);
        min->SetVariableLimits(0, 1, 100);
        min->SetVariable(1,"b(tt)", b, step);
        min->SetVariableLimits(1, 1, 150);


        min->Minimize();

        a = min->X()[0];
        b = min->X()[1];
    }


    void Run()
    {
        std::cout << ROOT::Math::MinimizerOptions::DefaultMinimizerType() << " "
                  << ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo() << std::endl;
        for(const auto& entry : samples)
        {
            if ( entry.id.IsSignal() && entry.spin != args.spin()) continue;
            auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
            auto tuple = ntuple::CreateEventTuple(args.tree_name(), input_file.get(), true, ntuple::TreeState::Skimmed);

            auto sampleid = entry.id.IsSignal() ? SampleId::MassTot() : SampleId::Bkg();

            for(const Event& event : *tuple) {
                auto eventInfoPtr =  analysis::MakeEventInfo(Parse<Channel>(args.tree_name()), event) ;
                EventInfoBase& eventbase = *eventInfoPtr;
//                if (!cuts::hh_bbtautau_2016::hh_tag::IsInsideMassWindow(eventbase.GetHiggsTTMomentum(true).mass(), eventbase.GetHiggsBB().GetMomentum().mass()))
//                    continue;
                const auto& Hbb = eventbase.GetHiggsBB().GetMomentum();
//                const auto& Htt_sv = eventbase.GetHiggsTTMomentum(true);
                const auto& Htt = eventbase.GetHiggsTTMomentum(false);
//                const auto& met = eventbase.GetMET().GetMomentum();
//                element[sampleid].emplace_back(Hbb.M(), (Htt_sv).M(), eventbase->weight_total);
//                element[sampleid].emplace_back(Hbb.M(), (Htt+met).M(), eventbase->weight_total);
                element[sampleid].emplace_back(Hbb.M(), Htt.M(), eventbase->weight_total);
                if (entry.id.IsBackground()) N_bkg += eventbase->weight_total;
                if (entry.id.IsSignal()) N_sgn += eventbase->weight_total;
            }
            std::cout << entry << " number of events: " << tuple->size() << std::endl;
        }
        TimeReport();
        auto histo_2d = std::make_shared<TH2D>("twod", "twod", 1000,0,500,1000,0,500);
        auto histo_mbb = std::make_shared<TH1D>("mbb", "mbb", 1000,0,500);
        auto histo_mtt = std::make_shared<TH1D>("mtt", "mtt", 1000,0,500);
        for (const auto& el : element[SampleId::MassTot()]){
            histo_mbb->Fill(el.mbb,el.weight);
            histo_mtt->Fill(el.mtt,el.weight);
            histo_2d->Fill(el.mtt,el.mbb,el.weight);
        }
        max_histo_mbb = histo_mbb->GetMaximumBin()*(500./histo_mbb->GetNbinsX());
        max_histo_mtt = histo_mtt->GetMaximumBin()*(500./histo_mtt->GetNbinsX());
        std::cout<<max_histo_mbb<<" "<<max_histo_mtt<<std::endl;
        auto outfile = root_ext::CreateRootFile("output_total.root");
        histo_2d->Write();
        histo_mbb->Write();
        histo_mtt->Write();

        std::cout<<element.size()<<std::endl;
        std::cout<<element[SampleId::MassTot()].size()<<"   "<<element[SampleId::Bkg()].size()<<std::endl;
        std::cout<<N_sgn<<"   "<<N_bkg<<std::endl;

        double a = histo_mbb->GetStdDev(), b = histo_mtt->GetStdDev();
        double x0 = histo_mbb->GetBinCenter(histo_mbb->GetMaximumBin());
        double y0 = histo_mtt->GetBinCenter(histo_mtt->GetMaximumBin());
        FindCenter(a, b, x0, y0, 0.1);
        std::cout << "x0: " << x0 << " y0: " << y0 << std::endl;
        /*a *= 2; */b *= 2;
        FindRadius(x0, y0, a, b, 0.01);
        std::cout << "a: " << a << " b: " << b << std::endl;

        std::cout << "Minimum:\n";
        double minimum[] = { a, b };
        BkgMin(x0, y0, minimum);

        std::cout << "AN:\n";
        const double AN[] = { 45., 35. };
        BkgMin(111., 116., AN);

        TimeReport(true);
    }

private:
    Arguments args;
    std::shared_ptr<TimeReporter> reporter;
    SampleEntryCollection samples;
};
}
}

PROGRAM_MAIN(analysis::mva_study::MassCutStudy, Arguments) // definition of the main program function
