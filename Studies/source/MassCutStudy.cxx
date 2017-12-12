/*! Study of elliptical mass cut
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <random>
#include "Math/Minimizer.h"
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

    double BkgMin (const double* param) const {
        static const double k = 1000;
        const double a = param[0];
        const double b = param[1];
        double n_b = 0, n_s = 0;

        for (const auto& entry: element.at(SampleId::MassTot())){
            auto temp = pow(entry.mbb-param[2],2)/pow(a,2) + pow(entry.mtt-param[3],2)/pow(b,2);
            if (temp<1) n_s += entry.weight;
        }
        for (const auto& entry: element.at(SampleId::Bkg())){
            auto temp = pow(entry.mbb-param[2],2)/pow(a,2) + pow(entry.mtt-param[3],2)/pow(b,2) ;
            if (temp<1) n_b += entry.weight;
        }
        double result = n_b / N_bkg;
        if( args.efficiency() > n_s / N_sgn)
            result += std::exp(k*std::pow(args.efficiency()-n_s/N_sgn, 2)) - 1;

        std::cout<<"a: "<<a<<" b: "<<b<<"  eff_s: "<<n_s/N_sgn<<"   eff_b"<<n_b / N_bkg<<"  result:"<<result<<std::endl;
        return result;
    }

    void Run()
    {
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
                const auto& Htt_sv = eventbase.GetHiggsTTMomentum(true);
                const auto& Htt = eventbase.GetHiggsTTMomentum(false);
                const auto& met = eventbase.GetMET().GetMomentum();
                element[sampleid].emplace_back(Hbb.M(), (Htt_sv).M(), eventbase->weight_total);
//                element[sampleid].emplace_back(Hbb.M(), (Htt+met).M(), eventbase->weight_total);
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

        ROOT::Math::Minimizer* min =ROOT::Math::Factory::CreateMinimizer("", "");
        min->SetMaxFunctionCalls(1000000);
        min->SetMaxIterations(100000);
        min->SetTolerance(0.01);
        min->SetPrintLevel(2);

        auto fn = std::bind(&MassCutStudy::BkgMin, this, std::placeholders::_1);
        ROOT::Math::Functor f(fn,4);
        double step[4] = {0.1,0.1, 0.01,0.01};
        double variable[4] = {150, 150, max_histo_mtt, max_histo_mbb};
        min->SetFunction(f);
        min->SetVariable(0,"a(bb)",variable[0], step[0]);
        min->SetVariable(1,"b(tt)",variable[1], step[1]);
        min->SetVariable(2,"max_mbb",variable[2], step[2]);
        min->SetVariable(3,"max_mtt",variable[3], step[3]);

        min->Minimize();

        const double *xs = min->X();
        std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << ","<< xs[2] << "," << xs[3] << "): "<< min->MinValue()  << std::endl;

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
