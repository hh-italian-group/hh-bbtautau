/*! Compute Mutual Information for each sample of mass or for ranges of masses.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <initializer_list>
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/StatEstimators.h"
#include "hh-bbtautau/Analysis/include/MvaVariables.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/ProgressReporter.h"
#include "AnalysisTools/Run/include/MultiThread.h"
#include "hh-bbtautau/Studies/include/MvaMethods.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "hh-bbtautau/Analysis/include/MvaConfigurationReader.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(std::string, optband_folder);
    REQ_ARG(unsigned, number_threads);
    REQ_ARG(bool, range);
    OPT_ARG(int, which_range, 0);
    OPT_ARG(Long64_t, number_events, 15000);
    OPT_ARG(bool, is_SM, false);
};

namespace analysis {
namespace mva_study{

using clock = std::chrono::system_clock;

class FindJSD {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    std::vector<ChannelSpin> set_SM{{"tauTau",SM_spin},{"muTau",SM_spin}, {"eTau",SM_spin},
                                 {"muTau",bkg_spin}, {"eTau",bkg_spin}, {"tauTau",bkg_spin}};
    std::vector<ChannelSpin> set_R{/*{"tauTau",0},*/ {"muTau",0}, {"eTau",0},
                                   {"tauTau",2}, {"muTau",2}, {"eTau",2},
                                 {"muTau",bkg_spin}, {"eTau",bkg_spin}, {"tauTau",bkg_spin}};
    std::vector<ChannelSpin> set;

    FindJSD(const Arguments& _args): args(_args), vars(1, 12345678,{}, {"channel", "mass", "spin"}), reporter(std::make_shared<TimeReporter>())
    {
        MvaSetupCollection setups;
        SampleEntryListCollection samples_list;

        ConfigReader configReader;
        MvaConfigReader setupReader(setups);
        configReader.AddEntryReader("SETUP", setupReader, true);
        SampleConfigReader sampleReader(samples_list);
        configReader.AddEntryReader("FILES", sampleReader, false);
        configReader.ReadConfig(args.cfg_file());

        set = args.is_SM() ? set_SM : set_R;
        samples = samples_list.at("Samples").files;
    }

    static bool IsInsideEllipse(double x, double y, double x0, double y0, double a, double b)
    {
        return pow(x - x0, 2) / pow(a, 2) + pow(y - y0, 2) / pow(b, 2) < 1.;
    }

    void LoadSkimmedData()
    {
        for (const auto& s: set){
            std::cout << s.channel << s.spin <<std::endl;
            for(const SampleEntry& entry:samples)
            {
                if ( entry.spin != s.spin) continue;
                auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
                auto tuple = ntuple::CreateEventTuple(s.channel, input_file.get(), true, ntuple::TreeState::Skimmed);
                Long64_t tot_entries = 0;
                for(const Event& event : *tuple) {
                    if(tot_entries >= args.number_events()) break;
                    LorentzVectorE_Float bb = event.jets_p4[0] + event.jets_p4[1];
                    if (!cuts::hh_bbtautau_2016::hh_tag::IsInsideMassWindow(event.SVfit_p4.mass(), bb.mass()))
                        continue;
                    if (entry.id == SampleType::Bkg_TTbar && event.file_desc_id>=2) continue;
                    if (entry.id == SampleType::Sgn_NonRes && event.file_desc_id!=0) continue;
                    auto eventInfoPtr =  analysis::MakeEventInfo(Parse<Channel>(s.channel) ,event) ;
                    EventInfoBase& eventbase = *eventInfoPtr;
//                    if (!IsInsideEllipse(eventbase.GetHiggsBB().GetMomentum().M(),eventbase.GetHiggsTTMomentum(false).M(),109.639, 87.9563, 43.0346,41.8451))
//                        continue;
                    vars.AddEvent(eventbase, entry.id, entry.spin, entry.weight);
                    tot_entries++;
                }
                std::cout << entry << " number of events: " << tot_entries << std::endl;
            }
            samples_mass[s] = vars.GetSampleVariables(s.channel, s.spin);
            TimeReport();
        }
    }

    SampleIdVarData LoadRangeData(const SampleIdVarData& samples_mass){
        SampleIdVarData samples_range;
        for (const auto& range : analysis::mva_study::ranges){
            if (!range.Contains(args.which_range())) continue;
            for (const auto& sample : samples_mass){
                if (!range.Contains(sample.first.mass)) continue;
                for (const auto& entry : sample.second){
                    for (const auto& value: entry.second){
                        samples_range[SampleId{SampleType::Sgn_Res, range.min()}][entry.first].push_back(value);
                    }
                }
            }
        }
        return samples_range;
    }


    void TimeReport(bool tot = false) const
    {
        reporter->TimeReport(tot);
    }

    void Run()
    {
        run::ThreadPull threads(args.number_threads());
        LoadSkimmedData();

        std::string file_name_prefix = "MutualInformationDistance";
        std::string file_name_prefix_bandwidth = "/OptimalBandwidth";
        if(args.range()) {
            std::cout<<"RANGES"<<std::endl;
            for (const auto& s: set){
                samples_range[s] = LoadRangeData(samples_mass[s]);
            }
            file_name_prefix += "Range";
            file_name_prefix_bandwidth += "Range";
        } else {
            samples_range = samples_mass;
        }

        for (const auto& s: set){
            std::cout<<"MID "<<s.channel<< "  " << s.spin<<std::endl;
            for (const auto& sample: samples_range[s]){
                std::cout<<"----"<<ToString(sample.first)<<"----"<<" entries: "<<sample.second.at("pt_l1").size()<<std::endl;
                std::stringstream ss;
                ss << std::fixed << std::setprecision(0) << s.spin;
                std::string spin = ss.str();
                bandwidth[s][sample.first] = Read_csvfile(args.optband_folder()+file_name_prefix_bandwidth+ToString(sample.first)+"_"+
                                                          s.channel+"_spin"+spin+"_newcut.csv");

                mutualmatrix[s][sample.first] = Mutual(sample.second, bandwidth.at(s).at(sample.first));
                std::ofstream ListMID(file_name_prefix+ToString(sample.first)+"_"+s.channel+"_spin"+spin+"_newcut.csv", std::ofstream::out);
                for(const auto& value: mutualmatrix[s][sample.first]){
                   for(const auto& var_name: value.first)
                   {
                       ListMID << var_name << "," ;
                   }
                   ListMID << value.second << std::endl;
                }
                TimeReport();
            }
        }

        TimeReport(true);
    }

private:
    Arguments args;
    SampleEntryCollection samples;
    MvaVariablesStudy vars;
    std::shared_ptr<TimeReporter> reporter;
    std::map<ChannelSpin,SampleIdVarData> samples_mass, samples_range;
    std::map<ChannelSpin,SampleIdNameElement> bandwidth, bandwidth_range, mutualmatrix, mutualmatrix_range;
    std::map<ChannelSpin, SamplePairNameNDElement> JSDivergenceSS, JSDivergenceSS_range;
};
}
}

PROGRAM_MAIN(analysis::mva_study::FindJSD, Arguments) // definition of the main program function
