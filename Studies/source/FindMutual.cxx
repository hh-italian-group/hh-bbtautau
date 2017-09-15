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
    OPT_ARG(Long64_t, number_events, 5000);
};

namespace analysis {
namespace mva_study{

using clock = std::chrono::system_clock;

namespace {
Range<int> low_mass(250, 320), medium_mass(340, 400), high_mass(450,900);
std::vector<Range<int>> ranges{low_mass, medium_mass, high_mass};
}

class FindJSD {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

        std::vector<ChannelSpin> set{{"muTau",0},{"eTau",0},{"tauTau",0},{"muTau",2},{"eTau",2},{"tauTau",2},{"muTau",1},{"eTau",1},{"tauTau",1},{"muTau",-1},{"eTau",-1},{"tauTau",-1}};

    std::map<ChannelSpin,SampleIdVarData> samples_mass, samples_range;
    std::map<ChannelSpin,SampleIdNameElement> bandwidth, bandwidth_range, mutualmatrix, mutualmatrix_range;
    std::map<ChannelSpin, SamplePairNameNDElement> JSDivergenceSS, JSDivergenceSS_range;

    FindJSD(const Arguments& _args): args(_args), vars(1, 12345678,{}, {"channel", "mass"}), reporter(std::make_shared<TimeReporter>())
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

    void LoadSkimmedData()
    {
        for (const auto& s: set){
            std::cout << s.first << s.second <<std::endl;
            for(const SampleEntry& entry:samples)
            {
                if ( entry.spin != s.second) continue;
                auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
                auto tuple = ntuple::CreateEventTuple(s.first, input_file.get(), true, ntuple::TreeState::Skimmed);
                Long64_t tot_entries = 0;
                for(const Event& event : *tuple) {
                    if(tot_entries >= args.number_events()) break;
                    if (entry.id == SampleType::Bkg_TTbar && event.file_desc_id>=2) continue;
                    if (entry.id == SampleType::Sgn_NonRes && event.file_desc_id!=0) continue;
                    vars.AddEvent(event, entry.id, entry.spin, s.first, entry.weight);
                    tot_entries++;
                }
                std::cout << entry << " number of events: " << tot_entries << std::endl;
            }
            samples_mass[s] = vars.GetSampleVariables(s.first, s.second);
            TimeReport();
        }
    }

    SampleIdVarData LoadRangeData(const SampleIdVarData& samples_mass){
        SampleIdVarData samples_range;
        for (const auto& range : ranges){
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

    void Run()
    {
        run::ThreadPull threads(args.number_threads());
        LoadSkimmedData();

        if (!args.range()){
            for (const auto& s: set){
                std::cout<<std::endl<<s.first<< "  " << s.second<<std::endl;
                for (const auto& sample: samples_mass[s]){
                    std::cout<<"----"<<ToString(sample.first)<<"----"<<" entries: "<<sample.second.at("pt_l1").size()<<std::endl;
                    std::stringstream ss;
                    ss << std::fixed << std::setprecision(0) << s.second;
                    std::string spin = ss.str();
                    bandwidth[s][sample.first] = Read_csvfile(args.optband_folder()+"/OptimalBandwidth"+ToString(sample.first)+"_"+s.first+"_spin"+spin+".csv");
                    mutualmatrix[s][sample.first] = Mutual(sample.second, bandwidth.at(s).at(sample.first));
                    std::ofstream ListMID("MutualInformationDistance"+ToString(sample.first)+"_"+s.first+"_spin"+spin+".csv", std::ofstream::out);
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
        }
        else{
            std::cout<<"RANGES"<<std::endl;
            for (const auto& s: set){
                samples_range[s]=LoadRangeData(samples_mass[s]);
            }
            for (const auto& s: samples_range){
                std::cout<<std::endl<<s.first.first<< "  " << s.first.second<<std::endl;
                for (const auto& sample: s.second){
                    std::cout<<"----Range"<<ToString(sample.first)<<"----"<<" entries: "<<sample.second.begin()->second.size()<<std::endl;
                    if (sample.second.size()<2) continue;
                    std::stringstream ss;
                    ss << std::fixed << std::setprecision(0) << s.first.second;
                    std::string spin = ss.str();
                    bandwidth_range[s.first][sample.first] = Read_csvfile(args.optband_folder()+"/OptimalBandwidthRange"+ToString(sample.first)+"_"+s.first.first+"_spin"+spin+".csv");
                    mutualmatrix_range[s.first][sample.first] = Mutual(sample.second,bandwidth_range.at(s.first).at(sample.first));

                    std::ofstream ListMID("MutualInformationDistanceRange"+ToString(sample.first)+"_"+s.first.first+"_spin"+spin+".csv", std::ofstream::out);
                    for(const auto& value: mutualmatrix_range[s.first][sample.first]){
                       for(const auto& var_name: value.first)
                       {
                           ListMID << var_name << "," ;
                       }
                       ListMID << value.second << std::endl;
                    }
                    TimeReport();
                }
            }
        }



        TimeReport(true);
    }


    void TimeReport(bool tot = false) const
    {
        reporter->TimeReport(tot);
    }

private:
    Arguments args;
    SampleEntryCollection samples;
    MvaVariablesStudy vars;
    std::shared_ptr<TimeReporter> reporter;
};
}
}

PROGRAM_MAIN(analysis::mva_study::FindJSD, Arguments) // definition of the main program function
