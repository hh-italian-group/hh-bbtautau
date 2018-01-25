/*! Compute Jensen Shannon Divergence for each sample of mass or for ranges of masses.
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
    REQ_ARG(std::string, suffix);
    REQ_ARG(unsigned, number_threads);
    REQ_ARG(bool, range);
    OPT_ARG(int, which_range, 0);
    OPT_ARG(Long64_t, number_events, 15000);
    OPT_ARG(bool, is_SM, false);
    OPT_ARG(bool, bkg_vs_sgn, true);
    REQ_ARG(std::string, channel);
    REQ_ARG(int, spin);

};

namespace analysis {
namespace mva_study{

using clock = std::chrono::system_clock;

class FindJSD {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    std::vector<ChannelSpin> set_SM{{"tauTau",SM_spin}, {"muTau",SM_spin}, {"eTau",SM_spin},
                                 {"muTau",bkg_spin}, {"eTau",bkg_spin}, {"tauTau",bkg_spin}};
    std::vector<ChannelSpin> set_R{{"tauTau",0}, {"muTau",0}, {"eTau",0},
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

        set.emplace_back(args.channel(), args.spin());
        set.emplace_back(args.channel(), bkg_spin);
        std::cout<<"set size "<<set.size()<<std::endl;
//        set = args.is_SM() ? set_SM : set_R;
        samples = samples_list.at("Samples").files;
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
                    if (args.suffix() == "_ANcut"){
                        if (!cuts::hh_bbtautau_2016::hh_tag::m_hh_window().IsInside(event.SVfit_p4.mass(),bb.mass())) continue;
                    }
                    if (entry.id == SampleType::Bkg_TTbar && event.file_desc_id>=2) continue;
                    if (entry.id == SampleType::Sgn_NonRes && event.file_desc_id!=0) continue;
                    auto eventInfoPtr =  analysis::MakeEventInfo(Parse<Channel>(s.channel) ,event) ;
                    EventInfoBase& eventbase = *eventInfoPtr;
                    if (args.suffix() == "_newcut"){
                        if (!cuts::hh_bbtautau_2016::hh_tag::new_m_hh_window().IsInside(eventbase.GetHiggsTTMomentum(false).M(),bb.mass()))
                            continue;

                    }
                    vars.AddEvent(eventbase, entry.id, entry.spin, entry.weight);
                    tot_entries++;
                }
                std::cout << entry << " number of events: " << tuple->size() << std::endl;
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

    SamplePairNameNDElement JensenDivergenceSS(const ChannelSpin& chsp){
       SamplePairNameNDElement JSDss_range;
       for(auto sample_mass1 = samples_mass.at(chsp).begin(); sample_mass1 != samples_mass.at(chsp).end(); ++sample_mass1) {
           if ( sample_mass1->first.IsBackground()) continue;
           for(auto sample_mass2 = sample_mass1; sample_mass2 != samples_mass.at(chsp).end(); ++sample_mass2) {
               if ( sample_mass2->first.IsBackground()) continue;
               SamplePair mass_pair(sample_mass1->first, sample_mass2->first);
               JSDss_range[mass_pair] = JensenDivergenceSamples(samples_mass.at(chsp).at(sample_mass1->first),
                                                                samples_mass.at(chsp).at(sample_mass2->first),
                                                                bandwidth.at(chsp).at(sample_mass1->first),
                                                                bandwidth.at(chsp).at(sample_mass2->first));
           }
       }
       return JSDss_range;
      }

     SamplePairNameNDElement JensenDivergenceSSRange(const ChannelSpin& chsp){
        SamplePairNameNDElement JSDss_range;
        for(auto const& range : ranges){
            std::cout<<range.min()<<std::endl;
            if (!range.Contains(args.which_range())) continue;
            for (const auto& sample_mass: samples_mass.at(chsp)){
                if (!range.Contains(sample_mass.first.mass)) continue;
                SampleId rangemin{SampleType::Sgn_Res, range.min()};
                SamplePair mass_pair(rangemin, sample_mass.first);
                JSDss_range[mass_pair] = JensenDivergenceSamples(samples_range.at(chsp).at(SampleId{SampleType::Sgn_Res, range.min()}),
                                                                 samples_mass.at(chsp).at(sample_mass.first),
                                                                 bandwidth_range.at(chsp).at(SampleId{SampleType::Sgn_Res, range.min()}),
                                                                 bandwidth.at(chsp).at(sample_mass.first));
            }
        }
        return JSDss_range;
     }

     void TimeReport(bool tot = false) const
     {
         reporter->TimeReport(tot);
     }

    void Run()
    {
        run::ThreadPull threads(args.number_threads());
        LoadSkimmedData();

        std::string file_name_prefix_JSDSS = "JensenShannonDivergenceSS";
        std::string file_name_prefix_JSDSB = "JensenShannonDivergenceSB";

        if(args.range()) {
            std::cout<<"RANGES"<<std::endl;
            for (const auto& s: set){
                samples_range[s] = LoadRangeData(samples_mass[s]);
            }
            file_name_prefix_JSDSB += "Range";
            file_name_prefix_JSDSS += "Range";
        } else {
            samples_range = samples_mass;
        }

        for (const auto& s: set){
            std::cout<<s.channel<< "  " << s.spin<<std::endl;
            std::stringstream ss;
            ss << std::fixed << std::setprecision(0) << s.spin;
            std::string spin = ss.str();
            for (const auto& sample: samples_mass[s]){
                std::cout<<"----"<<ToString(sample.first)<<"----"<<" entries: "<<sample.second.at("pt_l1").size()<<std::endl;
                std::cout<<args.optband_folder()+"/OptimalBandwidth"+ToString(sample.first)+"_"+s.channel+
                           "_spin"+spin+args.suffix()+".csv"<<std::endl;
                bandwidth[s][sample.first] = Read_csvfile(args.optband_folder()+"/OptimalBandwidth"+ToString(sample.first)+"_"+s.channel+
                                                          "_spin"+spin+args.suffix()+".csv", {});
                std::cout<<bandwidth[s][sample.first].size()<<std::endl;
                for (auto& el : bandwidth[s][sample.first])
                    if (el.second == 0) el.second = 0.0001;
                if(args.range()){
                    bandwidth_range[s][sample.first] = Read_csvfile(args.optband_folder()+"/OptimalBandwidthRange"+ToString(sample.first)+
                                                                    "_"+s.channel+"_spin"+spin+args.suffix()+".csv", {});
                    for (auto& el : bandwidth_range[s][sample.first])
                        if (el.second == 0) el.second = 0.0001;
                }
            }
            bandwidth[s][SampleType::Bkg_TTbar] = Read_csvfile(args.optband_folder()+"/OptimalBandwidthTT_"+s.channel+"_spin"+spin+args.suffix()+".csv", {});
            for (auto& el : bandwidth[s][SampleType::Bkg_TTbar])
                if (el.second == 0) el.second = 0.0001;

        }

        if (args.bkg_vs_sgn()){
            std::cout<<"SIGNAL-BACKGROUND"<<std::endl;
            for (const auto& s: set){
                std::cout<<"JSD"<<s.channel<< "  " << s.spin<<std::endl;
                for (const auto& sample: samples_range[s]){
                    if(args.range())
                        std::cout<<"----Range"<<ToString(sample.first)<<"----"<<" entries: "<<sample.second.begin()->second.size()<<std::endl;
                    else
                        std::cout<<"----"<<ToString(sample.first)<<"----"<<" entries: "<<sample.second.at("pt_l1").size()<<std::endl;

                    std::stringstream ss;
                    ss << std::fixed << std::setprecision(0) << s.spin;
                    std::string spin = ss.str();
                    ChannelSpin chsp_bkg(s.channel,-1);

                    const auto& sgn_band_ptr = args.range() ? bandwidth_range.at(s).at(sample.first) : bandwidth.at(s).at(sample.first);

                    JSDivergenceSB[s][sample.first] = JensenDivergenceSamples(sample.second, samples_mass.at(chsp_bkg).at(SampleType::Bkg_TTbar),
                                                                         sgn_band_ptr, bandwidth.at(chsp_bkg).at(SampleType::Bkg_TTbar));

                    std::ofstream ListJSD(file_name_prefix_JSDSB+ToString(sample.first)+"_"+s.channel+"_spin"+spin+args.suffix()+".csv", std::ofstream::out);
                    std::cout<<"list"<<std::endl;
                    for(const auto& value: JSDivergenceSB[s][sample.first]){
                       for(const auto& var_name: value.first)
                       {
                           ListJSD << var_name << "," ;
                       }
                       ListJSD << value.second << std::endl;
                    }
                    TimeReport();
                }
            }
        }
        else {
            std::cout<<std::endl<<"SIGNAL-SIGNAL"<<std::endl;
            for (const auto& s: set){
                std::cout<<"JSD"<<s.channel<< "  " << s.spin<<std::endl;
                JSDivergenceSS[s] = args.range() ? JensenDivergenceSSRange(s) : JensenDivergenceSS(s);
            }
            for (const auto& entry :  JSDivergenceSS){
                std::cout<<entry.first.channel<< "  " << entry.first.spin<<std::endl;
                for (const auto& sample: entry.second){
                    std::stringstream ss;
                    ss << std::fixed << std::setprecision(0) << entry.first.spin;
                    std::string spin = ss.str();
                    std::ofstream ListJSD(file_name_prefix_JSDSS+ToString(sample.first.first)+"_"+ToString(sample.first.second)+"_"+
                                          entry.first.channel+"_spin"+spin+args.suffix()+".csv", std::ofstream::out);
                    for (const auto& value: sample.second){
                        for (const auto& var: value.first){
                            ListJSD << var << "," ;
                        }
                        ListJSD << value.second << std::endl;
                    }
                }
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
    std::map<ChannelSpin,SampleIdNameElement> bandwidth, bandwidth_range, JSDivergenceSB, JSDivergenceSB_range;
    std::map<ChannelSpin, SamplePairNameNDElement> JSDivergenceSS, JSDivergenceSS_range;
};
}
}

PROGRAM_MAIN(analysis::mva_study::FindJSD, Arguments) // definition of the main program function
