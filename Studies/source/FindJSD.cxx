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
    std::map<ChannelSpin,SampleIdNameElement> bandwidth, bandwidth_range, JSDivergenceSB, JSDivergenceSB_range;
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
                std::cout << entry << " number of events: " << tuple->size() << std::endl;
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

    SamplePairNameNDElement JensenDivergenceSS(const ChannelSpin& chsp){
       SamplePairNameNDElement JSDss_range;
       std::map<SamplePair, std::map<Name_ND, std::future<double>>> JSDss_range_future;
       for(auto sample_mass1 = samples_mass.at(chsp).begin(); sample_mass1 != samples_mass.at(chsp).end(); ++sample_mass1) {
           if ( sample_mass1->first.IsBackground()) continue;
           for(auto sample_mass2 = sample_mass1; sample_mass2 != samples_mass.at(chsp).end(); ++sample_mass2) {
               if ( sample_mass2->first.IsBackground()) continue;
               SamplePair mass_pair(sample_mass1->first, sample_mass2->first);
               ChannelSpin chsp_bkg(chsp.first,-1);
               for (const auto& var: samples_mass.at(chsp_bkg).at(SampleType::Bkg_TTbar)){
                   Name_ND var_1{var.first};
                   std::vector<const DataVector*> sample_1, sample_2;
                   DataVector band_1, band_2;
                   sample_1.push_back(&samples_mass.at(chsp).at(sample_mass1->first).at(var.first));
                   band_1.push_back(bandwidth.at(chsp).at(sample_mass1->first).at(var_1));
                   sample_2.push_back(&samples_mass.at(chsp).at(sample_mass2->first).at(var.first));
                   band_2.push_back( bandwidth.at(chsp).at(sample_mass2->first).at(var_1));
                   JSDss_range_future[mass_pair][var_1]  = run::async(stat_estimators::JensenShannonDivergence_ND<double>,
                                                                         sample_1, sample_2, band_1, band_2);
                   for (const auto& other_var: samples_mass.at(chsp_bkg).at(SampleType::Bkg_TTbar)){
                       if (other_var.first <= var.first) continue;
                       Name_ND var2{other_var.first};
                       Name_ND var_pair{var.first, other_var.first};
                       sample_1.push_back(&samples_mass.at(chsp).at(sample_mass1->first).at(other_var.first));
                       band_1.push_back(bandwidth.at(chsp).at(sample_mass1->first).at(var2));
                       sample_2.push_back(&samples_mass.at(chsp).at(sample_mass2->first).at(other_var.first));
                       band_2.push_back( bandwidth.at(chsp).at(sample_mass2->first).at(var2));
                       JSDss_range_future[mass_pair][var_pair] = run::async(stat_estimators::JensenShannonDivergence_ND<double>,
                                                                                           sample_1, sample_2, band_1, band_2);

                       sample_1.erase(sample_1.end() - 1);
                       sample_2.erase(sample_2.end() - 1);
                       band_1.erase(band_1.end() - 1);
                       band_2.erase(band_2.end() - 1);
                   }
               }
           }
           TimeReport();
       }
       for(auto& mass_pair : JSDss_range_future) {
           std::cout<<ToString(mass_pair.first.first)<<"   "<<ToString(mass_pair.first.second)<<std::endl;
           for(auto& name : mass_pair.second) {
               JSDss_range[mass_pair.first][name.first] = name.second.get();
           }
       }
       return JSDss_range;
      }


     SamplePairNameNDElement JensenDivergenceSSRange(const ChannelSpin& chsp){
        SamplePairNameNDElement JSDss_range;
        std::map<SamplePair, std::map<Name_ND, std::future<double>>> JSDss_range_future;
        for(auto const& range : ranges){
            std::cout<<range.min()<<std::endl;
            for (const auto& sample_mass: samples_mass.at(chsp)){
                if (!range.Contains(sample_mass.first.mass)) continue;
                SampleId rangemin{SampleType::Sgn_Res, range.min()};
                SamplePair mass_pair(rangemin, sample_mass.first);
                ChannelSpin chsp_bkg(chsp.first,-1);
                for (const auto& var: samples_mass.at(chsp_bkg).at(SampleType::Bkg_TTbar)){
                    Name_ND var_1{var.first};
                    std::vector<const DataVector*> sample_1, sample_2;
                    DataVector band_1, band_2;
                    sample_1.push_back(&samples_range.at(chsp).at(SampleId{SampleType::Sgn_Res, range.min()}).at(var.first));
                    band_1.push_back(bandwidth_range.at(chsp).at(SampleId{SampleType::Sgn_Res, range.min()}).at(var_1));
                    sample_2.push_back(&samples_mass.at(chsp).at(sample_mass.first).at(var.first));
                    band_2.push_back( bandwidth.at(chsp).at(sample_mass.first).at(var_1));
                    JSDss_range_future[mass_pair][var_1]  = run::async(stat_estimators::JensenShannonDivergence_ND<double>,
                                                                          sample_1, sample_2, band_1, band_2);
                    for (const auto& other_var: samples_mass.at(chsp_bkg).at(SampleType::Bkg_TTbar)){
                        if (other_var.first <= var.first) continue;
                        Name_ND var2{other_var.first};
                        Name_ND var_pair{var.first, other_var.first};
                        sample_1.push_back(&samples_range.at(chsp).at(SampleId{SampleType::Sgn_Res, range.min()}).at(other_var.first));
                        band_1.push_back(bandwidth_range.at(chsp).at(SampleId{SampleType::Sgn_Res, range.min()}).at(var2));
                        sample_2.push_back(&samples_mass.at(chsp).at(sample_mass.first).at(other_var.first));
                        band_2.push_back( bandwidth.at(chsp).at(sample_mass.first).at(var2));
                        JSDss_range_future[mass_pair][var_pair] = run::async(stat_estimators::JensenShannonDivergence_ND<double>,
                                                                                            sample_1, sample_2, band_1, band_2);

                        sample_1.erase(sample_1.end() - 1);
                        sample_2.erase(sample_2.end() - 1);
                        band_1.erase(band_1.end() - 1);
                        band_2.erase(band_2.end() - 1);
                    }
                }
            }
            TimeReport();
        }
        for(auto& mass_pair : JSDss_range_future) {
            for(auto& name : mass_pair.second) {
                JSDss_range[mass_pair.first][name.first] = name.second.get();
            }
        }
        return JSDss_range;
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
                }
            }

            std::cout<<"SIGNAL-BACKGROUND"<<std::endl;
            for (const auto& s: set){
                std::cout<<std::endl<<s.first<< "  " << s.second<<std::endl;
                for (const auto& sample: samples_mass[s]){
                    std::cout<<"----"<<ToString(sample.first)<<"----"<<" entries: "<<sample.second.at("pt_l1").size()<<std::endl;
                    std::stringstream ss;
                    ss << std::fixed << std::setprecision(0) << s.second;
                    std::string spin = ss.str();
                    std::pair<std::string,double> chsp_bkg(s.first,-1);
                    JSDivergenceSB[s][sample.first] = JensenDivergenceSB(sample.second, samples_mass.at(chsp_bkg).at(SampleType::Bkg_TTbar), bandwidth.at(s).at(sample.first), bandwidth.at(chsp_bkg).at(SampleType::Bkg_TTbar));
                    std::ofstream ListJSD("JensenShannonDivergenceSB"+ToString(sample.first)+"_"+s.first+"_spin"+spin+".csv", std::ofstream::out);
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
                    std::pair<std::string,double> chsp_bkg(s.first.first,-1);
                    JSDivergenceSB_range[s.first][sample.first] = JensenDivergenceSB(sample.second, samples_mass.at(chsp_bkg).at(SampleType::Bkg_TTbar), bandwidth_range.at(s.first).at(sample.first), bandwidth.at(chsp_bkg).at(SampleType::Bkg_TTbar));
                    std::ofstream ListJSD("JensenShannonDivergenceSBRange"+ToString(sample.first)+"_"+s.first.first+"_spin"+spin+".csv", std::ofstream::out);
                    for(const auto& value: JSDivergenceSB_range[s.first][sample.first]){
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        std::cout<<std::endl<<"SIGNAL-SIGNAL"<<std::endl;

        if (!args.range()){
            for (const auto& s: set){
                JSDivergenceSS[s] = JensenDivergenceSS(s);
            }

            for (const auto& s: set){
                for(auto sample_mass1 = samples_mass.at(s).begin(); sample_mass1 != samples_mass.at(s).end(); ++sample_mass1) {
                     if ( sample_mass1->first.IsBackground()) continue;
                     for(auto sample_mass2 = sample_mass1; sample_mass2 != samples_mass.at(s).end(); ++sample_mass2) {
                         if ( sample_mass2->first.IsBackground()) continue;
                        SamplePair mass_pair(sample_mass1->first, sample_mass2->first);
                        std::stringstream ss;
                        ss << std::fixed << std::setprecision(0) << s.second;
                        std::string spin = ss.str();
                        std::ofstream ListJSD("JensenShannonDivergenceSS"+ToString(sample_mass1->first)+"_"+ToString(sample_mass2->first)+"_"+s.first+"_spin"+spin+".csv", std::ofstream::out);
                        for (const auto& value: JSDivergenceSS[s][mass_pair]){
                            for (const auto& var: value.first){
                                ListJSD << var << "," ;
                            }
                            ListJSD << value.second << std::endl;
                        }
                    }
                }
            }
        }
        else{
            std::cout<<"RANGES"<<std::endl;
            for (const auto& s: set){
                std::cout<<std::endl<<s.first<< "  " << s.second<<std::endl;
                JSDivergenceSS_range[s] = JensenDivergenceSSRange(s);
            }

            for (const auto& s: JSDivergenceSS_range){
                std::cout<<std::endl<<s.first.first<< "  " << s.first.second<<std::endl;
                for (const auto& sample: s.second){
                    std::cout<<"----Range"<<ToString(sample.first.first)<<"----"<<std::endl;
                    if (sample.second.size()<2) continue;
                    SamplePair mass_pair(sample.first.first, sample.first.second);
                    std::stringstream ss;
                    ss << std::fixed << std::setprecision(0) << s.first.second;
                    std::string spin = ss.str();
                    std::ofstream ListJSD("JensenShannonDivergenceSSRange"+ToString(sample.first.first)+"_"+ToString(sample.first.second)+"_"+s.first.first+"_spin"+spin+".csv", std::ofstream::out);
                    for (const auto& value: JSDivergenceSS_range[s.first][mass_pair]){
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
