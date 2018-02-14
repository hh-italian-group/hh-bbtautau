/*! Study of correlation matrix, mutual information and Jensen Shannon Divergence to search
  best variables for each range of mass.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <fstream>
#include <algorithm>
#include <functional>
#include <future>
#include <initializer_list>
#include <TCanvas.h>
#include <TPad.h>
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
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(std::string, optband_folder);
    REQ_ARG(std::string, jsd_folder);
    REQ_ARG(std::string, mutual_folder);
    REQ_ARG(std::string, suffix);
    REQ_ARG(unsigned, number_threads);
    REQ_ARG(size_t, number_variables);
    OPT_ARG(Long64_t, number_events, 1000000);
    OPT_ARG(size_t, number_sets, 0);
    OPT_ARG(size_t, set, 0);
    OPT_ARG(bool, is_SM, false);
    OPT_ARG(uint_fast32_t, seed, 10000);
};

namespace analysis {
namespace mva_study{

using clock = std::chrono::system_clock;

class VariableDistribution {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;
    using SummaryTuple = ntuple::SummaryTuple;


    std::vector<ChannelSpin> set_SM{{"tauTau",SM_spin}, {"muTau",SM_spin},{"eTau",SM_spin},
                                 {"muTau",bkg_spin},{"eTau",bkg_spin}, {"tauTau",bkg_spin}};
    std::vector<ChannelSpin> set_R{{"tauTau",0}, {"muTau",0},{"eTau",0},
                                   {"tauTau",2}, {"muTau",2},{"eTau",2},
                                 {"muTau",bkg_spin},{"eTau",bkg_spin}, {"tauTau",bkg_spin}};
    std::vector<ChannelSpin> set;

    std::map<ChannelSpin, SampleIdVarData> samples_mass;
    std::map<ChannelSpin, SampleIdNameElement> bandwidth, mutual_matrix, correlation_matrix, JSDivergenceSB, JSDvars_range_sb;

    std::map<ChannelSpin, SampleIdNameElement> bandwidth_range, mutual_matrix_range, correlation_matrix_range_signal, correlation_matrix_range_bkg;
    std::map<ChannelSpin,SampleIdVarData> samples_range;


    VariableDistribution(const Arguments& _args): args(_args),
        outfile(root_ext::CreateRootFile(args.output_file())), vars(args.number_sets(), args.seed(),{}, {"channel", "mass", "spin", "iso_2",
                                                                    "decayMode_1", "decayMode_2", "iso_1", "csv_1", "csv_2",
                                                                    "costheta_l1l2METhh_sv","costheta_htautau_svhhMET","costheta_htautau_svhh_sv",
                                                                    "costheta_htautau_svhh", "costheta_htautauhh_sv", "costheta_hbbhh_sv",
                                                                    "costheta_METhtautau_sv", "costheta_l2htautau_sv", "costheta_l1htautau_sv",
                                                                    "costheta_star_leptons_sv", "phi_2_sv", "phi_1_sv", "phi_sv", "mass_H_sv",
                                                                    "MT_htautau_sv", "mass_htautau_sv", "dR_l1l2_boosted_sv", "dR_l1l2Pt_htautau_sv",
                                                                    "dR_hbbhtautau_sv", "dR_METhtautau_sv", "deta_hbbhtautau_sv", "abs_deta_hbbhtautau_sv"
                                                                    "deta_METhtautau_sv", "abs_deta_METhtautau_sv","dphi_hbbhtautau_sv", "abs_dphi_hbbhatutau_sv",
                                                                    "dphi_hbbhtautau_sv", "abs_dphi_hbbhatutau_sv", "dphi_METhtautau_sv", "abs_dphi_METhtautau_sv","pt_htautau_sv"
}),
              reporter(std::make_shared<TimeReporter>())
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
        std::cout<<set.size()<<std::endl;

        samples = samples_list.at("Samples").files;
    }

    NameElement CorrelationSelected(const VarData& sample_vars, const SetNamesVar& selected){
        NameElement corr_matrix;
        NameElementFuture corr_matrix_future;
        for(const auto& var_1 : selected) {
            if (!sample_vars.count(var_1)) continue;
            for(const auto& var_2 : selected) {
                if (var_2 < var_1) continue;
                if (!sample_vars.count(var_2)) continue;
                corr_matrix_future[Name_ND{var_1, var_2}] = run::async(stat_estimators::Correlation<double>,
                                                                sample_vars.at(var_1), sample_vars.at(var_2), 1.);
            }
        }
        for(auto& var : corr_matrix_future) {
            corr_matrix[var.first] = var.second.get() *  100;
        }
        return corr_matrix;
    }

    SetNamesVar FindBestRangeVariables(const int& min, const int& max, const std::map<ChannelSpin, std::map<SampleId, double>>& max_distance,
                                      std::map<ChannelSpin,std::map<SampleId, VectorName_ND>>& JSDivergence_vector) const{

        static constexpr double threashold_mi = 0.8;
        SetNamesVar selected, not_corrected;
        std::ofstream best_entries_file("Best_entries_Range"+std::to_string(min)+"_"+std::to_string(max)+"_"+std::to_string(args.set())+"_"+std::to_string(args.number_sets())+args.suffix()+".csv", std::ofstream::out);
        best_entries_file << "," << "," << std::endl;
        const SampleId minSample = args.is_SM() ? SampleId{SampleType::Sgn_NonRes, 0} : SampleId{SampleType::Sgn_Res, min};
        static const std::map<std::pair<size_t, size_t>, std::string> prefixes = {
            { { 0, 0 }, "===" }, { { 0, 1 }, "*" }, { { 1, 0 }, "-" }
        };
        ChannelSpin s = set.front();
        auto& JSDivergence_min = JSDivergence_vector.at(s).at(minSample);
        while(selected.size() < args.number_variables() && JSDivergence_min.size()) {
            std::sort(JSDivergence_min.begin(), JSDivergence_min.end(), [](const auto& el1, const auto& el2){
                return el1.second > el2.second;
            });
            VectorName_ND distance;
            for (const auto& var_entry: JSDivergence_min){
                double d = 0;
                for (const auto& se : set){
                    for (const auto& entry : JSDivergenceSB.at(se)){
                        if (!args.is_SM()){
                            if (entry.first.IsSM() || entry.first.mass<min || entry.first.mass>max) continue;}
                        else
                            if (!entry.first.IsSM()) continue;
                        d+= entry.second.at(var_entry.first) / max_distance.at(se).at(entry.first);
                    }
                }
                distance.emplace_back(var_entry.first, d);
            }   
            std::sort(distance.begin(), distance.end(), [](const auto& el1, const auto& el2){
                return el1.second > el2.second;
            });

            const auto& best_entry = distance.front();
            for(const auto& name : best_entry.first) {
                const auto counts = std::make_pair(selected.count(name), not_corrected.count(name));
                const auto& prefix = prefixes.at(counts);
                best_entries_file << prefix << name << prefix << ",";
                if(counts.first || counts.second)
                    continue;

                for (const auto& se : set){
                    for (const auto& entry : JSDivergenceSB.at(se)){


                        if (!args.is_SM()){
                            if (entry.first.IsSM() || entry.first.mass<min || entry.first.mass>max) continue;}
                        else
                            if (!entry.first.IsSM()) continue;

                        const double JS_1d = entry.second.at(name);
                        for(auto& var : JSDivergence_vector.at(se).at(entry.first)) {
                            if(var.first.count(name))
                                var.second -= JS_1d;
                        }
                    }
                    for (const auto& entry : JSDivergenceSB.at(se)){
                        if (!args.is_SM()){
                            if (entry.first.IsSM() || entry.first.mass<min || entry.first.mass>max) continue;}
                        else
                            if (!entry.first.IsSM()) continue;
                        for(const auto& other_entry : samples_mass.at(se).at(entry.first)) {
                            if(other_entry.first == name) continue;
                            const Name_ND names({name, other_entry.first});
                            ChannelSpin chsp_bkg(se.channel, bkg_spin);
                            if (!mutual_matrix.at(se).at(entry.first).count(names) || !mutual_matrix.at(chsp_bkg).at(SampleType::Bkg_TTbar).count(names)) continue;
                            if(mutual_matrix.at(se).at(entry.first).at(names) < threashold_mi && mutual_matrix.at(chsp_bkg).at(SampleType::Bkg_TTbar).at(names) < threashold_mi){
                                not_corrected.insert(other_entry.first);
                                best_entries_file << "," <<name<<"-"<<other_entry.first << "," <<"MID(sgn) "<<mutual_matrix.at(se).at(entry.first).at(names) << "," <<"MID(bkg) "<<  mutual_matrix.at(chsp_bkg).at(SampleType::Bkg_TTbar).at(names) << "," << entry.first <<",";
                            }
                        }
                    }
                }
                selected.insert(name);
            }
            best_entries_file << std::endl;
            JSDivergence_min = CopySelectedVariables(JSDivergence_min, best_entry.first, not_corrected);
        }
       return selected;
    }

    void CreateMatrixIntersection(const SampleIdSetNamesVar& range_selected) const
    {
        std::cout << "Intersection ranges" << std::endl;
        auto matrix_intersection = std::make_shared<TH2D>("Intersection", "Intersection of variables", mva_study::ranges.size(), 0, mva_study::ranges.size(),
                                                          mva_study::ranges.size(), 0, mva_study::ranges.size());
        matrix_intersection->SetXTitle("range");
        matrix_intersection->SetYTitle("range");
        int c_label = 1;
        int k_row = 1;
        for (const auto& range1: analysis::mva_study::ranges){
            int j_column = 1;
            for (const auto& range2: analysis::mva_study::ranges){
                if (range1.min() == range2.min()) {
                    std::string label = std::to_string(range1.min()) + "-" + std::to_string(range1.max());
                    matrix_intersection->GetXaxis()->SetBinLabel(c_label, (label).c_str());
                    matrix_intersection->GetYaxis()->SetBinLabel(c_label, (label).c_str());
                    c_label++;
                }
                std::vector<std::string> intersection;
                std::set_intersection(range_selected.at(SampleId{SampleType::Sgn_Res, range1.min()}).begin(), range_selected.at(SampleId{SampleType::Sgn_Res, range1.min()}).end(),
                                      range_selected.at(SampleId{SampleType::Sgn_Res, range2.min()}).begin(), range_selected.at(SampleId{SampleType::Sgn_Res, range2.min()}).end(),
                                      std::back_inserter(intersection));
                matrix_intersection->SetBinContent(k_row, j_column, intersection.size());
                j_column++;
            }
            k_row++;
        }
        root_ext::WriteObject(*matrix_intersection, outfile.get());
    }

    SampleIdSetNamesVar VariablesSelection(){
        SampleIdSetNamesVar range_selected;
        std::map<ChannelSpin,std::map<SampleId, VectorName_ND>> JSDivergence_vector;
        std::map<ChannelSpin, std::map<SampleId, double>> max_distance;
        for (const auto& s: set){
            std::cout<<s.channel<<" "<<s.spin<<std::endl;
            for (const auto& entry: JSDivergenceSB.at(s)){
                VectorName_ND  vector(entry.second.begin(), entry.second.end());
                JSDivergence_vector[s][entry.first] = vector;
                std::sort(JSDivergence_vector.at(s).at(entry.first).begin(), JSDivergence_vector.at(s).at(entry.first).end(),
                          [](const auto& el1, const auto& el2){
                    return el1.second > el2.second;
                });
                max_distance[s][entry.first] = JSDivergence_vector.at(s).at(entry.first).front().second;
            }
        }

        if (!args.is_SM()){
            for (const auto& range: mva_study::ranges){
                std::ofstream list_variables(("Selected_range"+std::to_string(range.min())+"_"+std::to_string(range.max())+"_"+std::to_string(args.set())+"_"+std::to_string(args.number_sets())+args.suffix()+".txt").c_str(), std::ofstream::out);
                range_selected[SampleId{SampleType::Sgn_Res, range.min()}] = FindBestRangeVariables(range.min(), range.max(), max_distance, JSDivergence_vector);
                std::cout<<"range: "<<range.min()<<"-"<<range.max()<<std::endl;
                for (const auto& name : range_selected[SampleId{SampleType::Sgn_Res, range.min()}]){
                    list_variables << name << std::endl;
                }
                list_variables << std::endl;
            }
            TimeReport();
            CreateMatrixIntersection(range_selected);
        }
        else{
            std::ofstream list_variables(("Selected_SM_"+std::to_string(args.set())+"_"+std::to_string(args.number_sets())+args.suffix()+".txt").c_str(), std::ofstream::out);
            range_selected[SampleId{SampleType::Sgn_NonRes, 0}] = FindBestRangeVariables(0, 0,  max_distance, JSDivergence_vector);
            for (const auto& name : range_selected[SampleId{SampleType::Sgn_NonRes, 0}]){
                list_variables << name << std::endl;
            }
            list_variables << std::endl;
        }


        return range_selected;
    }

    void KolmogorovSignalPlotSelected(const SetNamesVar& selected, const Range<int>& range, TDirectory* directory, const ChannelSpin& set){
        std::map<std::string, std::shared_ptr<TGraph>> plot;
        int i = 0;
        for (const auto& mass_entry: samples_mass.at(set)){
            if (!range.Contains(mass_entry.first.mass))
                continue;
            for (const auto& var : selected){
                if (!samples_mass.at(set).at(mass_entry.first).count(var)) continue;
                std::vector<double> vector_signal = samples_mass.at(set).at(mass_entry.first).at(var);
                std::sort(vector_signal.begin(), vector_signal.end());
                std::vector<double> vector_signal_2 = samples_mass.at(set).at(SampleId{SampleType::Sgn_Res, range.min()}).at(var);
                std::sort(vector_signal_2.begin(), vector_signal_2.end());
                Double_t* v_s = vector_signal.data(), *v_s_2 = vector_signal_2.data();
                double k = TMath::KolmogorovTest(static_cast<int>(vector_signal.size()), v_s,
                                                 static_cast<int>(vector_signal_2.size()), v_s_2, "");
                if (plot.count(var) == 0) {
                    std::string name = var+"_Signal"+std::to_string(range.min())+"_"+std::to_string(range.max());
                    plot[var] = CreatePlot(name.c_str(), ("Ks_"+name).c_str(), "mass","KS Probability" );
                }
                plot[var]->SetPoint(i, mass_entry.first.mass, k);
            }
            i++;
        }
        for(const auto& var: selected){
            if (!plot.count(var)) continue;
            root_ext::WriteObject(*plot[var], directory);
        }
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
                        if (!cuts::hh_bbtautau_2016::hh_tag::new_m_hh_window().IsInside(eventbase.GetHiggsTTMomentum(false).M(),bb.mass())) continue;
                    }
                    vars.AddEvent(eventbase, entry.id, entry.spin, entry.weight);
                    tot_entries++;
                }
                std::cout << entry << " number of events: " << tot_entries << std::endl;
            }
            samples_mass[s] = vars.GetSampleVariables(s.channel, s.spin);
            TimeReport();
        }
    }

    void TimeReport(bool tot = false) const
    {
        reporter->TimeReport(tot);
    }

    void Run()
    {


        run::ThreadPull threads(args.number_threads());
        LoadSkimmedData();

        std::cout << "Bandwidth, Mutual Information, JensenShannon" << std::endl;
        for (const auto& s: set){
            std::cout<<s.channel<< "  " << s.spin<<std::endl;
            for (const auto& sample: samples_mass[s]){
                std::cout<<"----"<<ToString(sample.first)<<"----"<<" entries: "<<sample.second.at("pt_l1").size()<<std::endl;
                std::stringstream ss;
                ss << std::fixed << std::setprecision(0) << s.spin;
                std::string spin = ss.str();
                bandwidth[s][sample.first] = Read_csvfile(args.optband_folder()+"/OptimalBandwidth"+ToString(sample.first)+"_"+s.channel+"_spin"+spin+args.suffix()+".csv", vars.GetDisabledVars());
                mutual_matrix[s][sample.first] = Read_csvfile(args.mutual_folder()+"/MutualInformationDistance"+ToString(sample.first)+"_"+s.channel+"_spin"+spin+args.suffix()+".csv", vars.GetDisabledVars());
                JSDivergenceSB[s][sample.first] = Read_csvfile(args.jsd_folder()+"/JensenShannonDivergenceSB"+ToString(sample.first)+"_"+s.channel+"_spin"+spin+args.suffix()+".csv", vars.GetDisabledVars());
            }
        }

        std::map<ChannelSpin, std::map<Name_ND, std::shared_ptr<TGraph>>> plot_jsd;

        std::cout<<"Jensen Shannon Signal Background "<<std::endl;
        for (const auto& s: set){
            std::string spin;
            if (s.spin == 0) spin = "Radion";
            else if (s.spin == 2) spin = "Graviton";
            else if (s.spin == SM_spin) spin = "SM";
            else spin = "Bkg";
            std::cout << s.channel << "  " << spin << std::endl;
            int i = 0;
            for (const auto& sample: samples_mass.at(s)){
                if ( sample.first.IsBackground() ){
                    continue;
                }
                for (const auto& var : JSDivergenceSB.at(s).at(sample.first)){
                    if (!plot_jsd[s].count(var.first)) plot_jsd[s][var.first] = CreatePlot("","","","");
                    plot_jsd.at(s).at(var.first)->SetPoint(i, sample.first.mass, JSDivergenceSB.at(s).at(sample.first).at(var.first));
                }
                i++;
            }
            TimeReport();
        }


        std::cout<<"Selection variables"<<std::endl;
        SampleIdSetNamesVar range_selected = VariablesSelection();

        std::map<ChannelSpin, std::map<int, std::map<Name_ND, std::shared_ptr<TGraph>>>> plot_ss;
        for (const auto& s: set){
            std::string spin;
            if (s.spin == 0) spin = "Radion";
            else if (s.spin == 2) spin = "Graviton";
            else if (s.spin == SM_spin) spin = "SM";
            else spin = "Bkg";
            std::cout << s.channel << "  " << spin << std::endl;

            if (args.is_SM()) continue;

            auto directory_set = root_ext::GetDirectory(*outfile, s.channel+spin);
            auto directory_jensenshannon = root_ext::GetDirectory(*directory_set,"JensenShannonDivergence");
            std::cout<<"Signal Compatibility"<<std::endl;
            std::map<std::pair<int,int>, NameElement> JSDvars_range_ss;
            std::map<std::pair<int,int>, std::shared_ptr<TH1D>> histo_distribution;
            auto directory_ks = root_ext::GetDirectory(*directory_set,"Kolmogorov");
            auto directory_jenshan_ss = root_ext::GetDirectory(*directory_jensenshannon, "Range_SignalSignal");
            auto directory_plotselected = root_ext::GetDirectory(*directory_jenshan_ss, "Plot_RangeSelected");

            for (const auto& range : mva_study::ranges){
                KolmogorovSignalPlotSelected(range_selected.at(SampleId{SampleType::Sgn_Res, range.min()}), range, directory_ks, s);
                for (const auto& sample_mass: samples_mass.at(s)){
                    if (!range.Contains(sample_mass.first.mass)) continue;
                    std::pair<int, int> mass_pair(range.min(), sample_mass.first.mass);
                    std::stringstream ss;
                    ss << std::fixed << std::setprecision(0) << s.spin;
                    std::string spin = ss.str();
                    JSDvars_range_ss[mass_pair] = Read_csvfile(args.jsd_folder()+"/JensenShannonDivergenceSSRangeM"+ToString(range.min())+"_"+ToString(sample_mass.first)+"_"+s.channel+"_spin"+spin+args.suffix()+".csv", vars.GetDisabledVars());
                    histo_distribution[mass_pair] = std::make_shared<TH1D>(("JSDrange_"+std::to_string(range.min())+"_"+std::to_string(sample_mass.first.mass)).c_str(), ("JSDrange_"+std::to_string(range.min())+"_"+std::to_string(sample_mass.first.mass)).c_str(), 50,0,1);
                    histo_distribution[mass_pair]->SetXTitle("JSD");

                    for (const auto& value: JSDvars_range_ss[mass_pair]){
                        Name_ND var_pair{};
                        std::string name;
                        for(const auto& var: value.first){
                            var_pair.insert(var);
                            name=var+"_";
                        }
                        if (!plot_ss[s][range.min()].count(var_pair)) {
                            plot_ss[s][range.min()][var_pair] = CreatePlot(("JSD_"+name+"Range"+std::to_string(range.min())+"_"+std::to_string(range.max())).c_str(),("JSD_"+name+"Range"+std::to_string(range.min())+"_"+std::to_string(range.max())).c_str(),"mass","JSD");
                        }
                    }
                }
            }

            if (s.spin == bkg_spin) continue;
            int k = 0;
            for(auto& mass_pair : JSDvars_range_ss) {
                if (mass_pair.first.first == mass_pair.first.second) k = 0;
                for(auto& name : mass_pair.second) {
                    histo_distribution.at(mass_pair.first)->Fill(JSDvars_range_ss.at(mass_pair.first).at(name.first));
                    plot_ss.at(s).at(mass_pair.first.first).at(name.first)->SetPoint(k, mass_pair.first.second, JSDvars_range_ss.at(mass_pair.first).at(name.first));
                }
                k++;
            }
            for (const auto& range : mva_study::ranges){
                auto directory_plotrange = root_ext::GetDirectory(*directory_plotselected, ("Range"+std::to_string(range.min())+"_"+std::to_string(range.max())).c_str());
                std::cout<<range.min()<<std::endl;
                for (const auto& name : plot_ss.at(s).at(range.min())){
                    root_ext::WriteObject(*name.second, directory_plotrange);
                }
            }
            auto directory_distributionselected = root_ext::GetDirectory(*directory_jenshan_ss, "Distribution_RangeSelected");
            for (const auto& mass_pair : histo_distribution){
                root_ext::WriteObject(*mass_pair.second, directory_distributionselected);
            }
            TimeReport();
        }


        std::cout<<"Union sample"<<std::endl;
        for (const auto& s : set){
            std::cout<< s.channel << "    " << s.spin<< std::endl;
            std::string spin;
            if (s.spin == 0) spin = "Radion";
            else if (s.spin == 2) spin = "Graviton";
            else continue;
            std::cout<<"unire"<<std::endl;
            for (const auto& range : mva_study::ranges){
                for (const auto& var: range_selected.at(SampleId{SampleType::Sgn_Res, range.min()})){
                    for (const auto& sample_mass : samples_mass.at(s)){
                        if (!range.Contains(sample_mass.first.mass)) continue;
                        if (!sample_mass.second.count(var)) continue;
                        for (const auto& entry : sample_mass.second.at(var)){
                            samples_range[s][SampleId{SampleType::Sgn_Res, range.min()}][var].push_back(entry);
                        }
                    }
                }
            }

            std::cout<<"band mut JSD"<<std::endl;
            for (const auto& sample: samples_range.at(s)){
                std::stringstream ss;
                ss << std::fixed << std::setprecision(0) << s.spin;
                std::string spin = ss.str();
                bandwidth_range[s][sample.first] = Read_csvfile(args.optband_folder()+"/OptimalBandwidthRange"+ToString(sample.first)+"_"+s.channel+"_spin"+spin+args.suffix()+".csv", vars.GetDisabledVars());
                mutual_matrix_range[s][sample.first] = Read_csvfile(args.mutual_folder()+"/MutualInformationDistanceRange"+ToString(sample.first)+"_"+s.channel+"_spin"+spin+args.suffix()+".csv", vars.GetDisabledVars());
                JSDvars_range_sb[s][sample.first] = Read_csvfile(args.jsd_folder()+"/JensenShannonDivergenceSBRange"+ToString(sample.first)+"_"+s.channel+"_spin"+spin+args.suffix()+".csv", vars.GetDisabledVars());
            }

            ChannelSpin chsp_bkg(s.channel, bkg_spin);
            auto directory_set = root_ext::GetDirectory(*outfile, s.channel+spin);
            auto directory_mutinf = root_ext::GetDirectory(*directory_set, "MutualInformation");
            auto directory_correlation = root_ext::GetDirectory(*directory_set,"Correlation");

            for (const auto& sample: samples_range.at(s)){
                std::cout<<"----Range"<<ToString(sample.first)<<"----"<<" entries: "<<sample.second.begin()->second.size()<<std::endl;
                std::cout<<"correlation  "<<std::endl;
                correlation_matrix_range_signal[s][sample.first] = CorrelationSelected(sample.second, range_selected.at(sample.first));
                correlation_matrix_range_bkg[chsp_bkg][sample.first] = CorrelationSelected(samples_mass.at(chsp_bkg).at(SampleType::Bkg_TTbar), range_selected.at(sample.first));

                int bin = static_cast<int>(range_selected.at(SampleId{SampleType::Sgn_Res, sample.first.mass}).size());
                auto matrix = std::make_shared<TH2D>(("Bkg_"+std::to_string(sample.first.mass)).c_str(),("Bkg_"+std::to_string(sample.first.mass)).c_str(),bin,0,bin,bin,0,bin);
                int i = 1;

                for(const auto& var_1 : range_selected.at(sample.first)) {
                    int j = 1;
                    if (!correlation_matrix_range_bkg.at(chsp_bkg).at(sample.first).count(Name_ND{var_1, var_1})) continue;
                    matrix->GetXaxis()->SetBinLabel(i, (var_1).c_str());
                    matrix->GetYaxis()->SetBinLabel(i, (var_1).c_str());
                    for(const auto& var_2 : range_selected.at(sample.first)) {
                        if (!correlation_matrix_range_bkg.at(chsp_bkg).at(sample.first).count(Name_ND{var_1, var_2})) continue;
                        matrix->SetBinContent(i, j, correlation_matrix_range_bkg.at(chsp_bkg).at(sample.first).at(Name_ND{var_1, var_2}));
                        j++;
                    }
                    i++;
                }
                root_ext::WriteObject(*matrix, directory_correlation);

                if ( sample.first.sampleType == SampleType::Bkg_TTbar){
                    continue;
                }
                else {
                    std::cout<<"mutual plot  ";
                    MutualHisto(sample.first, mutual_matrix.at(s).at(sample.first), mutual_matrix.at(chsp_bkg).at(SampleType::Bkg_TTbar), directory_mutinf);
                }
                TimeReport();
            }
            std::cout <<"Mutual histos" << std::endl;
            directory_mutinf->mkdir("Matrix");
            auto directory_mutinf_matrix = root_ext::GetDirectory(*directory_mutinf, "Matrix");
            CreateMatrixHistos(samples_range.at(s), mutual_matrix_range.at(s), "MI", directory_mutinf_matrix);
            std::cout <<"Correlation histos" << std::endl;
            CreateMatrixHistos(samples_range.at(s), correlation_matrix_range_signal.at(s), "Signal", directory_correlation, true);
            TimeReport();

            std::cout<<"Jensen Shannon union Signal Backgound"<<std::endl;

            std::map<Name_ND, std::shared_ptr<TGraph>> plot_sb;
            std::map<SampleId, std::shared_ptr<TH1D>> histo_distribution_sb;

            for (const auto& range : mva_study::ranges){
                SampleId samplerangemin(SampleType::Sgn_Res, range.min());
                histo_distribution_sb[samplerangemin] = std::make_shared<TH1D>(("JSDrange_"+std::to_string(range.min())+"_SignalBkg").c_str(),
                                                                          ("JSDrange_"+std::to_string(range.min())+"_SignalBkg").c_str(),
                                                                          50,0,1);
                histo_distribution_sb[samplerangemin]->SetXTitle("JSD");
                for (const auto& value: JSDvars_range_sb.at(s).at(samplerangemin)){
                    if(value.first.size()==1)
                    {
                        auto var = value.first.begin();
                        if (!plot_sb.count(*var)){
                            plot_sb[*var] = CreatePlot(("JSD_"+*var+"_SignalBkg").c_str(), ("JSD_"+*var+"_SignalBkg").c_str(), "mass", "JSD");
                        }
                    }
                    else{
                        Name_ND var_pair{};
                        std::string name;
                        for(const auto& var: value.first){
                            var_pair.insert(var);
                            name=var+"_";
                        }
                        if (!plot_sb.count(var_pair)) {
                            plot_sb[var_pair] = CreatePlot(("JSD_"+name+"_SignalBkg").c_str(), ("JSD_"+name+"_SignalBkg").c_str(), "mass", "JSD");
                        }
                    }
                }
            }
            auto directory_jensenshannon = root_ext::GetDirectory(*directory_set,"JensenShannonDivergence");
            auto directory_rangebkg = root_ext::GetDirectory(*directory_jensenshannon, "Range_SignalBkg");
            auto directory_plotselected_sb = root_ext::GetDirectory(*directory_rangebkg, "Plot_RangeSelected");
            auto directory_distributionselected_sb = root_ext::GetDirectory(*directory_rangebkg, "Distribution_RangeSelected");
            int w = 0;

            for(auto& mass : JSDvars_range_sb.at(s)) {
                for(auto& name : mass.second) {
                    histo_distribution_sb.at(mass.first)->Fill(JSDvars_range_sb.at(s).at(mass.first).at(name.first));
                    plot_sb.at(name.first)->SetPoint(w, mass.first.mass, JSDvars_range_sb.at(s).at(mass.first).at(name.first));
                }
                w++;
            }
            for (const auto& name : plot_sb){
                root_ext::WriteObject(*name.second, directory_plotselected_sb);
            }
            for (const auto& mass_pair : histo_distribution_sb){
                root_ext::WriteObject(*mass_pair.second, directory_distributionselected_sb);
            }
            TimeReport();


            std::cout<<"Plot matrix"<<std::endl;
            auto directory_matrix = root_ext::GetDirectory(*directory_rangebkg, "Matrix");
            for (const auto& range: mva_study::ranges){
                SampleId samplerangemin(SampleType::Sgn_Res, range.min());
                int bin = static_cast<int>(range_selected.at(SampleId{SampleType::Sgn_Res, range.min()}).size());
                auto matrix_jsd_sb = std::make_shared<TH2D>(("JSD_Signal_Bkg_Range"+std::to_string(range.min())+"_"+std::to_string(range.min())).c_str(),("JSD_Signal_Bkg_Range"+std::to_string(range.min())+"_"+std::to_string(range.min())).c_str(), bin, 0, bin, bin, 0, bin);
                int i = 1;
                for (const auto& var1: range_selected.at(SampleId{SampleType::Sgn_Res, range.min()})){
                        matrix_jsd_sb->GetXaxis()->SetBinLabel(i, (var1).c_str());
                        matrix_jsd_sb->GetYaxis()->SetBinLabel(i, (var1).c_str());
                        int jj = 1;
                        for (const auto& var2: range_selected.at(SampleId{SampleType::Sgn_Res, range.min()})){
                            if (var1 == var2){
                                matrix_jsd_sb->SetBinContent(i, i, JSDvars_range_sb.at(s).at(samplerangemin).at(var1));
                            }
                            if (jj > i){
                                matrix_jsd_sb->SetBinContent(i, jj, JSDvars_range_sb.at(s).at(samplerangemin).at(Name_ND{var1, var2}));
                                matrix_jsd_sb->SetBinContent(jj, i, JSDvars_range_sb.at(s).at(samplerangemin).at(Name_ND{var1, var2}));
                            }
                            jj++;
                        }
                        i++;
                }
                root_ext::WriteObject(*matrix_jsd_sb, directory_matrix);
            }
            auto directory_compare = root_ext::GetDirectory(*directory_jensenshannon, "Comparison");
            std::map<Name_ND, TCanvas*> canvas;
            for (const auto & range: mva_study::ranges){
                for (const auto& selected_name : plot_ss.at(s).at(range.min())){
                    if (canvas.count(selected_name.first)) continue;
                    canvas[selected_name.first] = new TCanvas(selected_name.first.ToString().c_str(), selected_name.first.ToString().c_str(), 10,10,800,800);
                    canvas[selected_name.first]->DrawFrame(200,0,1000,1);
                    canvas[selected_name.first]->SetGrid();
                }
            }
            for (const auto& selected_name : plot_sb){
                canvas.at(selected_name.first)->cd();
                selected_name.second->Draw("P");
                for (const auto & range: mva_study::ranges){
                    if (!plot_ss.at(s).at(range.min()).count(selected_name.first)) continue;
                    plot_ss.at(s).at(range.min()).at(selected_name.first)->Draw();
                }
                plot_jsd.at(s).at(selected_name.first)->Draw();
                root_ext::WriteObject<TCanvas>(*canvas[selected_name.first], directory_compare);
                canvas[selected_name.first]->Close();
            }
        }
        TimeReport(true);
    }

private:
    Arguments args;
    SampleEntryCollection samples;
    std::shared_ptr<TFile> outfile;
    MvaVariablesStudy vars;
    std::shared_ptr<TimeReporter> reporter;
};
}
}

PROGRAM_MAIN(analysis::mva_study::VariableDistribution, Arguments) // definition of the main program function
