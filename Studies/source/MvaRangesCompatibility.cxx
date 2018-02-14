/*! Study of correlation matrix, mutual information and Jensen Shannon Divergence to search ranges of mass.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <fstream>
#include <functional>
#include <future>
#include <initializer_list>
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/StatEstimators.h"
#include "hh-bbtautau/Analysis/include/MvaVariables.h"
#include "AnalysisTools/Core/include/ProgressReporter.h"
#include "AnalysisTools/Run/include/MultiThread.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "hh-bbtautau/Studies/include/MvaMethods.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
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
    OPT_ARG(uint_fast32_t, seed, 10000);
};

namespace analysis {
namespace mva_study{

using SampleEntryCollection = std::vector<SampleEntry>;

class MvaClassification {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;
    using VecVariables =  std::vector<std::string>;

    std::map<ChannelSpin,SampleIdVarData> samples_mass;
    std::map<ChannelSpin,SampleIdNameElement> bandwidth, mutual_matrix, correlation_matrix, JSDivergenceSB;

    std::vector<ChannelSpin> set{{"muTau",0}, {"eTau",0}, {"tauTau",0},
                                 {"muTau",2}, {"eTau",2}, {"tauTau",2},
                                 {"muTau",SM_spin}, {"eTau",SM_spin}, {"tauTau",SM_spin},
                                 {"muTau",bkg_spin}, {"eTau",bkg_spin}, {"tauTau",bkg_spin}};

    MvaClassification(const Arguments& _args): args(_args),
        outfile(root_ext::CreateRootFile(args.output_file())), vars(args.number_sets(), args.seed(),{}, {"channel", "mass", "spin",
//                                                                    "decayMode_1", "decayMode_2", "iso_1", "iso_2", "csv_1", "csv_2",
                                                                    "costheta_l1l2METhh_sv","costheta_htautau_svhhMET","costheta_htautau_svhh_sv",
                                                                    "costheta_htautau_svhh", "costheta_htautauhh_sv", "costheta_hbbhh_sv",
                                                                    "costheta_METhtautau_sv", "costheta_l2htautau_sv", "costheta_l1htautau_sv",
                                                                    "costheta_star_leptons_sv", "phi_2_sv", "phi_1_sv", "phi_sv", "mass_H_sv",
                                                                    "MT_htautau_sv", "mass_htautau_sv", "dR_l1l2_boosted_sv", "dR_l1l2Pt_htautau_sv",
                                                                    "dR_hbbhtautau_sv", "dR_METhtautau_sv", "deta_hbbhtautau_sv", "abs_deta_hbbhtautau_sv"
                                                                    "deta_METhtautau_sv", "abs_deta_METhtautau_sv","dphi_hbbhtautau_sv", "abs_dphi_hbbhatutau_sv",
                                                                    "dphi_hbbhtautau_sv", "abs_dphi_hbbhatutau_sv", "dphi_METhtautau_sv", "abs_dphi_METhtautau_sv"
                                                                    "pt_htautau_sv"}),
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

        samples = samples_list.at("Samples").files;
    }

    void DistributionJSD_SB(const SampleId& mass_entry, TDirectory* directory, const ChannelSpin& set) const
    {
        auto mass = ToString(mass_entry);
        auto histo = std::make_shared<TH1D>(("JSD_"+mass+"_Background").c_str(), ("JensenShannonDivergence_Signal"+mass+"_Background").c_str(), 50,0,1);
        histo->SetXTitle("JSD");
        std::string spin;
        if (set.spin == 0) spin = "Radion";
        else if (set.spin == 2) spin = "Graviton";
        else if (set.spin == SM_spin) spin = "SM";
        else spin = "Bkg";
        std::ofstream of("InformationTable_"+ToString(mass)+"_"+std::to_string(args.set())+"_"+std::to_string(args.number_sets())
                         +set.channel+"_spin"+spin+args.suffix()+".csv", std::ofstream::out);
        of<<"Var_1"<<","<<"Var_2"<<","<<"JSD_ND"<<","<<"JSD_12-(JSD_1+JSD_2)"<<","<<"ScaledMI_Signal_12" <<","<<"ScaledMI_Bkg_12"<<","<<"selected"<<","<<","<<"eliminated by"<<std::endl;
        for(auto& entry : JSDivergenceSB.at(set).at(mass_entry)){
             histo->Fill(entry.second);
             root_ext::WriteObject(*histo, directory);
        }
    }

    void PlotJensenShannonSB(TDirectory* directory, const ChannelSpin& set) const
    {
        std::map<std::string, std::shared_ptr<TGraph>> plot;
        int i = 0;
        for (const auto& entry: JSDivergenceSB.at(set)){
            if ( entry.first.sampleType == SampleType::Bkg_TTbar )
                continue;
            for (const auto& var : entry.second){
                if (var.first.size()!=1)
                    continue;
                const std::string name = *var.first.begin();
                if (!plot.count(name))
                    plot[name] = CreatePlot("JSD_"+name+"_SignalBkg","JSD_"+name+"_SignalBkg","mass","JSD");
                plot[name]->SetPoint(i, entry.first.mass, var.second);
            }
            i++;
        }
        for(const auto& var: plot){
            root_ext::WriteObject(*var.second, directory);
        }
    }

    void JensenShannonSignalCompatibility(TDirectory* directory, const ChannelSpin& set) const
    {
        SamplePairNameNDElement JSDivergenceSS;
        std::stringstream ss;
        ss << std::fixed << std::setprecision(0) << set.spin;
        std::string spin = ss.str();

        for(auto mass_entry_1 = samples_mass.at(set).begin(); mass_entry_1 != samples_mass.at(set).end(); ++mass_entry_1) {
            if ( mass_entry_1->first.IsBackground())
                continue;
            for(auto mass_entry_2 = mass_entry_1; mass_entry_2 != samples_mass.at(set).end(); ++mass_entry_2) {
                if ( mass_entry_2->first.IsBackground())
                    continue;
                SamplePair mass_pair(mass_entry_1->first,mass_entry_2->first);
                JSDivergenceSS[mass_pair] = Read_csvfile(args.jsd_folder()+"/JensenShannonDivergenceSS"+ToString(mass_entry_1->first)
                                                         +"_"+ToString(mass_entry_2->first)+"_"+set.channel+"_spin"+spin+args.suffix()+".csv", vars.GetDisabledVars());
            }
        }
        std::map<Name_ND, std::map<SamplePair,double>> histovalue;
        for(const auto& entrypair: JSDivergenceSS){
            SamplePair masspair = entrypair.first;
            for (const auto& entryelement: entrypair.second){
                histovalue[entryelement.first][masspair] = entryelement.second;
            }
        }
        for(const auto& var: histovalue){
            CreateMatrixHistosCompatibility(samples_mass.at(set), var.second, ToString(var.first)+"_Signal_compatibility_JSD", directory);
        }
    }

    void KolmogorovSignalCompatibility(TDirectory* directory, const ChannelSpin& set) const
    {
        ChannelSpin chsp_bkg(set.channel, bkg_spin);
        for (const auto& var: samples_mass.at(chsp_bkg).at(SampleType::Bkg_TTbar)){
            std::map<SamplePair, double> kolmogorov;
            for(auto mass_entry_1 = samples_mass.at(set).begin(); mass_entry_1 != samples_mass.at(set).end(); ++mass_entry_1) {
                if ( mass_entry_1->first.IsBackground())
                    continue;
                for(auto mass_entry_2 = mass_entry_1; mass_entry_2 != samples_mass.at(set).end(); ++mass_entry_2) {
                    if ( mass_entry_2->first.IsBackground())
                        continue;
                    SamplePair mass_pair(mass_entry_1->first, mass_entry_2->first);
                    std::vector<double> vector_signal = samples_mass.at(set).at(mass_entry_1->first).at(var.first);
                    std::sort(vector_signal.begin(), vector_signal.end());
                    std::vector<double> vector_signal_2 = samples_mass.at(set).at(mass_entry_2->first).at(var.first);
                    std::sort(vector_signal_2.begin(), vector_signal_2.end());
                    Double_t* v_s = vector_signal.data(), *v_s_2 = vector_signal_2.data();
                    kolmogorov[mass_pair]  = TMath::KolmogorovTest(static_cast<int>(vector_signal.size()), v_s,
                                                                   static_cast<int>(vector_signal_2.size()), v_s_2, "");
                }
            }
            CreateMatrixHistosCompatibility(samples_mass.at(set), kolmogorov, var.first+"_Signal_compatibility_KS", directory);
        }
    }

    SetNamesVar FindBestVariables(SampleId mass, std::map<std::string, std::string>& eliminated, const ChannelSpin& set) const
    {
        std::string spin;
        if (set.spin == 0) spin = "Radion";
        else if (set.spin == 2) spin = "Graviton";
        else if (set.spin == SM_spin) spin = "SM";
        else spin = "Bkg";
        static constexpr double threashold_mi = 0.8;
        ChannelSpin chsp_bkg(set.channel, bkg_spin);
        const NameElement& JSDivergenceND = JSDivergenceSB.at(set).at(mass);
        const NameElement& mutual_matrix_signal = mutual_matrix.at(set).at(mass);
        const NameElement& mutual_matrix_bkg = mutual_matrix.at(chsp_bkg).at(SampleType::Bkg_TTbar);
        SetNamesVar selected, not_corrected;
        VectorName_ND JSDivergence_vector(JSDivergenceND.begin(), JSDivergenceND.end());

        std::ofstream best_entries_file("Best_entries_"+ToString(mass)+"_"+std::to_string(args.set())+"_"
                                        +std::to_string(args.number_sets())+set.channel+"_spin"+spin+args.suffix()+".csv", std::ofstream::out);
        best_entries_file<<","<<","<<"JSD_sb"<<","<<"MI_sgn"<<","<<"Mi_bkg"<<std::endl;

        static const std::map<std::pair<size_t, size_t>, std::string> prefixes = {
            { { 0, 0 }, "" }, { { 0, 1 }, "*" }, { { 1, 0 }, "-" }
        };
        while(selected.size() < args.number_variables() && JSDivergence_vector.size()) {
            std::sort(JSDivergence_vector.begin(), JSDivergence_vector.end(), [](const auto& el1, const auto& el2){
                return el1.second > el2.second;
            });
            auto best_entry = JSDivergence_vector.front();
            for(const auto& name : best_entry.first) {
                if (!mutual_matrix_bkg.count(name)) continue;
                const auto counts = std::make_pair(selected.count(name), not_corrected.count(name));
                const auto& prefix = prefixes.at(counts);
                best_entries_file << prefix << name << prefix << ",";
                if(counts.first || counts.second)
                    continue;
                const double JS_1d = JSDivergenceND.at(name);
                for(auto& entry : JSDivergence_vector) {
                    if(entry.first.count(name))
                        entry.second -= JS_1d;
                }
                for(const auto& other_entry : samples_mass.at(set).at(mass)) {
                    if(other_entry.first == name) continue;
                    if (!mutual_matrix_bkg.count(other_entry.first)) continue;
                    const Name_ND names({name, other_entry.first});
                    if(mutual_matrix_signal.at(names) < threashold_mi && mutual_matrix_bkg.at(names) < threashold_mi){
                        eliminated[other_entry.first] = "MI-"+name;
                        not_corrected.insert(other_entry.first);
                    }
                }
                selected.insert(name);
            }
            if (best_entry.first.IsSubset(selected)){
                    best_entries_file << JSDivergenceND.at(best_entry.first) << "," << mutual_matrix_signal.at(best_entry.first)
                           << "," << mutual_matrix_bkg.at(best_entry.first) ;
            }
            best_entries_file<<std::endl;
            JSDivergence_vector = CopySelectedVariables(JSDivergence_vector, best_entry.first, not_corrected);
        }
        best_entries_file.close();
        return selected;
    }

    void InformationTable(const SampleId& mass_entry, const SetNamesVar& selected, const std::map<std::string, std::string>& eliminated, const ChannelSpin& set) const
    {
        std::string spin;
        if (set.spin == 0) spin = "Radion";
        else if (set.spin == 2) spin = "Graviton";
        else if (set.spin == SM_spin) spin = "SM";
        else spin = "Bkg";
        ChannelSpin chsp_bkg(set.channel, bkg_spin);
        auto mass = ToString(mass_entry);
        std::ofstream of("InformationTable_"+ToString(mass)+"_"+std::to_string(args.set())+"_"+std::to_string(args.number_sets())
                         +set.channel+"_spin"+spin+args.suffix()+".csv", std::ofstream::out);
        of<<"Var_1"<<","<<"Var_2"<<","<<"JSD_ND"<<","<<"JSD_12-(JSD_1+JSD_2)"<<","<<"ScaledMI_Signal_12" <<","<<"ScaledMI_Bkg_12"<<","<<"selected"<<","<<","<<"eliminated by"<<std::endl;
        for(auto& entry : JSDivergenceSB.at(set).at(mass_entry)){
             if (entry.first.size() == 1){
                 auto name_1 = entry.first.at(0);
                 bool selected_1 = selected.count(name_1);
                 if (!JSDivergenceSB.at(set).at(mass_entry).count(name_1)) continue;
                 of<<name_1<<","<<"   "<<","<<JSDivergenceSB.at(set).at(mass_entry).at(name_1)<<","<<","<<","<<","<<selected_1;
                 if(eliminated.count(name_1)) of<<","<<","<<eliminated.at(name_1);
                 of<<std::endl;
             }
             else  if (entry.first.size() == 2){
                 auto name_1 = entry.first.at(0);
                 auto name_2 = entry.first.at(1);
                 bool selected_1 = selected.count(name_1);
                 bool selected_2 = selected.count(name_2);
                 if (!mutual_matrix.at(set).at(mass_entry).count(entry.first) || !mutual_matrix.at(chsp_bkg).at(SampleType::Bkg_TTbar).count(entry.first)) continue;
                 of<<name_1<<","<<name_2<<","<<JSDivergenceSB.at(set).at(mass_entry).at(entry.first)<<","
                  <<JSDivergenceSB.at(set).at(mass_entry).at(entry.first)-(JSDivergenceSB.at(set).at(mass_entry).at(name_1)+JSDivergenceSB.at(set).at(mass_entry).at(name_2))
                  <<","<<mutual_matrix.at(set).at(mass_entry).at(entry.first)<<","<<mutual_matrix.at(chsp_bkg).at(SampleType::Bkg_TTbar).at(entry.first)<<","<<selected_1<<","<<selected_2;
                 if(eliminated.count(name_1)) of<<","<<eliminated.at(name_1);
                 if(eliminated.count(name_2)) of<<","<<eliminated.at(name_2);
                 of<<std::endl;
             }
        }
        of.close();
    }

    void CreateMatrixIntersection(const SampleIdSetNamesVar& mass_selected, const ChannelSpin& set) const
    {
        std::string spin;
        if (set.spin == 0) spin = "Radion";
        else if (set.spin == 2) spin = "Graviton";
        else if (set.spin == SM_spin) spin = "SM";
        else if (set.spin == bkg_spin) spin = "Bkg";
        auto directory_set = root_ext::GetDirectory(*outfile, set.channel+spin);
        std::cout<< "Intersection" << std::endl;
        int bin = static_cast<int>(mass_selected.size());
        std::cout<<set.channel << " " << spin<< " "<< bin << std::endl;
        auto matrix_intersection = std::make_shared<TH2D>("Intersection", "Intersection of variables", bin, 0, bin, bin, 0, bin);
        matrix_intersection->SetXTitle("mass");
        matrix_intersection->SetYTitle("mass");
        int k = 1;
        for(const auto& mass_1: mass_selected){
            int j = 1;
            std::string mass = ToString(mass_1.first);
            matrix_intersection->GetXaxis()->SetBinLabel(k, (mass).c_str());
            matrix_intersection->GetYaxis()->SetBinLabel(k, (mass).c_str());
            for(const auto& mass_2: mass_selected){
                VecVariables intersection;
                std::set_intersection(mass_1.second.begin(), mass_1.second.end(), mass_2.second.begin(), mass_2.second.end(),
                                        std::back_inserter(intersection));
                matrix_intersection->SetBinContent(k, j, intersection.size());
                j++;
            }
            k++;
        }
        root_ext::WriteObject(*matrix_intersection, directory_set);;
    }

    void VariablesSelection(const ChannelSpin& set) const
    {
        SampleIdSetNamesVar mass_selected;
        std::string spin;
        if (set.spin == 0) spin = "Radion";
        else if (set.spin == 2) spin = "Graviton";
        else if (set.spin == SM_spin) spin = "SM";
        else spin = "Bkg";
        std::ofstream ofs("SelectedVariable_"+std::to_string(args.set())+"_"+std::to_string(args.number_sets())
                          +set.channel+"_spin"+spin+args.suffix()+".csv", std::ofstream::out);
        for (const auto& mass_entry : samples_mass.at(set)){
            if (mass_entry.first.IsBackground())
                continue;
            ofs << ToString(mass_entry.first) << std::endl;
            std::cout << ToString(mass_entry.first) << " " << std::flush;
            std::map<std::string, std::string> eliminated;

            mass_selected[mass_entry.first] =  FindBestVariables(mass_entry.first, eliminated, set);
            for(auto& entry_selected: mass_selected[mass_entry.first]){
                    ofs << entry_selected<<",";
            }
            ofs << std::endl;

            InformationTable(mass_entry.first, mass_selected.at(mass_entry.first), eliminated, set);
        }
        CreateMatrixIntersection(mass_selected, set);
    }

    void LoadSkimmedData(const ChannelSpin& set)
    {
        std::cout << set.channel << set.spin <<std::endl;
        for(const SampleEntry& entry:samples)
        {
            if ( entry.spin != set.spin) continue;
            auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
            auto tuple = ntuple::CreateEventTuple(set.channel, input_file.get(), true, ntuple::TreeState::Skimmed);
            Long64_t tot_entries = 0;
            for(const Event& event : *tuple) {
                if(tot_entries >= args.number_events()) break;
                LorentzVectorE_Float bb = event.jets_p4[0] + event.jets_p4[1];
                if (args.suffix() == "_ANcut"){
                    if (!cuts::hh_bbtautau_2016::hh_tag::m_hh_window().IsInside(event.SVfit_p4.mass(),bb.mass())) continue;
                }
                if (entry.id == SampleType::Bkg_TTbar && event.file_desc_id>=2) continue;
                if (entry.id == SampleType::Sgn_NonRes && event.file_desc_id!=0) continue;
                auto eventInfoPtr =  analysis::MakeEventInfo(Parse<Channel>(set.channel) ,event) ;
                EventInfoBase& eventbase = *eventInfoPtr;
                if (args.suffix() == "_newcut"){
                    if (!cuts::hh_bbtautau_2016::hh_tag::new_m_hh_window().IsInside(eventbase.GetHiggsTTMomentum(false).M(),bb.mass())) continue;
                }
                vars.AddEvent(eventbase, entry.id, entry.spin, entry.weight);
                tot_entries++;
            }
            std::cout << entry << " number of events: " << tot_entries << std::endl;
        }
        samples_mass[set] = vars.GetSampleVariables(set.channel, set.spin);

        TimeReport();
    }

    void TimeReport(bool tot = false) const
    {
        reporter->TimeReport(tot);
    }

    void Run()
    {
        run::ThreadPull threads(args.number_threads());

        std::cout << "Bandwidth, Mutual Information, JensenShannon" << std::endl;
        for (const auto& s: set){
            std::cout<<s.channel<< "  " << s.spin<<std::endl;
            LoadSkimmedData(s);
            for (const auto& sample: samples_mass[s]){
                std::cout<<"----"<<sample.first.sampleType<<"   "<<sample.first.mass<<"----"<<" entries: "<<sample.second.at("pt_l1").size()<<std::endl;
                std::stringstream ss;
                ss << std::fixed << std::setprecision(0) << s.spin;
                std::string spin = ss.str();
                bandwidth[s][sample.first] = Read_csvfile(args.optband_folder()+"/OptimalBandwidth"+ToString(sample.first)+"_"+s.channel+"_spin"+spin+args.suffix()+".csv", vars.GetDisabledVars());
                mutual_matrix[s][sample.first] = Read_csvfile(args.mutual_folder()+"/MutualInformationDistance"+ToString(sample.first)+"_"+s.channel+"_spin"+spin+args.suffix()+".csv", vars.GetDisabledVars());
                JSDivergenceSB[s][sample.first] = Read_csvfile(args.jsd_folder()+"/JensenShannonDivergenceSB"+ToString(sample.first)+"_"+s.channel+"_spin"+spin+args.suffix()+".csv", vars.GetDisabledVars());
            }
        }

        for (const auto& s: set){
            std::string spin;
            if (s.spin == 0) spin = "Radion";
            else if (s.spin == 2) spin = "Graviton";
            else if (s.spin == SM_spin) spin = "SM";
            else spin = "Bkg";
            ChannelSpin chsp_bkg(s.channel, bkg_spin);
            std::cout << s.channel << "  " << spin << std::endl;

            auto directory_set = root_ext::GetDirectory(*outfile, s.channel+spin);

            auto directory_mutinf = root_ext::GetDirectory(*directory_set, "MutualInformation");
            auto directory_jensenshannon =root_ext::GetDirectory(*directory_set, "JensenShannonDivergence");
            auto directory_jenshan_sgnlbkg = root_ext::GetDirectory(*directory_jensenshannon, "Signal_Background");
            auto directory_jenshan_distrib = root_ext::GetDirectory(*directory_jenshan_sgnlbkg,"Distribution");

            for (const auto& sample: samples_mass[s]){
                std::cout<<"----"<<ToString(sample.first)<<"----"<<" entries: "<<sample.second.at("pt_l1").size()<<std::endl;
                std::cout<<"correlation  " << std::flush;
                correlation_matrix[s][sample.first] = Correlation(sample.second);
                if (sample.first.sampleType == SampleType::Bkg_TTbar){
                    TimeReport();
                    continue;
                }
                std::cout<<"mutual plot  "<< std::flush;
                MutualHisto(sample.first, mutual_matrix.at(s).at(sample.first), mutual_matrix.at(chsp_bkg).at(SampleId{SampleType::Bkg_TTbar,0}), directory_mutinf);
                std::cout<<"JSD_sgn_bkg  "<< std::flush;
                DistributionJSD_SB(sample.first, directory_jenshan_distrib, s);
                TimeReport();
            }

            std::cout <<"Mutual histos" << std::endl;
            auto directory_mutinf_matrix = root_ext::GetDirectory(*directory_mutinf, "Matrix");
            CreateMatrixHistos(samples_mass.at(s), mutual_matrix.at(s), "MI", directory_mutinf_matrix);
            std::cout <<"Correlation histos" << std::endl;
            auto directory_correlation = root_ext::GetDirectory(*directory_set, "Correlation");
            CreateMatrixHistos(samples_mass.at(s), correlation_matrix.at(s), "correlation", directory_correlation, true);

            std::cout <<"MI-Correlation histos" << std::endl;
            auto directory_mi_corr = root_ext::GetDirectory(*directory_set, "MI_Correlation");
            for(const auto& sample : samples_mass.at(s)){
                std::string mass ="MI_Correlation_"+ToString(sample.first);
                auto histo = std::make_shared<TH2D>(mass.c_str(), mass.c_str(),50,0,1,50,0,1);
                histo->SetXTitle("MID");
                histo->SetYTitle("Correlation");
                for (const auto& pair: mutual_matrix.at(s).at(sample.first)){
                    if (!correlation_matrix.at(s).at(sample.first).count(pair.first)) continue;
                    histo->Fill(pair.second, std::abs(correlation_matrix.at(s).at(sample.first).at(pair.first)));
                }
                root_ext::WriteObject(*histo, directory_mi_corr);
            }
            auto directory_jenshan_matrix = root_ext::GetDirectory(*directory_jenshan_sgnlbkg, "Matrix");
            std::cout <<"JSD histos" << std::endl;
            CreateMatrixHistos(samples_mass.at(s), JSDivergenceSB.at(s), "JSD", directory_jenshan_matrix);
            std::cout<<"Plot Jensen Shannon"<<std::endl;
            auto directory_jenshan_allvars = root_ext::GetDirectory(*directory_jenshan_sgnlbkg, "All_variables");
            PlotJensenShannonSB(directory_jenshan_allvars,s);

            if (s.spin == 0 || s.spin == 2){
                std::cout<<"Signal Compatibility Jensen Shannon  "<<std::flush;
                auto directory_jenshan_sgnlsgnl = root_ext::GetDirectory(*directory_jensenshannon, "Signal_Signal");
                JensenShannonSignalCompatibility(directory_jenshan_sgnlsgnl,s);
                TimeReport();
                std::cout<<"Signal Compatibility Kolmogorov"<<std::flush;
                auto directory_ks = root_ext::GetDirectory(*directory_set, "Kolmogorov");
                KolmogorovSignalCompatibility(directory_ks,s);
                TimeReport();
                std::cout<<"Selection variables"<<std::endl;
                VariablesSelection(s);
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

PROGRAM_MAIN(analysis::mva_study::MvaClassification, Arguments) // definition of the main program function
