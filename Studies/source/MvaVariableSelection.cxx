/*! Study of correlation matrix, mutual information and Jensen Shannon Divergence to search
 * best varia bles for each range of mass.
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
#include "hh-bbtautau/Analysis/include/MvaMethods.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "hh-bbtautau/Analysis/include/MvaConfigurationReader.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(unsigned, number_threads);
    REQ_ARG(size_t, number_variables);
    OPT_ARG(Long64_t, number_events, 1000000);
    OPT_ARG(size_t, number_sets, 0);
    OPT_ARG(size_t, set, 0);
    OPT_ARG(uint_fast32_t, seed, 10000);
};

namespace analysis {
namespace mva_study{

using clock = std::chrono::system_clock;

namespace {
Range<int> low_mass(250, 280),  medium_mass1(300, 350), medium_mass2(400, 550), high_mass(600,900);
std::vector<Range<int>> ranges{low_mass, medium_mass1, medium_mass2, high_mass};
}

class VariableDistribution {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    SampleIdVarData samples_mass;
    SampleIdNameElement bandwidth, mutual_matrix, correlation_matrix, JSDivergenceSB;

    VariableDistribution(const Arguments& _args): args(_args),
              outfile(root_ext::CreateRootFile(args.output_file())), vars(args.number_sets(), args.seed()),
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

        samples = samples_list.at("inputs").files;
    }

    NameElement CorrelationSelected(const VarData& sample_vars, const SetNamesVar& selected){
        NameElement corr_matrix;
        NameElementFuture corr_matrix_future;
        for(const auto& var_1 : selected) {
            for(const auto& var_2 : selected) {
                if (var_2 < var_1) continue;
                corr_matrix_future[Name_ND{var_1, var_2}] = run::async(stat_estimators::Correlation<double>,
                                                                sample_vars.at(var_1), sample_vars.at(var_2), 1.);
            }
        }
        for(auto& var : corr_matrix_future) {
            corr_matrix[var.first] = var.second.get() *  100;
        }
        return corr_matrix;
    }


    SetNamesVar FindBestRangeVariables(const Range<int>& range, const std::map<SampleId, double>& max_distance,
                                      std::map<SampleId, VectorName_ND> JSDivergence_vector) const{

        static constexpr double threashold_mi = 0.8;
        SetNamesVar selected, not_corrected;
        std::ofstream best_entries_file("Best_entries_Range"+std::to_string(range.min())+"_"+std::to_string(range.max())+"_"+std::to_string(args.set())+"_"+std::to_string(args.number_sets())+args.tree_name()+".csv", std::ofstream::out);
        best_entries_file << "," << "," << "JSD_sb" << "," << "MI_sgn" << "," << "Mi_bkg" << std::endl;
        const SampleId minSample{SampleType::Sgn_Res, range.min()};
        auto& JSDivergence_min = JSDivergence_vector.at(minSample);
        static const std::map<std::pair<size_t, size_t>, std::string> prefixes = {
            { { 0, 0 }, "" }, { { 0, 1 }, "*" }, { { 1, 0 }, "-" }
        };
        while(selected.size() < args.number_variables() && JSDivergence_min.size()) {
            std::sort(JSDivergence_min.begin(), JSDivergence_min.end(), [](const auto& el1, const auto& el2){
                return el1.second > el2.second;
            });
            VectorName_ND distance;
            for (const auto& var_entry: JSDivergence_min){
                double d = 0;
                for (const auto& entry : JSDivergenceSB){
                    if (entry.first.IsSM() || !range.Contains(entry.first.mass)) continue;
                    d+= entry.second.at(var_entry.first) / max_distance.at(entry.first);
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
                for (const auto& entry : JSDivergenceSB){
                    if (entry.first.IsSM() || !range.Contains(entry.first.mass)) continue;
                    const double JS_1d = entry.second.at(name);
                    for(auto& var : JSDivergence_vector.at(entry.first)) {
                        if(var.first.count(name))
                            var.second -= JS_1d;
                    }
                }
                for (const auto& entry : JSDivergenceSB){
                    if (entry.first.IsSM() || !range.Contains(entry.first.mass)) continue;
                    for(const auto& other_entry : samples_mass.at(entry.first)) {
                        if(other_entry.first == name) continue;
                        const Name_ND names({name, other_entry.first});
                        if(mutual_matrix.at(entry.first).at(names) < threashold_mi && mutual_matrix.at(SampleType::Bkg_TTbar).at(names) < threashold_mi){
                            not_corrected.insert(other_entry.first);
                        }
                    }
                }
                selected.insert(name);
            }
            if (best_entry.first.IsSubset(selected)){
                    best_entries_file << JSDivergenceSB.at(minSample).at(best_entry.first) << "," <<
                           mutual_matrix.at(minSample).at(best_entry.first) << "," <<
                           mutual_matrix.at(SampleType::Bkg_TTbar).at(best_entry.first);;
            }
            best_entries_file << std::endl;
            JSDivergence_min = CopySelectedVariables(JSDivergence_min, best_entry.first, not_corrected);
        }
       return selected;
    }

    void CreateMatrixIntersection(const SampleIdSetNamesVar& range_selected) const
    {
        std::cout << std::endl << "Intersection ranges" << std::endl;
        auto matrix_intersection = std::make_shared<TH2D>("Intersection", "Intersection of variables", ranges.size(), 0, ranges.size(),
                                                          ranges.size(), 0, ranges.size());
        matrix_intersection->SetXTitle("range");
        matrix_intersection->SetYTitle("range");
        int c_label = 1;
        int k_row = 1;
        for (const auto& range1: ranges){
            int j_column = 1;
            for (const auto& range2: ranges){
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
        TimeReport();
    }

    SampleIdSetNamesVar VariablesSelection(){
        SampleIdSetNamesVar range_selected;
        std::map<SampleId, VectorName_ND> JSDivergence_vector;
        std::map<SampleId, double> max_distance;
        for (const auto& entry: JSDivergenceSB){
            VectorName_ND  vector(entry.second.begin(), entry.second.end());
            JSDivergence_vector[entry.first] = vector;
            std::sort(JSDivergence_vector.at(entry.first).begin(), JSDivergence_vector.at(entry.first).end(),
                      [](const auto& el1, const auto& el2){
                return el1.second > el2.second;
            });
            max_distance[entry.first] = JSDivergence_vector.at(entry.first).front().second;
        }
        for (const auto& range: ranges){
            std::ofstream list_variables(("Selected_range"+std::to_string(range.min())+"_"+std::to_string(range.max())+"_"+std::to_string(args.set())+"_"+std::to_string(args.number_sets())+args.tree_name()+".txt").c_str(), std::ofstream::out);
            range_selected[SampleId{SampleType::Sgn_Res, range.min()}] = FindBestRangeVariables(range, max_distance, JSDivergence_vector);
            std::cout<<" range: "<<range.min()<<"-"<<range.max();
            std::cout.flush();
            for (const auto& name : range_selected[SampleId{SampleType::Sgn_Res, range.min()}]){
                list_variables << name << std::endl;
            }
            list_variables << std::endl;
        }
        TimeReport();
        CreateMatrixIntersection(range_selected);
        return range_selected;
    }

    void KolmogorovSignalPlotSelected(const SetNamesVar& selected, const Range<int>& range, TDirectory* directory){
        std::map<std::string, std::shared_ptr<TGraph>> plot;
        int i = 0;
        for (const auto& mass_entry: samples_mass){
            if (!range.Contains(mass_entry.first.mass))
                continue;
            for (const auto& var : selected){
                std::vector<double> vector_signal = samples_mass.at(mass_entry.first).at(var);
                std::sort(vector_signal.begin(), vector_signal.end());
                std::vector<double> vector_signal_2 = samples_mass.at(SampleId{SampleType::Sgn_Res, range.min()}).at(var);
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
            root_ext::WriteObject(*plot[var], directory);
        }
    }

    void LoadData(){
        for(const SampleEntry& entry : samples)
        {
            if ( entry.channel != "" && args.tree_name() != entry.channel ) continue;
            auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
            EventTuple tuple(args.tree_name(), input_file.get(), true, {} , GetMvaBranches());
            Long64_t tot_entries = 0;
            for(Long64_t current_entry = 0; tot_entries < args.number_events() && current_entry < tuple.GetEntries(); ++current_entry) {
                tuple.GetEntry(current_entry);
                const Event& event = tuple.data();
                if (static_cast<EventEnergyScale>(event.eventEnergyScale) != EventEnergyScale::Central || (event.q_1+event.q_2) != 0 || event.jets_p4.size() < 2
                    || event.extraelec_veto == true || event.extramuon_veto == true || event.jets_p4[0].eta() > cuts::btag_2016::eta
                    || event.jets_p4[1].eta() > cuts::btag_2016::eta)
                    continue;
                auto bb = event.jets_p4[0] + event.jets_p4[1];
                if (!cuts::hh_bbtautau_2016::hh_tag::IsInsideEllipse(event.SVfit_p4.mass(), bb.mass()))
                    continue;
                tot_entries++;
                vars.AddEvent(event, entry.id, entry.weight);
            }
            std::cout << entry << " number of events: " << tot_entries << std::endl;
        }
        TimeReport();
        samples_mass = vars.GetSampleVariables(args.set());
    }

    void Run()
    {
        run::ThreadPull threads(args.number_threads());
        LoadData();
        std::cout << "n.variabili: " << samples_mass.at(SampleType::Bkg_TTbar).size() << std::endl;
        std::cout << "n.masse segnale: " << samples_mass.size() - 1 << " + " << 1 << " fondo."<<std::endl;

        std::map<Name_ND, std::shared_ptr<TGraph>> plot_jsd;
        int i = 0;
        auto directory_mutinf = root_ext::GetDirectory(*outfile, "MutualInformation");
        for (const auto& sample: samples_mass){
            std::cout<<"----"<<ToString(sample.first)<<"----"<<" entries: "<<sample.second.at("pt_l1").size()<<std::endl;
            std::cout<<"bandwidth  "<<std::flush;
            bandwidth[sample.first] = OptimalBandwidth(sample.second);
            std::cout<<"mutual  ";
            std::cout.flush();
            mutual_matrix[sample.first] = Mutual(sample.second, bandwidth.at(sample.first));
            if ( sample.first.IsBackground() ){
                TimeReport();
                continue;
            }
            std::cout<<"Jensen Shannon Signal Background ";
            JSDivergenceSB[sample.first] = JensenDivergenceSB(sample.second, samples_mass.at(SampleType::Bkg_TTbar), bandwidth.at(sample.first), bandwidth.at(SampleType::Bkg_TTbar));
            for (const auto& var : JSDivergenceSB.at(sample.first)){
                if (!plot_jsd.count(var.first)) plot_jsd[var.first] = CreatePlot("","","","");
                plot_jsd.at(var.first)->SetPoint(i, sample.first.mass, JSDivergenceSB.at(sample.first).at(var.first));
            }
            i++;
            TimeReport();
        }
        TimeReport();

        auto directory_jensenshannon = root_ext::GetDirectory(*outfile,"JensenShannonDivergence");
        std::cout<<"Selection variables"<<std::endl;
        SampleIdSetNamesVar range_selected = VariablesSelection();

        std::cout<<"Signal Compatibility"<<std::endl;
        std::map<std::pair<int,int>, std::map<Name_ND, double>> JSDvars_range_ss;
        std::map<int, std::map<Name_ND, std::shared_ptr<TGraph>>> plot_ss;
        std::map<std::pair<int,int>, std::shared_ptr<TH1D>> histo_distribution;
        auto directory_ks = root_ext::GetDirectory(*outfile,"Kolmogorov");
        auto directory_jenshan_ss = root_ext::GetDirectory(*directory_jensenshannon, "Range_SignalSignal");
        auto directory_plotselected = root_ext::GetDirectory(*directory_jenshan_ss, "Plot_RangeSelected");
        for (const auto& range : ranges){
            KolmogorovSignalPlotSelected(range_selected.at(SampleId{SampleType::Sgn_Res, range.min()}), range, directory_ks);
            std::map<std::pair<int,int>, std::map<Name_ND, std::future<double>>> JSDvars_range_ss_future;
            for (const auto& sample_mass: samples_mass){
                if (!range.Contains(sample_mass.first.mass)) continue;
                std::pair<int, int> mass_pair(range.min(), sample_mass.first.mass);
                histo_distribution[mass_pair] = std::make_shared<TH1D>(("JSDrange_"+std::to_string(range.min())+"_"+std::to_string(sample_mass.first.mass)).c_str(), ("JSDrange_"+std::to_string(range.min())+"_"+std::to_string(sample_mass.first.mass)).c_str(), 50,0,1);
                histo_distribution[mass_pair]->SetXTitle("JSD");
                for (const auto& var_1: range_selected.at(SampleId{SampleType::Sgn_Res, range.min()})){
                    if (!plot_ss[range.min()].count(var_1)){
                        plot_ss[range.min()][var_1] = CreatePlot(("JSD_"+var_1+"_Range"+std::to_string(range.min())+"_"+std::to_string(range.max())).c_str(),
                                                           ("JSD_"+var_1+"_Range"+std::to_string(range.min())+"_"+std::to_string(range.max())).c_str(), "mass", "JSD");

                        plot_ss[range.min()][var_1]->SetLineColor(static_cast<Color_t>(range.min()/100));
                        plot_ss[range.min()][var_1]->SetMarkerColor(1);
                        plot_ss[range.min()][var_1]->SetMarkerSize(1);
                        plot_ss[range.min()][var_1]->SetMarkerStyle(8);
                    }
                    std::vector<const DataVector*> sample_1, sample_2;
                    DataVector band_1, band_2;
                    sample_1.push_back(&samples_mass.at(SampleId{SampleType::Sgn_Res, range.min()}).at(var_1));
                    band_1.push_back(bandwidth.at(SampleId{SampleType::Sgn_Res, range.min()}).at(var_1));
                    sample_2.push_back(&samples_mass.at(sample_mass.first).at(var_1));
                    band_2.push_back( bandwidth.at(sample_mass.first).at(var_1));
                    JSDvars_range_ss_future[mass_pair][var_1]  = run::async(stat_estimators::JensenShannonDivergence_ND<double>,
                                                                          sample_1, sample_2, band_1, band_2);
                    for (const auto& var2: range_selected.at(SampleId{SampleType::Sgn_Res, range.min()})){
                        if (var2 < var_1) continue;
                        Name_ND var_pair{var_1, var2};
                        if (!plot_ss[range.min()].count(var_pair)) {
                            plot_ss[range.min()][var_pair] = std::make_shared<TGraph>();
                            plot_ss[range.min()][var_pair]->SetLineColor(4);
                            plot_ss[range.min()][var_pair]->SetLineWidth(1);
                            plot_ss[range.min()][var_pair]->SetMarkerColor(7);
                            plot_ss[range.min()][var_pair]->SetMarkerSize(1);
                            plot_ss[range.min()][var_pair]->SetMarkerStyle(8);
                            plot_ss[range.min()][var_pair]->SetTitle(("JSD_"+var_1+"_"+var2+"_Range"+std::to_string(range.min())+"_"+std::to_string(range.max())).c_str());
                            plot_ss[range.min()][var_pair]->SetName(("JSD_"+var_1+"_"+var2+"_Range"+std::to_string(range.min())+"_"+std::to_string(range.max())).c_str());
                            plot_ss[range.min()][var_pair]->GetHistogram()->GetXaxis()->SetTitle("mass");
                            plot_ss[range.min()][var_pair]->GetHistogram()->GetYaxis()->SetTitle("JSD");
                        }
                        sample_1.push_back(&samples_mass.at(SampleId{SampleType::Sgn_Res, range.min()}).at(var2));
                        band_1.push_back(bandwidth.at(SampleId{SampleType::Sgn_Res, range.min()}).at(var2));
                        sample_2.push_back(&samples_mass.at(sample_mass.first).at(var2));
                        band_2.push_back( bandwidth.at(sample_mass.first).at(var2));
                        JSDvars_range_ss_future[mass_pair][var_pair] = run::async(stat_estimators::JensenShannonDivergence_ND<double>,
                                                                                            sample_1, sample_2, band_1, band_2);

                        sample_1.erase(sample_1.end() - 1);
                        sample_2.erase(sample_2.end() - 1);
                        band_1.erase(band_1.end() - 1);
                        band_2.erase(band_2.end() - 1);
                    }
                }
            }
            for(auto& mass_pair : JSDvars_range_ss_future) {
                for(auto& name : mass_pair.second) {
                    JSDvars_range_ss[mass_pair.first][name.first] = name.second.get();
                }
            }
        }
        auto directory_distributionselected = root_ext::GetDirectory(*directory_jenshan_ss, "Distribution_RangeSelected");
        int k = 0;
        for(auto& mass_pair : JSDvars_range_ss) {
            if (mass_pair.first.first == mass_pair.first.second) k = 0;
            for(auto& name : mass_pair.second) {
                histo_distribution.at(mass_pair.first)->Fill(JSDvars_range_ss.at(mass_pair.first).at(name.first));
                plot_ss.at(mass_pair.first.first).at(name.first)->SetPoint(k, mass_pair.first.second, JSDvars_range_ss.at(mass_pair.first).at(name.first));
            }
            k++;
        }
        for (const auto& range : ranges){
            auto directory_plotrange = root_ext::GetDirectory(*directory_plotselected, ("Range"+std::to_string(range.min())+"_"+std::to_string(range.max())).c_str());
            for (const auto& name : plot_ss.at(range.min())){
                root_ext::WriteObject(*name.second, directory_plotrange);
            }
        }
        for (const auto& mass_pair : histo_distribution){
            root_ext::WriteObject(*mass_pair.second, directory_distributionselected);
        }
        TimeReport();

        std::cout<<"Union sample"<<std::endl;
        SampleIdNameElement bandwidth_range, mutual_matrix_range, correlation_matrix_range_signal, correlation_matrix_range_bkg;
        SampleIdVarData samples_range;
        for (const auto& range : ranges){
            for (const auto& var: range_selected.at(SampleId{SampleType::Sgn_Res, range.min()})){
                for (const auto& sample_mass : samples_mass){
                    if (!range.Contains(sample_mass.first.mass)) continue;

                    for (const auto& entry : sample_mass.second.at(var)){
                        samples_range[SampleId{SampleType::Sgn_Res, range.min()}][var].push_back(entry);
                    }
                }
            }
        }
        auto directory_correlation = root_ext::GetDirectory(*outfile,"Correlation");
        for (const auto& sample: samples_range){
            std::cout<<"----Range"<<ToString(sample.first)<<"----"<<" entries: "<<sample.second.begin()->second.size()<<std::endl;
            std::cout<<"correlation  ";
            std::cout.flush();
            correlation_matrix_range_signal[sample.first] = CorrelationSelected(sample.second, range_selected.at(sample.first));
            correlation_matrix_range_bkg[sample.first] = CorrelationSelected(samples_mass.at(SampleType::Bkg_TTbar), range_selected.at(sample.first));
            int bin = static_cast<int>(range_selected.at(SampleId{SampleType::Sgn_Res, sample.first.mass}).size());
            auto matrix = std::make_shared<TH2D>(("Bkg_"+std::to_string(sample.first.mass)).c_str(),("Bkg_"+std::to_string(sample.first.mass)).c_str(),bin,0,bin,bin,0,bin);
            int i = 1;
            for(const auto& var_1 : range_selected.at(sample.first)) {
                int j = 1;
                matrix->GetXaxis()->SetBinLabel(i, (var_1).c_str());
                matrix->GetYaxis()->SetBinLabel(i, (var_1).c_str());
                for(const auto& var_2 : range_selected.at(sample.first)) {
                    matrix->SetBinContent(i, j, correlation_matrix_range_bkg.at(sample.first).at(Name_ND{var_1, var_2}) * 100);
                    j++;
                }
                i++;
            }
            root_ext::WriteObject(*matrix, directory_correlation);
            std::cout<<"bandwidth  ";
            std::cout.flush();;
            bandwidth_range[sample.first] = OptimalBandwidth(sample.second);
            std::cout<<"mutual  ";
            std::cout.flush();;
            mutual_matrix_range[sample.first] = Mutual(sample.second, bandwidth.at(sample.first));
            if ( sample.first.sampleType == SampleType::Bkg_TTbar){
                std::cout<<std::endl;
                continue;
            }
            else {
                std::cout<<"mutual plot  ";
                MutualHisto(sample.first, mutual_matrix.at(sample.first), mutual_matrix.at(SampleType::Bkg_TTbar), directory_mutinf);
            }
            TimeReport();
        }
        std::cout <<"Mutual histos" << std::endl;
        directory_mutinf->mkdir("Matrix");
        auto directory_mutinf_matrix = root_ext::GetDirectory(*directory_mutinf, "Matrix");
        CreateMatrixHistos(samples_range, mutual_matrix_range, "MI", directory_mutinf_matrix);
        std::cout <<"Correlation histos" << std::endl;
        CreateMatrixHistos(samples_range, correlation_matrix_range_signal, "Signal", directory_correlation);
        TimeReport();

        std::cout<<"Jensen Shannon union Signal Backgound"<<std::endl;
        std::map<int, std::map<Name_ND, double>> JSDvars_range_sb;
        std::map<Name_ND, std::shared_ptr<TGraph>> plot_sb;
        std::map<int, std::shared_ptr<TH1D>> histo_distribution_sb;
        for (const auto& range : ranges){
            std::map<int, std::map<Name_ND, std::future<double>>> JSDvars_range_sb_future;
            for (const auto& sample: samples_range){
                if (!range.Contains(sample.first.mass)) continue;
                histo_distribution_sb[range.min()] = std::make_shared<TH1D>(("JSDrange_"+std::to_string(range.min())+"_SignalBkg").c_str(),
                                                                          ("JSDrange_"+std::to_string(range.min())+"_SignalBkg").c_str(),
                                                                          50,0,1);
                histo_distribution_sb[range.min()]->SetXTitle("JSD");
                for (const auto& var: range_selected.at(SampleId{SampleType::Sgn_Res, range.min()})){
                    if (!plot_sb.count(var)){
                        plot_sb[var] = CreatePlot(("JSD_"+var+"_SignalBkg").c_str(), ("JSD_"+var+"_SignalBkg").c_str(), "mass", "JSD");
                        plot_sb[var]->SetLineColor(kRed+2);
                        plot_sb[var]->SetMarkerColor(2);
                        plot_sb[var]->SetMarkerSize(1);
                        plot_sb[var]->SetMarkerStyle(5);
                    }
                    std::vector<const DataVector*> sample_1, sample_2;
                    DataVector band_1, band_2;
                    sample_1.push_back(&sample.second.at(var));
                    band_1.push_back(bandwidth_range.at(SampleId{SampleType::Sgn_Res, range.min()}).at(var));
                    sample_2.push_back(&samples_mass.at(SampleType::Bkg_TTbar).at(var));
                    band_2.push_back(bandwidth.at(SampleType::Bkg_TTbar).at(var));
                    JSDvars_range_sb_future[range.min()][var]  = run::async(stat_estimators::JensenShannonDivergence_ND<double>,
                                                                          sample_1, sample_2, band_1, band_2);
                    for (const auto& var2: range_selected.at(SampleId{SampleType::Sgn_Res, range.min()})){
                        if (var2 <= var) continue;
                        Name_ND var_pair{var, var2};
                        if (!plot_sb.count(Name_ND{var, var2})) {
                            plot_sb[var_pair] = std::make_shared<TGraph>();
                            plot_sb[var_pair]->SetLineColor(kRed+2);
                            plot_sb[var_pair]->SetLineWidth(1);
                            plot_sb[var_pair]->SetMarkerColor(2);
                            plot_sb[var_pair]->SetMarkerSize(1);
                            plot_sb[var_pair]->SetMarkerStyle(5);
                            plot_sb[var_pair]->SetTitle(("JSD_"+var+"_"+var2+"_SignalBkg").c_str());
                            plot_sb[var_pair]->SetName(("JSD_"+var+"_"+var2+"_SignalBkg").c_str());
                            plot_sb[var_pair]->GetHistogram()->GetXaxis()->SetTitle("mass");
                            plot_sb[var_pair]->GetHistogram()->GetYaxis()->SetTitle("JSD");
                        }
                        sample_1.push_back(&sample.second.at(var2));
                        band_1.push_back(bandwidth_range.at(SampleId{SampleType::Sgn_Res, range.min()}).at(var2));
                        sample_2.push_back(&samples_mass.at(SampleType::Bkg_TTbar).at(var2));
                        band_2.push_back(bandwidth.at(SampleType::Bkg_TTbar).at(var2));
                        JSDvars_range_sb_future[range.min()][var_pair] = run::async(stat_estimators::JensenShannonDivergence_ND<double>,
                                                                                            sample_1, sample_2, band_1, band_2);

                        sample_1.erase(sample_1.end() - 1);
                        sample_2.erase(sample_2.end() - 1);
                        band_1.erase(band_1.end() - 1);
                        band_2.erase(band_2.end() - 1);
                    }
                }
            }
            for(auto& mass_pair : JSDvars_range_sb_future) {
                for(auto& name : mass_pair.second) {
                    JSDvars_range_sb[mass_pair.first][name.first] = name.second.get();
                }
            }
        }
        auto directory_rangebkg = root_ext::GetDirectory(*directory_jensenshannon, "Range_SignalBkg");
        auto directory_plotselected_sb = root_ext::GetDirectory(*directory_rangebkg, "Plot_RangeSelected");
        auto directory_distributionselected_sb = root_ext::GetDirectory(*directory_rangebkg, "Distribution_RangeSelected");
        int w = 0;
        for(auto& mass : JSDvars_range_sb) {
            for(auto& name : mass.second) {
                histo_distribution_sb.at(mass.first)->Fill(JSDvars_range_sb.at(mass.first).at(name.first));
                plot_sb.at(name.first)->SetPoint(w, mass.first, JSDvars_range_sb.at(mass.first).at(name.first));
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
        for (const auto& range: ranges){
            int bin = static_cast<int>(range_selected.at(SampleId{SampleType::Sgn_Res, range.min()}).size());
            auto matrix_jsd_sb = std::make_shared<TH2D>(("JSD_Signal_Bkg_Range"+std::to_string(range.min())+"_"+std::to_string(range.min())).c_str(),("JSD_Signal_Bkg_Range"+std::to_string(range.min())+"_"+std::to_string(range.min())).c_str(), bin, 0, bin, bin, 0, bin);
            int i = 1;
            for (const auto& var1: range_selected.at(SampleId{SampleType::Sgn_Res, range.min()})){
                    matrix_jsd_sb->GetXaxis()->SetBinLabel(i, (var1).c_str());
                    matrix_jsd_sb->GetYaxis()->SetBinLabel(i, (var1).c_str());
                    int jj = 1;
                    for (const auto& var2: range_selected.at(SampleId{SampleType::Sgn_Res, range.min()})){
                        if (var1 == var2){
                            matrix_jsd_sb->SetBinContent(i, i, JSDvars_range_sb.at(range.min()).at(var1));
                        }
                        if (jj > i){
                            matrix_jsd_sb->SetBinContent(i, jj, JSDvars_range_sb.at(range.min()).at(Name_ND{var1, var2}));
                            matrix_jsd_sb->SetBinContent(jj, i, JSDvars_range_sb.at(range.min()).at(Name_ND{var1, var2}));
                        }
                        jj++;
                    }
                    i++;
            }
            root_ext::WriteObject(*matrix_jsd_sb, directory_matrix);
        }

        auto directory_compare = root_ext::GetDirectory(*directory_jensenshannon, "Comparison");
        std::map<Name_ND, TCanvas*> canvas;
        for (const auto & range: ranges){
            for (const auto& selected_name : plot_ss.at(range.min())){
                if (canvas.count(selected_name.first)) continue;
                canvas[selected_name.first] = new TCanvas(selected_name.first.ToString().c_str(), selected_name.first.ToString().c_str(), 10,10,800,800);
                canvas[selected_name.first]->DrawFrame(200,0,1000,1);
                canvas[selected_name.first]->SetGrid();
            }
        }
        for (const auto& selected_name : plot_sb){
            canvas.at(selected_name.first)->cd();
            selected_name.second->Draw("P");
            for (const auto & range: ranges){
                if (!plot_ss.at(range.min()).count(selected_name.first)) continue;
                plot_ss.at(range.min()).at(selected_name.first)->Draw();
            }
            plot_jsd.at(selected_name.first)->Draw();
            root_ext::WriteObject<TCanvas>(*canvas[selected_name.first], directory_compare);
            canvas[selected_name.first]->Close();
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
    std::shared_ptr<TFile> outfile;
    MvaVariablesStudy vars;
    std::shared_ptr<TimeReporter> reporter;
};
}
}

PROGRAM_MAIN(analysis::mva_study::VariableDistribution, Arguments) // definition of the main program function
