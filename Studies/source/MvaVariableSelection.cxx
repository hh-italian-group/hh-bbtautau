/*! Study of correlation matrix and mutual information of BDT variables
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <fstream>
#include <random>
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
#include "hh-bbtautau/Analysis/include/MvaConfiguration.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(unsigned long, number_threads);
    REQ_ARG(size_t, number_variables);
    OPT_ARG(Long64_t, number_events, 1000000);
    OPT_ARG(bool, split, false);
    OPT_ARG(bool, istraining, false);
};

namespace analysis {
namespace Mva_Study{

using clock = std::chrono::system_clock;

Range<int> low_mass(250, 280),  medium_mass(300, 400), high_mass(450,900);
std::vector<Range<int>> ranges{low_mass, medium_mass, high_mass};

NameElement CorrelationSelected(const VarData& sample_vars, const VarNameSet& selected){
    NameElement corr_matrix;
    for(const auto& var_1 : selected) {
        for(const auto& var_2 : selected) {
            if (var_2 < var_1) continue;
            corr_matrix[Name_ND{var_1, var_2}] = stat_estimators::Correlation(sample_vars.at(var_1), sample_vars.at(var_2));
        }
    }
    return corr_matrix;
}


static constexpr double threashold_mi = 0.7;
VarNameSet CheckList_2(int min, int max, const std::map<int, double>& max_distance, const MassVar& samples,
                       const MassNameElement& JSDivergenceND, std::map<int, VectorName_ND> JSDivergence_vector,
                       const MassNameElement& mutual_matrix, size_t number_variables){

    VarNameSet selected, not_corrected;
    std::ofstream ofb("Best_entries_Range"+std::to_string(min)+"_"+std::to_string(max)+".csv", std::ofstream::out);
    while(selected.size() < number_variables && JSDivergence_vector.at(min).size()) {
        std::sort(JSDivergence_vector.at(min).begin(), JSDivergence_vector.at(min).end(), [](const std::pair<Name_ND, double>& el1, const std::pair<Name_ND, double>& el2){
            return el1.second > el2.second;
        });
        VectorName_ND distance;
        for (const auto& var_entry: JSDivergence_vector.at(min)){
            double d = 0;
            for (const auto& mass_entry: JSDivergenceND){
                if (mass_entry.first < min || mass_entry.first > max) continue;
                d+= (max_distance.at(mass_entry.first)-mass_entry.second.at(var_entry.first)) / max_distance.at(mass_entry.first);
            }
            distance.emplace_back(var_entry.first,d);
        }
        std::sort(distance.begin(), distance.end(), [](const std::pair<Name_ND, double>& el1, const std::pair<Name_ND, double>& el2){
        return el1.second < el2.second;
        });
        const auto& best_entry = distance.front();
        for(const auto& name : best_entry.first.names) {
            if(selected.count(name)) {
                ofb<<"''"<<name<<"''"<<",";
                continue;
            }
            for (const auto& mass : JSDivergenceND){
                if (mass.first ==  Signal_SM) continue;
                const double JS_1d = mass.second.at(Name_ND({name}));
                for(auto& entry : JSDivergence_vector.at(mass.first)) {
                    if(entry.first.names.count(name))
                        entry.second -= JS_1d;
                }
            }
            if (not_corrected.count(name)) continue;
            for (const auto& mass : JSDivergenceND){
                if (mass.first ==  Signal_SM) continue;
                for(const auto& other_entry : samples.at(mass.first)) {
                    if(other_entry.first == name) continue;
                    const Name_ND names({name, other_entry.first});
                    if(mutual_matrix.at(mass.first).at(names) < threashold_mi && mutual_matrix.at(Bkg).at(names) < threashold_mi){
                        not_corrected.insert(other_entry.first);
                    }
                }
            }
            selected.insert(name);
            ofb<<name<<",";
        }
        ofb<<std::endl;
        VectorName_ND copy;
        for (const auto& entry : JSDivergence_vector.at(min)){
            if(entry.second <= 0) continue;
            bool consider = true;
            for(const auto& other : entry.first) {
                consider = consider && !not_corrected.count(other);
            }
            if(consider)
                copy.push_back(entry);
        }
        JSDivergence_vector.at(min) = std::move(copy);
        for(auto mass : JSDivergence_vector){
            int c = 0;
            for (auto entry : mass.second){
                if (entry.first.IsEqual(best_entry.first)){
                    JSDivergence_vector.at(mass.first).erase(JSDivergence_vector.at(mass.first).begin() + c);
                }
                c++;
            }
        }
    }
   return selected;
}


MassVarNameSet VariablesSelection(const MassVar& samples, const MassNameElement& JSDivergenceSB,
                        const MassNameElement& mutual_matrix, size_t number_variables){
    MassVarNameSet range_selected;
    auto start = clock::now();
    std::cout<<"Checklist"<<std::endl;
    std::map<int, VectorName_ND> JSDivergence_vector;
    std::map<int, double> max_distance;
    for (const auto& entry_JSD: JSDivergenceSB){
        VectorName_ND  vector(entry_JSD.second.begin(), entry_JSD.second.end());
        JSDivergence_vector[entry_JSD.first] = vector;
        std::sort(JSDivergence_vector.at(entry_JSD.first).begin(), JSDivergence_vector.at(entry_JSD.first).end(), [](const std::pair<Name_ND, double>& el1, const std::pair<Name_ND, double>& el2){
            return el1.second > el2.second;
        });
        max_distance[entry_JSD.first] = JSDivergence_vector.at(entry_JSD.first).front().second;
    }
    std::ofstream ofr("Selected_range.csv", std::ofstream::out);
    for (const auto& range: ranges){
        range_selected[range.min()] = CheckList_2(range.min(), range.max(), max_distance, samples, JSDivergenceSB, JSDivergence_vector, mutual_matrix, number_variables);
        std::cout<<" range: "<<range.min()<<"-"<<range.max()<<std::endl;
        ofr<<"range: "<<range.min()<<"-"<<range.max()<<std::endl;
        for (const auto& name : range_selected[range.min()]){
            ofr<<name<<",";
        }
        ofr<<std::endl;   
    }
    ofr.close();
    auto stop = clock::now();
    std::cout<<std::endl<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;
    return range_selected;
}

void KolmogorovSignalPlotSelected(const MassVar& sample_signal, const VarNameSet& selected,
                          int min, int max, TDirectory* directory){
    std::map<std::string, std::shared_ptr<TGraph>> plot;
     int i = 0;
    for (const auto& mass_entry: sample_signal){
        if (mass_entry.first < min || mass_entry.first > max)
            continue;
        for (const auto& var : selected){
            std::vector<double> vector_signal = sample_signal.at(mass_entry.first).at(var);
            std::sort(vector_signal.begin(), vector_signal.end());
            std::vector<double> vector_signal_2 = sample_signal.at(min).at(var);
            std::sort(vector_signal_2.begin(), vector_signal_2.end());
            Double_t* v_s = vector_signal.data(), *v_s_2 = vector_signal_2.data();
            double k = TMath::KolmogorovTest(vector_signal.size(), v_s, vector_signal_2.size(), v_s_2, "");
            if (plot.count(var) == 0) {
                plot[var] = CreatePlot(("KolmogorovSmirnov_"+var+"_Signal"+std::to_string(min)+"_"+std::to_string(max)).c_str(),
                                       ("Ks_"+var+"_Signal"+std::to_string(min)+"_"+std::to_string(max)).c_str(),
                                       "mass","KS Probability" );
            }
            plot[var]->SetPoint(i,mass_entry.first,k);
        }
        i++;
    }

    for(const auto& var: selected){
        root_ext::WriteObject(*plot[var], directory);
    }
}


class VariableDistribution {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    static const VarNameSet& GetEnabledBranches()
    {

        static const VarNameSet EnabledBranches_read = {
            "eventEnergyScale", "q_1", "q_2", "jets_p4", "extraelec_veto", "extramuon_veto ", "SVfit_p4",
            "pfMET_p4", "p4_1", "p4_2"
        };
        return EnabledBranches_read;
    }

//    static const VarNameSet& GetDisabledBranches()
//    {
//        static const VarNameSet DisabledBranches_read = {
//            "dphi_mumet", "dphi_metsv", "dR_taumu", "mT1", "mT2", "dphi_bbmet", "dphi_bbsv", "dR_bb", "m_bb", "n_jets",
//            "btag_weight", "ttbar_weight",  "PU_weight", "shape_denominator_weight", "trigger_accepts", "trigger_matches",
//            "event.tauId_keys_1","event.tauId_keys_2","event.tauId_values_1","event.tauId_values_2"
//        };
//        return DisabledBranches_read;
//    }

    VariableDistribution(const Arguments& _args): args(_args), samples(ReadConfig(args.cfg_file())),
        outfile(root_ext::CreateRootFile(args.output_file())), split_training_testing(args.split()),
        vars(split_training_testing,UINT_FAST32_MAX, enabled_vars, disabled_vars)
    {
    }

    void Run()
    {
        auto start_tot = clock::now();
        run::ThreadPull threads(args.number_threads());

        MassVar samples_mass;
        auto start = clock::now();
        for(const SampleEntry& entry : samples)
        {
            if ( entry.channel != "" && args.tree_name() != entry.channel ) continue;
            auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
            EventTuple tuple(args.tree_name(), input_file.get(), true, {} , GetEnabledBranches());
            Long64_t tot_entries = 0, current_entry = 0;
            while(tot_entries < args.number_events() && current_entry < tuple.GetEntries()) {
                tuple.GetEntry(current_entry);
                const Event& event = tuple.data();

                if (event.eventEnergyScale != 0 || (event.q_1+event.q_2) != 0 || event.jets_p4.size() < 2
                    || event.extraelec_veto == true || event.extramuon_veto == true || event.jets_p4[0].eta() > cuts::btag_2016::eta
                    || event.jets_p4[1].eta() > cuts::btag_2016::eta){
                    current_entry++;
                    continue;
                }

                LorentzVectorE_Float bb = event.jets_p4[0] + event.jets_p4[1];
                double ellipse_cut = pow(event.SVfit_p4.mass()-116,2)/pow(35.,2) + pow(bb.mass()-111,2)/pow(45.,2);
                if (ellipse_cut>1){
                    current_entry++;
                    continue;
                }

                current_entry++;
                tot_entries++;
                vars.AddEvent(event, entry.mass, entry.weight);
            }
            std::cout << entry << " number of events: " << tot_entries << std::endl;
        }

        auto stop = clock::now();
        std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;

        samples_mass = vars.GetSampleVariables(args.istraining());

        std::cout << "n.variabili: " << samples_mass.at(Bkg).size() << std::endl;
        std::cout << "n.masse segnale: " << samples_mass.size() - 1 << " + " << 1 << " fondo."<<std::endl;

        MassNameElement bandwidth, mutual_matrix;

        std::cout << "Bandwidth and Mutual Information" << std::endl;
        outfile->mkdir("MutualInformation");
        auto directory_mutinf = outfile->GetDirectory("MutualInformation");
        start = clock::now();
        for (const auto& sample_mass: samples_mass){
            std::cout<<"----"<<sample_mass.first<<"----"<<" entries: "<<sample_mass.second.at("pt_l1").size()<<std::endl;
            std::cout<<"bandwidth  ";
            bandwidth[sample_mass.first] = OptimalBandwidth(sample_mass.second);
            std::cout<<"mutual  ";
            mutual_matrix[sample_mass.first] = Mutual(sample_mass.second, bandwidth.at(sample_mass.first));
            std::cout<<std::endl;
        }
        stop = clock::now();
        std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;

        MassNameElement JSDivergenceSB;
        std::cout<<"Jensen Shannon Signal Background"<<std::endl;
        start = clock::now();
        std::map<Name_ND, std::shared_ptr<TGraph>> plot_jsd;
        int i = 0;
        for (const auto& sample_mass: samples_mass){
            if (sample_mass.first == Bkg)
                continue;
            JSDivergenceSB[sample_mass.first] = JensenDivergenceSB(sample_mass.second, samples_mass.at(Bkg), bandwidth.at(sample_mass.first), bandwidth.at(Bkg));
            for (const auto& var : JSDivergenceSB.at(sample_mass.first)){
                if (!plot_jsd.count(var.first)) plot_jsd[var.first] = CreatePlot("","","","");
                plot_jsd.at(var.first)->SetPoint(i, sample_mass.first, JSDivergenceSB.at(sample_mass.first).at(var.first));
            }
            i++;
        }
        stop = clock::now();
        std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;
        outfile->mkdir("JensenShannonDivergence");
        auto directory_jensenshannon = outfile->GetDirectory("JensenShannonDivergence");

        std::cout<<"Selection variables"<<std::endl;
        MassVarNameSet range_selected = VariablesSelection(samples_mass, JSDivergenceSB, mutual_matrix, args.number_variables());

        std::cout<<"Intersection ranges"<<std::endl;
        auto matrix_intersection = std::make_shared<TH2D>("Intersection_ranges","Intersection of variables", ranges.size(), 0, ranges.size(), ranges.size(), 0, ranges.size());
        matrix_intersection->SetXTitle("range");
        matrix_intersection->SetYTitle("range");
        int c = 1;
        for (const auto& range1: ranges){
            for (const auto& range2: ranges){
                if (range1.min() != range2.min()) continue;
                std::string label;
                label = std::to_string(range1.min()) + "-" + std::to_string(range1.max());
                matrix_intersection->GetXaxis()->SetBinLabel(c, (label).c_str());
                matrix_intersection->GetYaxis()->SetBinLabel(c, (label).c_str());
                c++;
            }
        }
        int kk = 1;
        for (const auto& range1 : ranges){
            int jj = 1;
             for (const auto& range2 : ranges){
                 std::vector<std::string> intersection;
                 std::set_intersection(range_selected.at(range1.min()).begin(), range_selected.at(range1.min()).end(), range_selected.at(range2.min()).begin(), range_selected.at(range2.min()).end(),
                                         std::back_inserter(intersection));
                 matrix_intersection->SetBinContent(kk, jj, intersection.size());
                 jj++;
             }
             kk++;
        }
        root_ext::WriteObject(*matrix_intersection, directory_jensenshannon);

        std::map<std::pair<int,int>, std::map<Name_ND, double>> JSDvar_range_ss, JSDpair_range_ss;
        std::map<int, std::map<Name_ND, std::shared_ptr<TGraph>>> plot_ss;
        std::map<std::pair<int,int>, std::shared_ptr<TH1D>> histo_distribution;
        outfile->mkdir("Kolmogorov");
        auto directory_ks = outfile->GetDirectory("Kolmogorov");
        directory_jensenshannon->mkdir("Range_SignalSignal");
        auto directory_jenshan_ss = directory_jensenshannon->GetDirectory("Range_SignalSignal");
        auto directory_plotselected = directory_jenshan_ss->GetDirectory("Plot_RangeSelected");
        if (directory_plotselected == nullptr) {
            directory_jenshan_ss->mkdir("Plot_RangeSelected");
            directory_plotselected = directory_jenshan_ss->GetDirectory("Plot_RangeSelected");
        }
        for (const auto& range : ranges){
            KolmogorovSignalPlotSelected(samples_mass, range_selected.at(range.min()), range.min(), range.max(), directory_ks);
            std::map<std::pair<int,int>, std::map<Name_ND, std::future<double>>> JSDvar_range_ss_future, JSDpair_range_ss_future ;
            for (const auto& sample_mass: samples_mass){
                if (!range.Contains(sample_mass.first)) continue;
                std::pair<int, int> mass_pair(range.min(), sample_mass.first);
                histo_distribution[mass_pair] = std::make_shared<TH1D>(("JSDrange_"+std::to_string(range.min())+"_"+std::to_string(sample_mass.first)).c_str(), ("JSDrange_"+std::to_string(range.min())+"_"+std::to_string(sample_mass.first)).c_str(), 50,0,1);
                histo_distribution[mass_pair]->SetXTitle("JSD");
                for (const auto& var: range_selected.at(range.min())){
                    if (!plot_ss[range.min()].count(Name_ND{var})){
                        plot_ss[range.min()][Name_ND{var}] = CreatePlot(("JSD_"+var+"_Range"+std::to_string(range.min())+"_"+std::to_string(range.max())).c_str(),
                                                           ("JSD_"+var+"_Range"+std::to_string(range.min())+"_"+std::to_string(range.max())).c_str(), "mass", "JSD");

                        plot_ss[range.min()][Name_ND{var}]->SetLineColor(4);
                        plot_ss[range.min()][Name_ND{var}]->SetMarkerColor(7);
                        plot_ss[range.min()][Name_ND{var}]->SetMarkerSize(1);
                        plot_ss[range.min()][Name_ND{var}]->SetMarkerStyle(8);
                    }
                    std::vector<const DataVector*> sample_1, sample_2;
                    DataVector band_1, band_2;
                    sample_1.push_back(&samples_mass.at(range.min()).at(var));
                    band_1.push_back(bandwidth.at(range.min()).at(Name_ND{var}));
                    sample_2.push_back(&samples_mass.at(sample_mass.first).at(var));
                    band_2.push_back( bandwidth.at(sample_mass.first).at(Name_ND{var}));
                    JSDvar_range_ss_future[mass_pair][Name_ND{var}]  = run::async(stat_estimators::JensenShannonDivergence_ND<double>,
                                                                          sample_1, sample_2, band_1, band_2);
                    for (const auto& var2: range_selected.at(range.min())){
                        if (var2 < var) continue;
                        if (!plot_ss[range.min()].count(Name_ND{var, var2})) {
                            plot_ss[range.min()][Name_ND{var, var2}] = CreatePlot(("JSD_"+var+"_"+var2+"_Range"+std::to_string(range.min())+"_"+std::to_string(range.max())).c_str(),
                                                                                  ("JSD_"+var+"_"+var2+"_Range"+std::to_string(range.min())+"_"+std::to_string(range.max())).c_str(),"mass", "JSD");
                            plot_ss[range.min()][Name_ND{var, var2}]->SetLineColor(4);
                            plot_ss[range.min()][Name_ND{var, var2}]->SetMarkerColor(7);
                            plot_ss[range.min()][Name_ND{var, var2}]->SetMarkerSize(1);
                            plot_ss[range.min()][Name_ND{var, var2}]->SetMarkerStyle(8);
                        }
                        sample_1.push_back(&samples_mass.at(range.min()).at(var2));
                        band_1.push_back(bandwidth.at(range.min()).at(Name_ND{var2}));
                        sample_2.push_back(&samples_mass.at(sample_mass.first).at(var2));
                        band_2.push_back( bandwidth.at(sample_mass.first).at(Name_ND{var2}));
                        JSDpair_range_ss_future[mass_pair][Name_ND{var, var2}] = run::async(stat_estimators::JensenShannonDivergence_ND<double>,
                                                                                            sample_1, sample_2, band_1, band_2);

                        sample_1.erase(sample_1.end() - 1);
                        sample_2.erase(sample_2.end() - 1);
                        band_1.erase(band_1.end() - 1);
                        band_2.erase(band_2.end() - 1);
                    }
                }
            }
            for(auto& mass_pair : JSDvar_range_ss_future) {
                for(auto& name : mass_pair.second) {
                    if(!name.second.valid())
                        throw exception("future not valid");
                    JSDvar_range_ss[mass_pair.first][name.first] = name.second.get();
                }
            }
            for(auto& mass_pair : JSDpair_range_ss_future) {
                for(auto& name : mass_pair.second) {
                    if(!name.second.valid())
                        throw exception("future not valid");
                    JSDpair_range_ss[mass_pair.first][name.first] = name.second.get();
                }
            }
        }

        auto directory_distributionselected = directory_jenshan_ss->GetDirectory("Distribution_RangeSelected");
        if (directory_distributionselected == nullptr) {
            directory_jenshan_ss->mkdir("Distribution_RangeSelected");
            directory_distributionselected = directory_jenshan_ss->GetDirectory("Distribution_RangeSelected");
        }
        std::cout<<"directory"<<std::endl;

        int k = 0;
        for(auto& mass_pair : JSDvar_range_ss) {
            if (mass_pair.first.first == mass_pair.first.second) k = 0;
            for(auto& name : mass_pair.second) {
                histo_distribution.at(mass_pair.first)->Fill(JSDvar_range_ss.at(mass_pair.first).at(name.first));
                plot_ss.at(mass_pair.first.first).at(name.first)->SetPoint(k, mass_pair.first.second, JSDvar_range_ss.at(mass_pair.first).at(name.first));
            }
            k++;
        }
        int y = 0;
        for(auto& mass_pair : JSDpair_range_ss) {
            if (mass_pair.first.first == mass_pair.first.second) y = 0;
            for(auto& name : mass_pair.second) {
                histo_distribution.at(mass_pair.first)->Fill(JSDpair_range_ss.at(mass_pair.first).at(name.first));
                plot_ss.at(mass_pair.first.first).at(name.first)->SetPoint(y, mass_pair.first.second, JSDpair_range_ss.at(mass_pair.first).at(name.first));
            }
            y++;
        }
        for (const auto& range : ranges){
            directory_plotselected->mkdir(("Range"+std::to_string(range.min())+"_"+std::to_string(range.max())).c_str());
            auto directory_plotrange = directory_plotselected->GetDirectory(("Range"+std::to_string(range.min())+"_"+std::to_string(range.max())).c_str());
            for (const auto& name : plot_ss.at(range.min())){
                root_ext::WriteObject(*name.second, directory_plotrange);
            }
        }
        for (const auto& mass_pair : histo_distribution){
            root_ext::WriteObject(*mass_pair.second, directory_distributionselected);
        }

        std::cout<<"unione sample"<<std::endl;
        MassNameElement bandwidth_range, mutual_matrix_range, correlation_matrix_range_signal, correlation_matrix_range_bkg;
        MassVar samples_range;
        for (const auto& range : ranges){
            for (const auto& var: range_selected.at(range.min())){
                for (const auto& sample_mass : samples_mass){
                    if (!range.Contains(sample_mass.first)) continue;
                    for (const auto& entry : sample_mass.second.at(var)){
                        samples_range[range.min()][var].push_back(entry);
                    }
                }
            }
        }
        outfile->mkdir("Correlation");
        auto directory_correlation = outfile->GetDirectory("Correlation");
        for (const auto& sample: samples_range){
            std::cout<<"----Range"<<sample.first<<"----"<<" entries: "<<sample.second.at("pt_l1").size()<<std::endl;
            std::cout<<"correlation  ";
            correlation_matrix_range_signal[sample.first] = Correlation(sample.second);
            correlation_matrix_range_bkg[sample.first] = CorrelationSelected(samples_mass.at(Bkg), range_selected.at(sample.first));
            int bin = range_selected.at(sample.first).size();
            auto matrix = std::make_shared<TH2D>(("Bkg_"+std::to_string(sample.first)).c_str(),("Bkg_"+std::to_string(sample.first)).c_str(),bin,0,bin,bin,0,bin);
            int i = 1;
            for(const auto& var_1 : range_selected.at(sample.first)) {
                int j = 1;
                matrix->GetXaxis()->SetBinLabel(i, (var_1).c_str());
                matrix->GetYaxis()->SetBinLabel(i, (var_1).c_str());
                for(const auto& var_2 : range_selected.at(sample.first)) {
                    matrix->SetBinContent(i, j, correlation_matrix_range_bkg.at(sample.first).at(Name_ND{var_1, var_2})*100);
                    j++;
                }
                i++;
            }
            root_ext::WriteObject(*matrix, directory_correlation);
            std::cout<<"bandwidth  ";
            bandwidth_range[sample.first] = OptimalBandwidth(sample.second);
            std::cout<<"mutual  ";
            mutual_matrix_range[sample.first] = Mutual(sample.second, bandwidth.at(sample.first));
            if (sample.first == Bkg){
                std::cout<<std::endl;
                continue;
            }
            else {
                std::cout<<"mutual plot  ";
                MutualHisto(sample.first, mutual_matrix.at(sample.first), mutual_matrix.at(Bkg), directory_mutinf);
            }
            std::cout<<std::endl;

        }
        std::cout <<"Mutual histos" << std::endl;
        directory_mutinf->mkdir("Matrix");
        auto directory_mutinf_matrix = directory_mutinf->GetDirectory("Matrix");
        CreateMatrixHistos(samples_range, mutual_matrix_range, "MI", directory_mutinf_matrix);
        std::cout <<"Correlation histos" << std::endl;

        CreateMatrixHistos(samples_range, correlation_matrix_range_signal, "correlation_signal", directory_correlation);
//        CreateMatrixHistos(correlation_matrix_range_bkg, "correlation_background", directory_correlation);

        stop = clock::now();
        std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;


        std::map<int, std::map<Name_ND, double>> JSDvar_range_sb, JSDpair_range_sb;
        std::map<Name_ND, std::shared_ptr<TGraph>> plot_sb;
        std::map<int, std::shared_ptr<TH1D>> histo_distribution_sb;
        for (const auto& range : ranges){
            std::map<int, std::map<Name_ND, std::future<double>>> JSDvar_range_sb_future, JSDpair_range_sb_future ;
            for (const auto& sample_mass: samples_range){
                if (!range.Contains(sample_mass.first)) continue;
                histo_distribution_sb[range.min()] = std::make_shared<TH1D>(("JSDrange_"+std::to_string(range.min())+"_SignalBkg").c_str(),
                                                                          ("JSDrange_"+std::to_string(range.min())+"_SignalBkg").c_str(),
                                                                          50,0,1);
                histo_distribution_sb[range.min()]->SetXTitle("JSD");
                for (const auto& var: range_selected.at(range.min())){
                    if (!plot_sb.count(Name_ND{var})){
                        plot_sb[Name_ND{var}] = CreatePlot(("JSD_"+var+"_SignalBkg").c_str(), ("JSD_"+var+"_SignalBkg").c_str(), "mass", "JSD");
                        plot_sb[Name_ND{var}]->SetLineColor(kRed+2);
                        plot_sb[Name_ND{var}]->SetMarkerColor(2);
                        plot_sb[Name_ND{var}]->SetMarkerSize(1);
                        plot_sb[Name_ND{var}]->SetMarkerStyle(5);
                    }
                    std::vector<const DataVector*> sample_1, sample_2;
                    DataVector band_1, band_2;
                    sample_1.push_back(&sample_mass.second.at(var));
                    band_1.push_back(bandwidth_range.at(range.min()).at(Name_ND{var}));
                    sample_2.push_back(&samples_mass.at(Bkg).at(var));
                    band_2.push_back(bandwidth.at(Bkg).at(Name_ND{var}));
                    JSDvar_range_sb_future[range.min()][Name_ND{var}]  = run::async(stat_estimators::JensenShannonDivergence_ND<double>,
                                                                          sample_1, sample_2, band_1, band_2);
                    for (const auto& var2: range_selected.at(range.min())){
                        if (var2 <= var) continue;
                        if (!plot_sb.count(Name_ND{var, var2})) {
                            plot_sb[Name_ND{var, var2}] = CreatePlot(("JSD_"+var+"_"+var2+"_SignalBkg").c_str(),
                                                                     ("JSD_"+var+"_"+var2+"_SignalBkg").c_str(),"mass", "JSD");
                            plot_sb[Name_ND{var, var2}]->SetLineColor(kRed+2);
                            plot_sb[Name_ND{var, var2}]->SetMarkerColor(2);
                            plot_sb[Name_ND{var, var2}]->SetMarkerSize(1);
                            plot_sb[Name_ND{var, var2}]->SetMarkerStyle(5);
                        }
                        sample_1.push_back(&sample_mass.second.at(var2));
                        band_1.push_back(bandwidth_range.at(range.min()).at(Name_ND{var2}));
                        sample_2.push_back(&samples_mass.at(Bkg).at(var2));
                        band_2.push_back(bandwidth.at(Bkg).at(Name_ND{var2}));
                        JSDpair_range_sb_future[range.min()][Name_ND{var, var2}] = run::async(stat_estimators::JensenShannonDivergence_ND<double>,
                                                                                            sample_1, sample_2, band_1, band_2);

                        sample_1.erase(sample_1.end() - 1);
                        sample_2.erase(sample_2.end() - 1);
                        band_1.erase(band_1.end() - 1);
                        band_2.erase(band_2.end() - 1);
                    }
                }
            }
            for(auto& mass_pair : JSDvar_range_sb_future) {
                for(auto& name : mass_pair.second) {
                    if(!name.second.valid())
                        throw exception("future not valid");
                    JSDvar_range_sb[mass_pair.first][name.first] = name.second.get();
                }
            }
            for(auto& mass_pair : JSDpair_range_sb_future) {
                for(auto& name : mass_pair.second) {
                    if(!name.second.valid())
                        throw exception("future not valid");
                    JSDpair_range_sb[mass_pair.first][name.first] = name.second.get();
                }
            }
        }
        directory_jensenshannon->mkdir("Range_SignalBkg");
        auto directory_rangebkg = directory_jensenshannon->GetDirectory("Range_SignalBkg");
        auto directory_plotselected_sb = directory_rangebkg->GetDirectory("Plot_RangeSelected");
        if (directory_plotselected_sb == nullptr) {
            directory_rangebkg->mkdir("Plot_RangeSelected");
            directory_plotselected_sb = directory_rangebkg->GetDirectory("Plot_RangeSelected");
        }
        auto directory_distributionselected_sb = directory_rangebkg->GetDirectory("Distribution_RangeSelected");
        if (directory_distributionselected_sb == nullptr) {
            directory_rangebkg->mkdir("Distribution_RangeSelected");
            directory_distributionselected_sb = directory_rangebkg->GetDirectory("Distribution_RangeSelected");
        }
        int w = 0;
        for(auto& mass : JSDvar_range_sb) {
            for(auto& name : mass.second) {
                histo_distribution_sb.at(mass.first)->Fill(JSDvar_range_sb.at(mass.first).at(name.first));
                plot_sb.at(name.first)->SetPoint(w, mass.first, JSDvar_range_sb.at(mass.first).at(name.first));
            }
            w++;
        }
        int j = 0;
        for(auto& mass : JSDpair_range_sb) {
            for(auto& name : mass.second) {
                histo_distribution_sb.at(mass.first)->Fill(JSDpair_range_sb.at(mass.first).at(name.first));
                plot_sb.at(name.first)->SetPoint(j, mass.first, JSDpair_range_sb.at(mass.first).at(name.first));
            }
            j++;
        }
        for (const auto& name : plot_sb){
            root_ext::WriteObject(*name.second, directory_plotselected_sb);
        }
        for (const auto& mass_pair : histo_distribution_sb){
            root_ext::WriteObject(*mass_pair.second, directory_distributionselected_sb);
        }

        start = clock::now();
        std::cout<<"Plot matrix";
        directory_rangebkg->mkdir("Matrix");
        auto directory_matrix = directory_rangebkg->GetDirectory("Matrix");
        for (const auto& range: ranges){
            int bin = range_selected.at(range.min()).size();
            auto matrix_jsd_sb = std::make_shared<TH2D>(("JSD_Signal_Bkg_Range"+std::to_string(range.min())+"_"+std::to_string(range.min())).c_str(),("JSD_Signal_Bkg_Range"+std::to_string(range.min())+"_"+std::to_string(range.min())).c_str(), bin, 0, bin, bin, 0, bin);
            int i = 1;
            for (const auto& var1: range_selected.at(range.min())){
                    matrix_jsd_sb->GetXaxis()->SetBinLabel(i, (var1).c_str());
                    matrix_jsd_sb->GetYaxis()->SetBinLabel(i, (var1).c_str());
                    int jj = 1;
                    for (const auto& var2: range_selected.at(range.min())){
                        if (var1 == var2){
                            matrix_jsd_sb->SetBinContent(i, i, JSDvar_range_sb.at(range.min()).at(Name_ND{var1}));
                        }
                        if (jj > i){
                            matrix_jsd_sb->SetBinContent(i, jj, JSDpair_range_sb.at(range.min()).at(Name_ND{var1, var2}));
                            matrix_jsd_sb->SetBinContent(jj, i, JSDpair_range_sb.at(range.min()).at(Name_ND{var1, var2}));
                        }
                        jj++;
                    }
                    i++;
            }
            root_ext::WriteObject(*matrix_jsd_sb, directory_matrix);
        }

        directory_jensenshannon->mkdir("Comparison");
        auto directory_compare = directory_jensenshannon->GetDirectory("Comparison");
        std::map<Name_ND, TCanvas*> canvas;
        for (const auto & range: ranges){
            for (const auto& selected_name : plot_ss.at(range.min())){
                if (canvas.count(selected_name.first)) continue;
                if (selected_name.first.size() == 1)
                    canvas[selected_name.first] = new TCanvas((*selected_name.first.names.begin()).c_str(), (*selected_name.first.names.begin()).c_str(), 10,10,800,800);
                else
                    canvas[selected_name.first] = new TCanvas((*selected_name.first.names.begin()+"_"+*selected_name.first.names.rbegin()).c_str(), (*selected_name.first.names.begin()+"_"+*selected_name.first.names.rbegin()).c_str(), 10,10,800,800);
                canvas[selected_name.first]->DrawFrame(200,0,1000,1);
                canvas[selected_name.first]->SetGrid();
            }
        }
        for (const auto& selected_name : plot_sb){
            canvas.at(selected_name.first)->cd();
            selected_name.second->Draw("P");
            for (const auto & range: ranges){
                if (!plot_ss.at(range.min()).count(selected_name.first)) continue;
                plot_ss.at(range.min()).at(selected_name.first)->Draw("same");
            }
            plot_jsd.at(selected_name.first)->Draw("same");
            root_ext::WriteObject<TCanvas>(*canvas[selected_name.first], directory_compare);
            canvas[selected_name.first]->Close();
        }
        auto stop_tot = clock::now();
        std::cout<<"secondi totali: "<<std::chrono::duration_cast<std::chrono::seconds>(stop_tot - start_tot).count()<<std::endl;
    }
private:
    Arguments args;
    SampleEntryCollection samples;
    std::shared_ptr<TFile> outfile;
    bool split_training_testing;
    MvaVariables::VarNameSet enabled_vars, disabled_vars;
    MvaVariablesStudy vars;
};
}
}

PROGRAM_MAIN(analysis::Mva_Study::VariableDistribution, Arguments) // definition of the main program function

