/*! Study of correlation matrix, mutual information and Jensen Shannon Divergence to search range of masses.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <fstream>
#include <random>
#include <algorithm>
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

static constexpr double threashold_mi = 0.7;
VarNameSet CheckList(int mass,  MassVar sample, const NameElement& JSDivergenceND, const NameElement& mutual_matrix_signal,
                                const NameElement& mutual_matrix_bkg, std::map<std::string, std::string>& eliminated, const MassNameElement& band,
                                size_t number_variables, TDirectory* directory){

    VarNameSet selected, not_corrected;
    VectorName_ND JSDivergence_vector(JSDivergenceND.begin(), JSDivergenceND.end());

    while(selected.size() < number_variables && JSDivergence_vector.size()) {
        std::sort(JSDivergence_vector.begin(), JSDivergence_vector.end(),
                  [](const std::pair<Name_ND, double>& el1, const std::pair<Name_ND, double>& el2){
            return el1.second > el2.second;
        });
        auto best_entry = JSDivergence_vector.front();
        for(const auto& name : best_entry.first.names) {
            if(selected.count(name))
                continue;
            const double JS_1d = JSDivergenceND.at(Name_ND({name}));
            for(auto& entry : JSDivergence_vector) {
                if(entry.first.names.count(name))
                    entry.second -= JS_1d;
            }
            if (not_corrected.count(name))
                continue;

            for(const auto& other_entry : sample.at(mass)) {
                if(other_entry.first == name) continue;
                const Name_ND names({name, other_entry.first});
                if(mutual_matrix_signal.at(names) < threashold_mi && mutual_matrix_bkg.at(names) < threashold_mi){
                    eliminated[other_entry.first] = "MI-"+name;
                    not_corrected.insert(other_entry.first);
                }
            }
            selected.insert(name);
        }
        JSDivergence_vector.pop_front();
        VectorName_ND copy;
        for (const auto& entry : JSDivergence_vector){
            if(entry.second <= 0) continue;
            bool consider = true;
            for(const auto& other : entry.first) {
                consider = consider && !not_corrected.count(other);
            }
            if(consider)
                copy.push_back(entry);
        }
        JSDivergence_vector = std::move(copy);
    }
    auto directory_v = directory->GetDirectory("Plot_Selected_variables");
    if (directory_v == nullptr) {
        directory->mkdir("Plot_Selected_variables");
        directory_v = directory->GetDirectory("Plot_Selected_variables");
    }
    TDirectory* directory_m;
    if (mass==Signal_SM){
        directory_v->mkdir("Signal_SM");
        directory_m = directory_v->GetDirectory("Signal_SM");
    }
    else{
    directory_v->mkdir(("Signal_mass"+std::to_string(mass)).c_str());
    directory_m = directory_v->GetDirectory(("Signal_mass"+std::to_string(mass)).c_str());
    }

    std::map<std::string, std::map<int, std::future<double>>> JSDss_selected_future;
    std::map<std::string, std::map<int, double>> JSDss_selected;
    std::map<std::string, std::shared_ptr<TGraph>> plot;
    for(const auto& selected_name : selected){
        plot[selected_name] = CreatePlot(("JSD_"+selected_name).c_str(),("JSD_"+selected_name).c_str(),"mass","JSD");
        std::vector<const DataVector*> sample_1;
        DataVector band_1;
        sample_1.push_back(&sample.at(mass).at(selected_name));
        band_1.push_back(band.at(mass).at(Name_ND{selected_name}));
        for (const auto& mass_entry: sample){
            if (mass_entry.first == Bkg) continue;
            std::vector<const DataVector*> sample_2;
            DataVector band_2;
            sample_2.push_back(&sample.at(mass_entry.first).at(selected_name));
            band_2.push_back(band.at(mass_entry.first).at(Name_ND{selected_name}));
            JSDss_selected_future[selected_name][mass_entry.first] = run::async(stat_estimators::JensenShannonDivergence_ND<double>, sample_1, sample_2, band_1, band_2);
        }
    }
    for(auto& var : JSDss_selected_future) {
        int i = 0;
        for(auto& mass : var.second) {
            mass.second.wait();
            if(!mass.second.valid())
                throw exception("future not valid");
            JSDss_selected[var.first][mass.first] = mass.second.get();
            plot.at(var.first)->SetPoint(i,mass.first,JSDss_selected.at(var.first).at(mass.first));
            i++;
        }
        root_ext::WriteObject(*plot.at(var.first), directory_m);
    }
    return selected;
}

void PlotJensenShannon(const MassNameElement& JSDivergenceND, TDirectory* directory){
    std::map<std::string, std::shared_ptr<TGraph>> plot;
    int i = 0;
    for (const auto& mass_entry: JSDivergenceND){
        if (mass_entry.first == Bkg)
            continue;
        for (const auto& var : mass_entry.second){
            if (var.first.names.size()!=1)
                continue;
            const std::string name = *var.first.names.begin();
            if (plot.count(name) == 0) {
                plot[name] = CreatePlot("JSD_"+name+"_SignalBkg","JSD_"+name+"_SignalBkg","mass","JSD");
            }
            plot[name]->SetPoint(i,mass_entry.first,JSDivergenceND.at(mass_entry.first).at(Name_ND({name})));
        }
        i++;
    }
    for(const auto& var: plot){
        root_ext::WriteObject(*var.second, directory);
    }
}

void KolmogorovSignalCompatibility(int mass, const MassVar& sample_var, const VarNameSet& selected, TDirectory* directory){
    auto directory_v = directory->GetDirectory("Plot_Selected_variables");
    if (directory_v == nullptr) {
        directory->mkdir("Plot_Selected_variables");
        directory_v = directory->GetDirectory("Plot_Selected_variables");
    }
    TDirectory* directory_m;
    if (mass == Signal_SM){
        directory_v->mkdir("Signal_SM");
        directory_m = directory_v->GetDirectory("Signal_SM");
    }
    else{
    directory_v->mkdir(("Signal_mass"+std::to_string(mass)).c_str());
    directory_m = directory_v->GetDirectory(("Signal_mass"+std::to_string(mass)).c_str());
    }
    std::map<std::string, std::shared_ptr<TGraph>> plot;
    MassNameElement kolmogorov;
    int i = 0;
    for (const auto& mass_entry: sample_var){
        if ( mass_entry.first == Bkg )
            continue;
        for (const auto& var : selected){
            std::vector<double> vector_signal = sample_var.at(mass_entry.first).at(var);
            std::sort(vector_signal.begin(), vector_signal.end());\
            std::vector<double> vector_signal_2 = sample_var.at(mass).at(var);
            std::sort(vector_signal_2.begin(), vector_signal_2.end());
            Double_t* v_s = vector_signal.data(), *v_s_2 = vector_signal_2.data();
            kolmogorov[mass_entry.first][Name_ND{var}] = TMath::KolmogorovTest(vector_signal.size(), v_s, vector_signal_2.size(), v_s_2, "");
            if (plot.count(var) == 0) {
                plot[var] = CreatePlot(("KolmogorovSmirnov_"+var).c_str(), ("KolmogorovSmirnov_"+var).c_str(), "mass","KS Probability" );
            }
            plot[var]->SetPoint(i,mass_entry.first,kolmogorov.at(mass_entry.first).at(Name_ND{var}));
        }
        i++;
    }
    for(const auto& var: selected){
        root_ext::WriteObject(*plot[var], directory_m);
    }
}

using VecVariables =  std::vector<std::string>;
void VariablesSelection(const MassVar& sample, const MassNameElement& JSDivergenceSB, const MassNameElement& bandwidth, const MassNameElement& mutual_matrix,
                   std::shared_ptr<TFile> outfile, size_t number_variables){

    outfile->mkdir("Kolmogorov");
    auto directory_ks = outfile->GetDirectory("Kolmogorov");
    auto directory_jensenshannon = outfile->GetDirectory("JensenShannonDivergence");
    auto directory_jenshan_sgnlbkg = directory_jensenshannon->GetDirectory("Signal_Background");
    std::cout<<"Checklist1"<<std::endl;
    directory_jenshan_sgnlbkg->mkdir("Distribution");
    auto directory_jenshan_distrib = directory_jenshan_sgnlbkg->GetDirectory("Distribution");
    directory_jensenshannon->mkdir("Signal_Signal");
    auto directory_jenshan_sgnlsgnl = directory_jensenshannon->GetDirectory("Signal_Signal");
    directory_jenshan_sgnlsgnl->mkdir("Check_Mass");
    auto directory_jenshan_mass = directory_jenshan_sgnlsgnl->GetDirectory("Check_Mass");
    MassVarNameSet mass_selected;
    std::ofstream ofs("SelectedVariable.csv", std::ofstream::out);
    for (const auto& mass_entry : sample){
        ofs << mass_entry.first << ",";
        std::cout<<mass_entry.first<<" ";
        std::cout.flush();
        if (mass_entry.first == Bkg)
            continue;
        std::map<std::string, std::string> eliminated;
        mass_selected[mass_entry.first] =  CheckList(mass_entry.first, sample,JSDivergenceSB.at(mass_entry.first),
                                                     mutual_matrix.at(mass_entry.first), mutual_matrix.at(Bkg),
                                                     eliminated, bandwidth, number_variables, directory_jenshan_mass);

        KolmogorovSignalCompatibility(mass_entry.first, sample, mass_selected[mass_entry.first], directory_ks);
        ofs << mass_entry.first << std::endl;
        for(auto& entry_selected: mass_selected[mass_entry.first]){
                ofs << entry_selected<<",";
        }
        ofs << std::endl;
        std::string  mass = std::to_string(mass_entry.first);

        if (mass_entry.first == Signal_SM) mass = "SM";
        auto histo = std::make_shared<TH1D>(("JSD_Signal"+mass+"_Background").c_str(), ("JensenShannonDivergence_Signal"+mass+"_Background").c_str(), 50,0,1);
        histo->SetXTitle("JSD");
        std::ofstream of("InformationTable_mass"+(std::to_string(mass_entry.first))+".csv", std::ofstream::out);
        of<<"Var_1"<<","<<"Var_2"<<","<<"JSD_ND"<<","<<"JSD_12-(JSD_1+JSD_2)"<<","<<"ScaledMI_Signal_12" <<","<<"ScaledMI_Bkg_12"<<","<<"selected"<<","<<","<<"eliminated by"<<std::endl;
        for(auto& entry : JSDivergenceSB.at(mass_entry.first)){
             histo->Fill(entry.second);
             if ((entry.first).names.size() == 1){
                 auto name_1 = ((entry.first).names).begin();
                 bool selected_1 = mass_selected[mass_entry.first].count(*name_1);
                 of<<*name_1<<","<<"   "<<","<<JSDivergenceSB.at(mass_entry.first).at(Name_ND({*name_1}))<<","<<","<<","<<","<<selected_1;
                 if(eliminated.count(*name_1)) of<<","<<","<<eliminated.at(*name_1);
                 of<<std::endl;
             }
             else  if ((entry.first).names.size() == 2){
                 auto name_1 = ((entry.first).names).begin();
                 auto name_2 = ((entry.first).names).rbegin();
                 bool selected_1 = mass_selected[mass_entry.first].count(*name_1);
                 bool selected_2 = mass_selected[mass_entry.first].count(*name_2);
                 of<<*name_1<<","<<*name_2<<","<<JSDivergenceSB.at(mass_entry.first).at(Name_ND({*name_1,*name_2}))<<","<<JSDivergenceSB.at(mass_entry.first).at(Name_ND({*name_1,*name_2}))-(JSDivergenceSB.at(mass_entry.first).at(Name_ND({*name_1}))+JSDivergenceSB.at(mass_entry.first).at(Name_ND({*name_2})))<<","<<mutual_matrix.at(mass_entry.first).at(Name_ND({*name_1,*name_2}))<<","<<mutual_matrix.at(Bkg).at(Name_ND({*name_1,*name_2}))<<","<<selected_1<<","<<selected_2;
                 if(eliminated.count(*name_1)) of<<","<<eliminated.at(*name_1);
                 if(eliminated.count(*name_2)) of<<","<<eliminated.at(*name_2);
                 of<<std::endl;
             }
        }
        of.close();
        root_ext::WriteObject(*histo, directory_jenshan_distrib);
    }
    std::cout<<std::endl<<"Intersection and Union"<<std::endl;
    int bin = mass_selected.size();
    auto matrix_intersection = std::make_shared<TH2D>("Intersection","Intersection of variables",bin, 0, bin, bin, 0, bin);
    matrix_intersection->SetXTitle("mass");
    matrix_intersection->SetYTitle("mass");
    int k = 1;
    for(const auto& mass_1: mass_selected){
        int j = 1;
        std::string mass;
        if (mass_1.first == Signal_SM) mass = "SM";
        else mass = std::to_string(mass_1.first);
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
    root_ext::WriteObject(*matrix_intersection, outfile.get());;
}


using SampleEntryCollection = std::vector<SampleEntry>;

class MvaClassification {
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

    MvaClassification(const Arguments& _args): args(_args), samples(ReadConfig(args.cfg_file())),
        outfile(root_ext::CreateRootFile(args.output_file())), split_training_testing(args.split()),
        vars(split_training_testing, UINT_FAST32_MAX, enabled_vars, disabled_vars)
    {
    }

    void Run()
    {
        auto start_tot = clock::now();
        run::ThreadPull threads(args.number_threads());

        MassVar sample_vars;
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

        sample_vars = vars.GetSampleVariables(args.istraining());

        std::cout << "n.variabili: " << sample_vars.at(Bkg).size() << std::endl;
        std::cout << "n.masse segnale: " << sample_vars.size() - 1 << " + " << 1 << " fondo."<<std::endl;

        MassNameElement bandwidth, mutual_matrix, correlation_matrix;

        std::cout << "Bandwidth and Mutual Information" << std::endl;
        outfile->mkdir("MutualInformation");
        auto directory_mutinf = outfile->GetDirectory("MutualInformation");
        start = clock::now();
        for (const auto& sample: sample_vars){
            std::cout<<"----"<<sample.first<<"----"<<" entries: "<<sample.second.at("pt_l1").size()<<std::endl;
            std::cout<<"correlation  ";
            correlation_matrix[sample.first] = Correlation(sample.second);
            std::cout<<"bandwidth  ";
            bandwidth[sample.first] = OptimalBandwidth(sample.second);
            std::cout<<"mutual  ";
            mutual_matrix[sample.first] = Mutual(sample.second, bandwidth.at(sample.first));
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
        CreateMatrixHistos(sample_vars, mutual_matrix, "MI", directory_mutinf_matrix);

        std::cout <<"Correlation histos" << std::endl;
        outfile->mkdir("Correlation");
        auto directory_correlation = outfile->GetDirectory("Correlation");
        CreateMatrixHistos(sample_vars, correlation_matrix, "correlation", directory_correlation);
        stop = clock::now();
        std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;

        outfile->mkdir("MI_Correlation");
        auto directory_mi_corr = outfile->GetDirectory("MI_Correlation");
        for(const auto& sample : sample_vars){
            std::string mass= "";
            if (sample.first == Signal_SM) mass = "_SM";
            else if (sample.first == Bkg) mass = "_Background";
            else mass = "_mass"+std::to_string(sample.first);
            auto histo = std::make_shared<TH2D>(("MI_Correlation_"+mass).c_str(),("MI_Correlation_"+mass).c_str(),50,0,1,50,0,1);
            histo->SetXTitle("MI");
            histo->SetYTitle("Correlation");
            for (const auto pair: mutual_matrix.at(sample.first)){
                histo->Fill(pair.second, std::abs(correlation_matrix.at(sample.first).at(pair.first)));
            }
            root_ext::WriteObject(*histo, directory_mi_corr);
            }

        MassNameElement JSDivergenceSB;
        std::cout<<"Jensen Shannon Signal Background";
        start = clock::now();
        for (const auto& mass_entry: sample_vars){
            if (mass_entry.first == Bkg)
                continue;
            JSDivergenceSB[mass_entry.first] = JensenDivergenceSB(mass_entry.second, sample_vars.at(Bkg),bandwidth.at(mass_entry.first), bandwidth.at(Bkg));
        }
        std::cout<<" - Plot Jensen Shannon"<<std::endl;
        outfile->mkdir("JensenShannonDivergence");
        auto directory_jensenshannon = outfile->GetDirectory("JensenShannonDivergence");
        directory_jensenshannon->mkdir("Signal_Background");
        auto directory_jenshan_sgnlbkg = directory_jensenshannon->GetDirectory("Signal_Background");
        directory_jenshan_sgnlbkg->mkdir("All_variables");
        auto directory_jenshan_allvars = directory_jenshan_sgnlbkg->GetDirectory("All_variables");
        PlotJensenShannon(JSDivergenceSB, directory_jenshan_allvars);
        directory_jenshan_sgnlbkg->mkdir("Matrix");
        auto directory_jenshan_matrix = directory_jenshan_sgnlbkg->GetDirectory("Matrix");
        CreateMatrixHistos(sample_vars, JSDivergenceSB, "JSD", directory_jenshan_matrix);
        stop = clock::now();
        std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;

        std::cout<<"Selection variables"<<std::endl;
        VariablesSelection(sample_vars, JSDivergenceSB, bandwidth, mutual_matrix, outfile, args.number_variables());
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

PROGRAM_MAIN(analysis::Mva_Study::MvaClassification, Arguments) // definition of the main program function

//./run.sh MvaRangesCompatibility --input_path ~/Desktop/tuples --output_file Checklist1_muTau.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name eTau  --number_threads 4 --number_variables 20 --number_events 50 --split true --istraining true
