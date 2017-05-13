/*! Variable distribution(1D and 2D) for all the variables and all the pairs for each mass
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <random>
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/StatEstimators.h"
#include "hh-bbtautau/Analysis/include/MvaVariables.h"
#include "hh-bbtautau/Analysis/include/MvaConfiguration.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(std::string, tree_name);
    OPT_ARG(bool, split, false);
    OPT_ARG(bool, istraining, false);
};

namespace analysis {

using clock = std::chrono::system_clock;

using DataVector = std::vector<double>;
using VarData = std::map<std::string, DataVector>;
using MassVar = std::map<int, VarData>;

void Histo1D(const VarData& sample, TDirectory* directory){
    for (const auto & var : sample){
        auto histo = std::make_shared<TH1D>((var.first).c_str(), (var.first).c_str(), 50,0,0);
        histo->SetCanExtend(TH1::kAllAxes);
        histo->SetXTitle((var.first).c_str());
        for(const auto& entry : var.second ){
                histo->Fill(entry);
        }
        root_ext::WriteObject(*histo, directory);
    }
}

void Histo2D(const VarData& sample, TDirectory* directory){
    for (const auto & var_1 : sample){
        for (const auto & var_2 : sample){
            auto histo = std::make_shared<TH2D>((var_1.first+"_"+var_2.first).c_str(), (var_1.first+"_"+var_2.first).c_str(), 50,0,0,50,0,0);
            histo->SetCanExtend(TH1::kAllAxes);
            histo->SetXTitle((var_1.first).c_str());
            histo->SetYTitle((var_2.first).c_str());
            for(size_t i = 0; i < var_1.second.size(); i++){
                if (var_1.second.size() != var_2.second.size())
                    throw analysis::exception("Different vector size.");
                histo->Fill(var_1.second[i],var_2.second[i]);
            }
            root_ext::WriteObject(*histo, directory);
        }
    }
}

class VariablesDistribution {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    static const std::set<std::string>& GetEnabledBranches()
    {
        static const std::set<std::string> EnabledBranches_read = {
            "eventEnergyScale", "q_1", "q_2", "jets_p4", "extraelec_veto", "extramuon_veto ", "SVfit_p4",
            "pfMET_p4", "p4_1", "p4_2"
        };
        return EnabledBranches_read;
    }

    VariablesDistribution(const Arguments& _args): args(_args), samples(ReadConfig(args.cfg_file())),
        outfile(root_ext::CreateRootFile(args.output_file())), split_training_testing(args.split()),
        vars(split_training_testing, UINT_FAST32_MAX, enabled_vars, disabled_vars)
    {
    }
    void Run()
    {
        auto start_tot = clock::now();

        MassVar sample_vars;
        auto start = clock::now();
        for(const SampleEntry& entry:samples)
        {
            if ( entry.channel != "" && args.tree_name() != entry.channel ) continue;
            auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
            EventTuple tuple(args.tree_name(), input_file.get(), true, {} , GetEnabledBranches());
            Long64_t tot_entries = 0;
            for (Long64_t  current_entry = 0; current_entry < tuple.GetEntries(); current_entry++) {
                tuple.GetEntry(current_entry);
                const Event& event = tuple.data();

                if (event.eventEnergyScale != 0 || (event.q_1+event.q_2) != 0 || event.jets_p4.size() < 2
                    || event.extraelec_veto == true || event.extramuon_veto == true || event.jets_p4[0].eta() > cuts::btag_2016::eta
                    || event.jets_p4[1].eta() > cuts::btag_2016::eta)
                    continue;

                LorentzVectorE_Float bb = event.jets_p4[0] + event.jets_p4[1];
                double ellipse_cut = pow(event.SVfit_p4.mass()-116,2)/pow(35.,2) + pow(bb.mass()-111,2)/pow(45.,2);
                if (ellipse_cut>1)
                    continue;

                tot_entries++;
                vars.AddEvent(event, entry.mass, entry.weight);
            }
            std::cout << entry << " number of events: " << tot_entries << std::endl;
        }
        auto stop = clock::now();
        std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;

        sample_vars = vars.GetSampleVariables(args.istraining());

        outfile->mkdir("Distribution1D");
        auto directory_distribution1d = outfile->GetDirectory("Distribution1D");
        outfile->mkdir("Distribution2D");
        auto directory_distribution2d = outfile->GetDirectory("Distribution2D");
        for (const auto& sample: sample_vars){
            start = clock::now();
            std::string mass;
            if (sample.first == Bkg) mass = "Bkg";
            else if (sample.first == Signal_SM) mass = "SM";
            else mass = std::to_string(sample.first);
            std::cout<<"-----"<<mass<<"-----"<<std::endl;
            std::cout<<"Distribution1D"<<std::endl;
            directory_distribution1d->mkdir((mass).c_str());
            auto directory_mass1d = directory_distribution1d->GetDirectory((mass).c_str());
            Histo1D(sample.second, directory_mass1d);
            stop = clock::now();
            std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;
            start = clock::now();
            std::cout<<"Distribution2D"<<std::endl;
            directory_distribution2d->mkdir((mass).c_str());
            auto directory_mass2d = directory_distribution2d->GetDirectory((mass).c_str());
            Histo2D(sample.second, directory_mass2d);
            stop = clock::now();
            std::cout<<"secondi: "<<std::chrono::duration_cast<std::chrono::seconds>(stop - start).count()<<std::endl;
        }

        auto stop_tot = clock::now();
        std::cout<<"secondi totali: "<<std::chrono::duration_cast<std::chrono::seconds>(stop_tot - start_tot).count()<<std::endl;
    }
private:
    Arguments args;
    SampleEntryCollection samples;
    std::shared_ptr<TFile> outfile;
    bool split_training_testing;
    std::mt19937 gen;
    std::uniform_int_distribution<> test_vs_training;
    MvaVariables::VarNameSet enabled_vars, disabled_vars;
    MvaVariablesStudy vars;
};
}

PROGRAM_MAIN(analysis::VariablesDistribution, Arguments) // definition of the main program function

//./run.sh VariableDistribution --input_path ~/Desktop/tuples --output_file VariableDistribution_muTau.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name muTau
