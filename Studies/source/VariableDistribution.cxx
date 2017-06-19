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
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "hh-bbtautau/Analysis/include/MvaMethods.h"
#include "hh-bbtautau/Analysis/include/MvaConfigurationReader.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(std::string, tree_name);
};

namespace analysis {
namespace mva_study{

using clock = std::chrono::system_clock;

class VariablesDistribution {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    SampleIdVarData sample_vars;

    VariablesDistribution(const Arguments& _args): args(_args),
        outfile(root_ext::CreateRootFile(args.output_file())), reporter(std::make_shared<TimeReporter>())
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

    void Histo1D(const VarData& sample, TDirectory* directory)
    {
        for (const auto & var : sample){
            auto histo_var = std::make_shared<TH1D>((var.first).c_str(), (var.first).c_str(), 50,0,0);
            histo_var->SetCanExtend(TH1::kAllAxes);
            histo_var->SetXTitle((var.first).c_str());
            for(const auto& entry : var.second ){
                    histo_var->Fill(entry);
            }
            root_ext::WriteObject(*histo_var, directory);
        }
    }

    void Histo2D(const VarData& sample, TDirectory* directory)
    {
        for(auto var_1 = sample.begin(); var_1 != sample.end(); ++var_1){
            for(auto var_2 = std::next(var_1); var_2 != sample.end(); ++var_2) {
                auto histo_pair = std::make_shared<TH2D>((var_1->first+"_"+var_2->first).c_str(), (var_1->first+"_"+var_2->first).c_str(), 50,0,0,50,0,0);
                histo_pair->SetCanExtend(TH1::kAllAxes);
                histo_pair->SetXTitle((var_1->first).c_str());
                histo_pair->SetYTitle((var_2->first).c_str());
                for(size_t i = 0; i < var_1->second.size(); i++){
                    if (var_1->second.size() != var_2->second.size())
                        throw analysis::exception("Different vector size.");
                    histo_pair->Fill(var_1->second[i], var_2->second[i]);
                }
                root_ext::WriteObject(*histo_pair, directory);
            }
        }
    }

    void LoadData()
    {
        for(const SampleEntry& entry:samples)
        {
            if ( entry.channel != "" && args.tree_name() != entry.channel ) continue;
            auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
            EventTuple tuple(args.tree_name(), input_file.get(), true, {} , GetMvaBranches());
            Long64_t tot_entries = 0;
            for (Long64_t  current_entry = 0; current_entry < tuple.GetEntries(); current_entry++) {
                tuple.GetEntry(current_entry);
                const Event& event = tuple.data();

                if (static_cast<EventEnergyScale>(event.eventEnergyScale) != EventEnergyScale::Central || (event.q_1+event.q_2) != 0 || event.jets_p4.size() < 2
                    || event.extraelec_veto == true || event.extramuon_veto == true || event.jets_p4[0].eta() > cuts::btag_2016::eta
                    || event.jets_p4[1].eta() > cuts::btag_2016::eta)
                    continue;

                LorentzVectorE_Float bb = event.jets_p4[0] + event.jets_p4[1];

                if (!cuts::hh_bbtautau_2016::hh_tag::IsInsideEllipse(event.SVfit_p4.mass(), bb.mass()))
                    continue;

                tot_entries++;
                vars.AddEvent(event, entry.id, entry.weight);
            }
            std::cout << entry << " number of events: " << tot_entries << std::endl;
        }
        sample_vars = vars.GetSampleVariables();
        TimeReport();
    }


    void Run()
    {
        LoadData();

        auto directory_distribution1d = root_ext::GetDirectory(*outfile, "Distribution1D");
        auto directory_distribution2d = root_ext::GetDirectory(*outfile, "Distribution2D");
        for (const auto& sample: sample_vars){
            std::string mass = ToString(sample.first);
            std::cout<<"-----"<<mass<<"-----"<<std::endl;

            std::cout<<"Distribution1D"<<std::endl;
            auto directory_mass1d = root_ext::GetDirectory(*directory_distribution1d, mass);
            Histo1D(sample.second, directory_mass1d);

            std::cout<<"Distribution2D"<<std::endl;
            auto directory_mass2d = root_ext::GetDirectory(*directory_distribution2d, (mass).c_str());
            Histo2D(sample.second, directory_mass2d);
            TimeReport();
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

PROGRAM_MAIN(analysis::mva_study::VariablesDistribution, Arguments) // definition of the main program function

//./run.sh VariableDistribution --input_path ~/Desktop/tuples --output_file VariableDistribution_muTau.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name muTau
