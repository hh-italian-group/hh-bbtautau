/*! Variable distribution(1D and 2D) for all the variables and all the pairs for each mass
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <random>
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/StatEstimators.h"
#include "hh-bbtautau/Analysis/include/MvaVariables.h"
#include "hh-bbtautau/Studies/include/MvaConfiguration.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"
#include "hh-bbtautau/Studies/include/MvaMethods.h"
#include "hh-bbtautau/Analysis/include/MvaConfigurationReader.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(int, spin);
    REQ_ARG(bool, skimmed);
    REQ_ARG(std::string, suffix);
};

namespace analysis {
namespace mva_study{

using clock = std::chrono::system_clock;

class VariablesDistribution {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

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

        samples = samples_list.at("Samples").files;
    }

    void Histo1D(const VarData& sample, TDirectory* directory)
    {
        for (const auto & var : sample){
            auto histo_var = std::make_shared<TH1D>((var.first).c_str(), (var.first).c_str(), 25,0,200);
            histo_var->SetCanExtend(TH1::kAllAxes);
            histo_var->SetXTitle((var.first).c_str());
            for(const auto& entry : var.second ){
                    histo_var->Fill(entry);
            }
            histo_var->Scale(1/histo_var->Integral());
            root_ext::WriteObject(*histo_var, directory);
        }
    }

    void Histo2D(const VarData& sample, TDirectory* directory)
    {
        for(auto var_1 = sample.begin(); var_1 != sample.end(); ++var_1){
            for(auto var_2 = var_1; var_2 != sample.end(); ++var_2) {
                auto histo_pair = std::make_shared<TH2D>((var_1->first+"_"+var_2->first).c_str(), (var_1->first+"_"+var_2->first).c_str(), 60,0,300,60,0,300);
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

    void LoadSkimmedData()
    {
        for(const SampleEntry& entry:samples)
        {
            if ( entry.id.IsSignal() && entry.spin != args.spin()) continue;
            auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
            auto tuple = ntuple::CreateEventTuple(args.tree_name(), input_file.get(), true, ntuple::TreeState::Skimmed);
            for(const Event& event : *tuple) {
                LorentzVectorE_Float bb = event.jets_p4[0] + event.jets_p4[1];
                if (args.suffix() == "_ANcut"){
                    if (!cuts::hh_bbtautau_2016::hh_tag::m_hh_window().IsInside(event.SVfit_p4.mass(),bb.mass())) continue;
                }
                auto eventInfoPtr =  analysis::MakeEventInfo(Parse<Channel>(args.tree_name()) ,event) ;
                EventInfoBase& eventbase = *eventInfoPtr;
                if (args.suffix() == "_newcut"){
                    if (!cuts::hh_bbtautau_2016::hh_tag::new_m_hh_window().IsInside(eventbase.GetHiggsTTMomentum(false).M(),bb.mass())) continue;
                }
                vars.AddEvent(eventbase, entry.id, entry.spin, entry.weight);
            }
            std::cout << entry << " number of events: " << tuple->size() << std::endl;
        }
        sample_vars = vars.GetSampleVariables(args.tree_name(), args.spin());
        TimeReport();
    }

    void LoadData()
    {
        for(const SampleEntry& entry:samples)
        {
            if ( entry.id.IsSignal() && entry.spin != args.spin()) continue;
            auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
            auto tuple = ntuple::CreateEventTuple(args.tree_name(), input_file.get(), true, ntuple::TreeState::Full);
            for(const Event& event : *tuple) {

                if (static_cast<EventEnergyScale>(event.eventEnergyScale) != EventEnergyScale::Central || (event.q_1+event.q_2) != 0 || event.jets_p4.size() < 2
                    || event.extraelec_veto == true || event.extramuon_veto == true || event.jets_p4[0].eta() > cuts::btag_2016::eta
                    || event.jets_p4[1].eta() > cuts::btag_2016::eta)
                    continue;

                LorentzVectorE_Float bb = event.jets_p4[0] + event.jets_p4[1];
                if (!cuts::hh_bbtautau_2016::hh_tag::m_hh_window().IsInside(event.SVfit_p4.mass(),bb.mass())) continue;
                if (entry.id == SampleType::Bkg_TTbar && event.file_desc_id>=2) continue;
                if (entry.id == SampleType::Sgn_NonRes && event.file_desc_id!=0) continue;
                auto eventInfoPtr =  analysis::MakeEventInfo(Parse<Channel>(args.tree_name()) ,event) ;
                EventInfoBase& eventbase = *eventInfoPtr;
                vars.AddEvent(eventbase, entry.id, entry.spin, entry.weight);
            }
            std::cout << entry << " number of events: " << tuple->size() << "  spin:" << entry.spin << "    " << entry.weight << std::endl;
        }
        sample_vars = vars.GetSampleVariables(args.tree_name(), args.spin());
        TimeReport();
    }

    void TimeReport(bool tot = false) const
    {
        reporter->TimeReport(tot);
    }

    void Run()
    {
        if (args.skimmed())
            LoadSkimmedData();
        else LoadData();

        auto directory_distribution1d = root_ext::GetDirectory(*outfile, "Distribution1D");
        auto directory_distribution2d = root_ext::GetDirectory(*outfile, "Distribution2D");

        for (const auto& sample: sample_vars){

            std::string mass = ToString(sample.first);
            std::cout<<"-----"<<mass<<"-----"<<std::endl;

            std::cout<<"Distribution1D"<<std::endl;
            auto directory_mass1d = root_ext::GetDirectory(*directory_distribution1d, mass.c_str());
            Histo1D(sample.second, directory_mass1d);

            std::cout<<"Distribution2D"<<std::endl;
            auto directory_mass2d = root_ext::GetDirectory(*directory_distribution2d, mass.c_str());
            Histo2D(sample.second, directory_mass2d);
            TimeReport();
        }

        TimeReport(true);
    }

private:
    Arguments args;
    SampleEntryCollection samples;
    std::shared_ptr<TFile> outfile;
    MvaVariablesStudy vars;
    std::shared_ptr<TimeReporter> reporter;
    SampleIdVarData sample_vars;

};
}
}

PROGRAM_MAIN(analysis::mva_study::VariablesDistribution, Arguments) // definition of the main program function
