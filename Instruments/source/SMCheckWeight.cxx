/*! Check SM weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "Instruments/include/SMFileConfigEntryReader.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "Instruments/include/SampleDescriptor.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "McCorrections/include/HH_BMStoSM_weight.h"


struct Arguments {
    run::Argument<std::string> weight_file{"weight_file", "file to reweight BSM samples"};
    run::Argument<std::string> input_path{"input_path", "Input path of the samples"};
    run::Argument<std::string> file_cfg_name{"file_cfg_name", "SM file cfg"};
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> output_file{"output_file", "Output root file"};
};

namespace analysis {

namespace sample_merging{

class SMCheckWeightData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    TH1D_ENTRY(m_sv, 20, 0, 400)
    TH1D_ENTRY(m_bb, 20, 0, 400)
    TH1D_ENTRY(m_kinFit, 25, 200, 1200)
};


class SMCheckWeight {
public:
    using HH_BMStoSM_weightPtr = std::shared_ptr<analysis::mc_corrections::HH_BMStoSM_weight>;

    SMCheckWeight(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.output_file())), anaData(output)
    {
        LoadInputs();
    }

    void Run() {}

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    SMCheckWeightData anaData;
    HH_BMStoSM_weightPtr smWeight;

    static const std::set<std::string>& GetEnabledBranches()
    {
        static const std::set<std::string> EnabledBranches_read = {
            "eventEnergyScale", "lhe_hh_m", "lhe_hh_cosTheta", "jets_p4", "SVfit_p4", "kinFit_m"
        };
        return EnabledBranches_read;
    }

    void LoadInputs()
    {

        smWeight = HH_BMStoSM_weightPtr(new analysis::mc_corrections::HH_BMStoSM_weight(args.weight_file(),"weight_node_BSM"));

        analysis::ConfigReader config_reader;

        SMBinDescriptorCollection file_descriptors;
        SMFileConfigEntryReader file_entry_reader(file_descriptors);
        config_reader.AddEntryReader("FILE", file_entry_reader, true);

        config_reader.ReadConfig(args.file_cfg_name());
        std::string name_sm;

        for (auto file_descriptor : file_descriptors){ //loop n file_descriptors
            const SMBinDescriptor file_descriptor_element = file_descriptor.second;
            const std::string& name = file_descriptor_element.name;

            for (auto single_file_path : file_descriptor_element.file_paths){ //loop on files

                std::cout << "File descriptor characteristics: " << file_descriptor.first << ", " <<
                             single_file_path << ", " << file_descriptor_element.fileType
                          << std::endl;
                std::string filename = args.input_path()  + "/" + single_file_path;
                auto inputFile = root_ext::OpenRootFile(filename);
                ntuple::EventTuple eventTuple(args.tree_name(), inputFile.get(), true, {}, GetEnabledBranches());
                const Long64_t n_entries = eventTuple.GetEntries();

                for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on entries
                    eventTuple.GetEntry(current_entry);
                    const ntuple::Event& event = eventTuple.data();
                    if (static_cast<EventEnergyScale>(event.eventEnergyScale) != analysis::EventEnergyScale::Central ||
                            event.jets_p4.size() < 2 /*|| (event.q_1+event.q_2)!=0*/)
                        continue;

                    LorentzVectorE_Float jet_1 = event.jets_p4.at(0);
                    LorentzVectorE_Float jet_2 = event.jets_p4.at(1);
                    if (!(std::abs(jet_1.eta()) < 2.1 && std::abs(jet_2.eta()) < 2.1)) continue;

                    LorentzVectorE_Float bb = jet_1 + jet_2;

//                    double elliptical_mass_cut=(pow(event.SVfit_p4.mass()-116.,2)/pow(35,2))+(pow(bb.M()-111,2)/pow(45,2));
//                    if (!(elliptical_mass_cut<1)) continue;

                    double sm_weight = smWeight->Get(event);

                    if (file_descriptor_element.fileType == FileType::sm){
                        name_sm = name;
                        sm_weight = 1;
                    }

                    anaData.m_bb(name).Fill(bb.M(),sm_weight);
                    anaData.m_sv(name).Fill(event.SVfit_p4.M(),sm_weight);
                    //for (unsigned n = 0; n < event.kinFit_m.size(); ++n)
                    anaData.m_kinFit(name).Fill(event.kinFit_m.at(0),sm_weight);

                } //end loop on entries
            } // end loop on files

//            const double scale_m_bb = 1 / anaData.m_bb(name).Integral(0,anaData.m_bb(name).GetNbinsX()+1);//PUT OVERFLOW - USE RENORMALIZE
//            anaData.m_bb(name).Scale(scale_m_bb);
//            const double scale_m_sv = 1 / anaData.m_sv(name).Integral(0,anaData.m_sv(name).GetNbinsX()+1);
//            anaData.m_sv(name).Scale(scale_m_sv);
//            const double scale_m_kinFit = 1 / anaData.m_kinFit(name).Integral(0,anaData.m_kinFit(name).GetNbinsX()+1);
//            anaData.m_kinFit(name).Scale(scale_m_kinFit);

            RenormalizeHistogram(anaData.m_bb(name), 1, true);
            RenormalizeHistogram(anaData.m_sv(name), 1, true);
            RenormalizeHistogram(anaData.m_kinFit(name), 1, true);

        } //end loop n file_descriptors

        for (const auto& file_descriptor : file_descriptors) {
            if (file_descriptor.second.fileType == FileType::sm) continue;
            const std::string& name = file_descriptor.second.name;

            std::cout << "KS m_bb: " << anaData.m_bb(name).KolmogorovTest(&anaData.m_bb(name_sm)) << std::endl;
            std::cout << "KS m_sv: " << anaData.m_sv(name).KolmogorovTest(&anaData.m_sv(name_sm)) << std::endl;
            std::cout << "KS m_kinFit: " << anaData.m_kinFit(name).KolmogorovTest(&anaData.m_kinFit(name_sm)) << std::endl;

        }

    }


};

} //namespace sample_merging

} //namespace analysis

PROGRAM_MAIN(analysis::sample_merging::SMCheckWeight, Arguments)
