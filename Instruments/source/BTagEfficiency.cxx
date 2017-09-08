// BTagEfficiency.cxx
#include <boost/format.hpp>
#include <vector>
#include <set>
#include <map>

#include "AnalysisTools/Run/include/program_main.h" // definition of wrappers for the program main and program arguments.
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"
#include "AnalysisTools/Core/include/Tools.h"
#include "AnalysisTools/Core/include/TextIO.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::vector<std::string>, input_file); // required argument "input_file"
    REQ_ARG(std::string, output_file); // required argument "output_file"
    OPT_ARG(bool, flag, false); // optional argument "flag" with the default value = false
};

class BTagData : public root_ext::AnalyzerData {
    public:
        explicit BTagData(std::shared_ptr<TFile> _outputFile, const std::string& directoryName = "") :
            AnalyzerData(_outputFile, directoryName)
            {

            }
        TH2D_ENTRY(h2,701,-0.5,700.5, 600,-3,3)
};

class BTagEfficiency { // simple analyzer definition
    public:
        using Event = ntuple::Event;
        using EventTuple = ntuple::EventTuple;
        using EventEnergyScale = analysis::EventEnergyScale;
        using Channel = analysis::Channel; 



        BTagEfficiency(const Arguments& _args) : args(_args), outfile(root_ext::CreateRootFile(args.output_file())), anaData(outfile)
    {
        // Analyzer initialization (e.g. open input/output files, parse configs...)
    }
        void Run()
        {
            // analyzer code
            static const std::set<Channel> channels = { Channel::ETau, Channel::MuTau };
            static const std::map<std::string, double> btag_working_points = { { "L", cuts::btag_2016::CSVv2L },
                { "M", cuts::btag_2016::CSVv2M } };
            static const std::map<int, std::string> flavours = { { 5, "b" }, { 4, "c" }, { 0, "udsg" } };
            static const std::string flavour_all = "all";
            static const std::string num = "Num", denom = "Denom", eff = "Eff" ;

            for(const auto& channel : channels) {
                const std::string treeName = ToString(channel);
                for (const auto& name : args.input_file()){
                    std::shared_ptr<TFile> in_file(root_ext::OpenRootFile(name));
                    std::shared_ptr<EventTuple> tuple;
                    try {
                        tuple = ntuple::CreateEventTuple(treeName, in_file.get(), true, ntuple::TreeState::Full,true);
                    } catch(std::exception&) {
                        std::cerr << "WARNING: tree "<<treeName<<" not found in file '"<< std::endl;
                    }

                    for(const Event& event : *tuple){

                        //auto event_ptr = std::make_shared<Event>(event);

                        const EventEnergyScale es = static_cast<EventEnergyScale>(event.eventEnergyScale);
                        if (es != EventEnergyScale::Central || event.jets_p4.size() < 2 || event.extraelec_veto
                                || event.extramuon_veto
                                || std::abs(event.jets_p4.at(0).eta()) >= cuts::btag_2016::eta
                                || std::abs(event.jets_p4.at(1).eta()) >= cuts::btag_2016::eta) continue;

                        auto bb = event.jets_p4.at(0) + event.jets_p4.at(1);
                        if (!cuts::hh_bbtautau_2016::hh_tag::IsInsideEllipse(event.SVfit_p4.mass(),bb.mass())) continue;

                        if ((event.q_1+event.q_2) != 0) continue;

                        static const std::string tau_iso_disc = "byMediumIsolationMVArun2v1DBoldDMwLT";
                        static const uint32_t tau_iso_disc_hash = analysis::tools::hash(tau_iso_disc);
                        static const double tau_iso_cut = 0.5;
                        if(!PassTauIdCut(event.tauId_keys_1, event.tauId_values_1, tau_iso_disc_hash, tau_iso_cut) ||
                                !PassTauIdCut(event.tauId_keys_2, event.tauId_values_2, tau_iso_disc_hash, tau_iso_cut)) continue;   

                        for (int i=0; i<2; i++){
                            const auto& jet = event.jets_p4.at(i);
                            double jet_csv = event.jets_csv.at(i);
                            double jet_hadronFlavour = event.jets_hadronFlavour.at(i);
                            const std::string& jet_flavour = flavours.at(jet_hadronFlavour);
                  
                            anaData.h2(denom, flavour_all, "all", channel).Fill(jet.Pt(), jet.Eta());
                            anaData.h2(denom, jet_flavour, "all", channel).Fill(jet.Pt(), jet.Eta());

                            for(const auto& btag_wp : btag_working_points) {
                                if(btag_wp.second > jet_csv){
                                    anaData.h2(num, flavour_all, btag_wp.first, channel).Fill(jet.pt(), jet.eta());
                                    anaData.h2(num, jet_flavour, btag_wp.first, channel).Fill(jet.pt(), jet.eta());
                                }
                            }
                        }//end loop on jets 
                    }//end loop on events
                }// end loop on files
                for(const auto& btag_wp : btag_working_points) {
                    for(const auto& hadron_flavour : flavours){
                        const std::string& jet_flavour = hadron_flavour.second;
                        anaData.h2(eff,jet_flavour,btag_wp.first,channel).CopyContent(anaData.h2(num,jet_flavour,btag_wp.first,channel));
                        anaData.h2(eff,jet_flavour,btag_wp.first,channel).Divide(&anaData.h2(denom,jet_flavour,"all",channel));
                    }
                    anaData.h2(eff,flavour_all,btag_wp.first,channel).CopyContent(anaData.h2(num,flavour_all,btag_wp.first,channel));
                    anaData.h2(eff,flavour_all,btag_wp.first,channel).Divide(&anaData.h2(denom,flavour_all,"all",channel));
                }
            }//end loop on channel
        }
    private:
        Arguments args;
        std::shared_ptr<TFile> outfile;
        BTagData anaData;



        static bool PassTauIdCut(const std::vector<uint32_t>& tauId_keys, const std::vector<float>& tauId_values,
                uint32_t discriminator_name_hash, float cut_value)
        {
            for(size_t n = 0; n < tauId_keys.size(); ++n) {
                if(tauId_keys.at( n ) == discriminator_name_hash) return tauId_values.at( n ) > cut_value;
            }
            return true;
        }



};


PROGRAM_MAIN(BTagEfficiency, Arguments) // definition of the main program function
