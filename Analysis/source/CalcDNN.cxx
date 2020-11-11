/*! Computes DNN scores for AnaTuple events.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "AnalysisTools/Core/include/ProgressReporter.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "cms_hh_tf_inference/inference/interface/inf_wrapper.hh"
#include "cms_hh_proc_interface/processing/interface/evt_proc.hh"
#include "cms_hh_proc_interface/processing/interface/feat_comp.hh"
#include "hh-bbtautau/Analysis/include/AnaTuple.h"
#include "hh-bbtautau/Analysis/include/AnaTupleReader.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"


namespace analysis {

struct Arguments {
    REQ_ARG(std::string, input);
    REQ_ARG(Channel, channel);
    REQ_ARG(Period, period);
    REQ_ARG(std::string, model);
    REQ_ARG(std::string, output);
    OPT_ARG(std::string, kls, "");
    OPT_ARG(std::string, spins, "");
    OPT_ARG(std::string, masses, "");
    OPT_ARG(unsigned, n_threads, 1);
    OPT_ARG(Long64_t, max_events, std::numeric_limits<Long64_t>::max());
    OPT_ARG(Long64_t, begin_entry_index, 0);
    OPT_ARG(Long64_t, end_entry_index, std::numeric_limits<Long64_t>::max());
};

class CalcDNN {
public:
    CalcDNN(const Arguments& _args) :
        args(_args), inputFile(root_ext::OpenRootFile(args.input())),
        outputFile(root_ext::CreateRootFile(args.output(), ROOT::kLZMA, 9)),
        anaTuple(ToString(args.channel()), inputFile.get(), true)
    {
        ROOT::EnableThreadSafety();
        if(args.n_threads() > 1)
            ROOT::EnableImplicitMT(args.n_threads());

        if(!args.kls().empty()) {
            if(!args.spins().empty() || !args.masses().empty())
                throw exception("Incompatible arguments.");
            const auto kls = SplitValueListT<double>(args.kls(), false, ", \t");
            for(double kl : kls) {
                bbtautau::HyperPoint p;
                p.kl = kl;
                points.push_back(p);
            }
        } else if(!args.spins().empty() && !args.masses().empty()) {
            const auto spins = SplitValueListT<int>(args.spins(), false, ", \t");
            const auto masses = SplitValueListT<double>(args.masses(), false, ", \t");
            for(int spin : spins) {
                for(double mass : masses) {
                    bbtautau::HyperPoint p;
                    p.spin = spin;
                    p.mass = mass;
                    points.push_back(p);
                }
            }
        } else {
            throw exception("No hyperparameters were provided.");
        }

        const std::string feat_file = args.model() + "/features.txt";
        const std::vector<std::string> requested = GetRequestedVariables(feat_file);
        evt_proc = std::make_unique<EvtProc>(false, requested, true);
        const bool verbose = false;
        wrapper = std::make_unique<InfWrapper>(args.model() + "/ensemble", args.n_threads(), verbose);
    }

    void Run()
    {
        const std::string channelName = ToString(args.channel());

        outputFile->cd();
        auto outputTree = std::make_unique<TTree>(channelName.c_str(), channelName.c_str());

        std::vector<float> outputs(points.size());
        for(size_t point_index = 0; point_index < points.size(); ++point_index) {
            const std::string branchName = "dnn_score_" + points.at(point_index).ToString();
            outputTree->Branch(branchName.c_str(), &outputs.at(point_index), (branchName + "/F").c_str());
        }

        const Long64_t endEntry = std::min(anaTuple.GetEntries(), args.end_entry_index());
        const Long64_t nEntries = std::min(endEntry - args.begin_entry_index(), args.max_events());
        tools::ProgressReporter progressReporter(10, std::cout);
        progressReporter.SetTotalNumberOfEvents(static_cast<unsigned long>(nEntries));

        for (Long64_t n = 0; n < nEntries; ++n) {
            const Long64_t entry_index = args.begin_entry_index() + n;
            anaTuple.GetEntry(entry_index);

            for(size_t point_index = 0; point_index < points.size(); ++point_index) {
                outputs.at(point_index) = ComputeScore(anaTuple.data(), points.at(point_index));
            }

            outputTree->Fill();
            if (n % 100 == 0)
                progressReporter.Report(static_cast<unsigned long>(n + 1), false);
        }
        progressReporter.Report(static_cast<unsigned long>(nEntries), true);

        std::cout << "Writing output... " << std::flush;
        outputFile->Write();
        std::cout << "done" << std::endl;
    }

private:

    static std::vector<std::string> GetRequestedVariables(const std::string& feat_file) {
        std::ifstream infile(feat_file);
        std::cout << "Reading features from file: " << feat_file << "\nFeatures:";
        std::string line;
        std::vector<std::string> requested;
        while (std::getline(infile, line)) {
            std::cout << line << " ";
            requested.push_back(line);
        }
        std::cout << "\n";
        infile.close();
        return requested;
    }

    float ComputeScore(const bbtautau::AnaEvent& event, const bbtautau::HyperPoint& point)
    {
        using LorentzVectorPEP = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;
        using LorentzVector    = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;

        static constexpr double E_MASS  = 0.0005109989; //GeV
        static constexpr double MU_MASS = 0.1056583715; //GeV

        float tau1_mass;
        if(args.channel() == Channel::ETau)
            tau1_mass = E_MASS;
        else if(args.channel() == Channel::MuTau)
            tau1_mass = MU_MASS;
        else
            tau1_mass = event.tau1_m;

        const auto createLVec = [](float pt, float eta, float phi, float m) {
            const LorentzVectorPEP pep(pt, eta, phi, m);
            LorentzVector v(pep);
            return v;
        };

        const LorentzVector l_1 = createLVec(event.tau1_pt, event.tau1_eta, event.tau1_phi, tau1_mass),
                            l_2 = createLVec(event.tau2_pt, event.tau2_eta, event.tau2_phi, event.tau2_m),
                            met = createLVec(event.MET_pt, 0, event.MET_phi, 0),
                            b_1 = createLVec(event.b1_pt, event.b1_eta, event.b1_phi, event.b1_m),
                            b_2 = createLVec(event.b2_pt, event.b2_eta, event.b2_phi, event.b2_m),
                            vbf_1 = createLVec(event.VBF1_pt, event.VBF1_eta, event.VBF1_phi, event.VBF1_m),
                            vbf_2 = createLVec(event.VBF2_pt, event.VBF2_eta, event.VBF2_phi, event.VBF2_m),
                            svfit = createLVec(event.SVfit_pt, event.SVfit_eta, event.SVfit_phi, event.SVfit_m);

        const bool is_boosted = event.is_boosted;

        static const std::map<::analysis::Channel, ::Channel> ch_map = {
            { Channel::ETau, eTau }, { Channel::MuTau, muTau }, { Channel::TauTau, tauTau },
        };

        static const std::map<Period, ::Year> year_map = {
            { Period::Run2016, y16 }, { Period::Run2017, y17 }, { Period::Run2018, y18 },
        };

        const float res_mass = point.mass ? *point.mass : 125;

        Spin spin;
        if(point.spin) {
            if(*point.spin == 0)
                spin = radion;
            else if(*point.spin == 2)
                spin = graviton;
            else
                throw exception("Spin not supported.");
        } else {
            spin = nonres;
        }

        const float klambda = point.kl ? *point.kl : 1;

        const int n_vbf = event.has_VBF_pair ? 2 : 0;

        const float cv = 1, c2v = 1, c3 = 0;
        const bool cut_passed = true;

        const std::vector<float> feat_vals = evt_proc->process_as_vec(
                b_1, b_2, l_1, l_2, met, svfit, vbf_1, vbf_2, event.kinFit_m, event.kinFit_chi2, event.MT2,
                is_boosted, event.b1_DeepFlavour, event.b2_DeepFlavour, ch_map.at(args.channel()),
                year_map.at(args.period()), res_mass, spin, klambda, n_vbf,
                event.SVfit_valid > 0, event.kinFit_convergence > 0,
                event.b1_HHbtag, event.b2_HHbtag, event.VBF1_HHbtag, event.VBF2_HHbtag, event.b1_DeepFlavour_CvsL,
                event.b2_DeepFlavour_CvsL, event.VBF1_DeepFlavour_CvsL, event.VBF2_DeepFlavour_CvsL,
                event.b1_DeepFlavour_CvsB, event.b2_DeepFlavour_CvsB, event.VBF1_DeepFlavour_CvsB,
                event.VBF2_DeepFlavour_CvsB, cv, c2v, c3, cut_passed);
        return wrapper->predict(feat_vals, event.evt);
    }

private:
    Arguments args;
    std::shared_ptr<TFile> inputFile, outputFile;
    bbtautau::AnaTuple anaTuple;
    std::vector<bbtautau::HyperPoint> points;
    std::unique_ptr<EvtProc> evt_proc;
    std::unique_ptr<InfWrapper> wrapper;
};

} // namespace analysis

PROGRAM_MAIN(analysis::CalcDNN, analysis::Arguments)
