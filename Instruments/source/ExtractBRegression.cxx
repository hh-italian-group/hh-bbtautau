/*! Extract b-regression values to be aligned with AnaTuple.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "AnalysisTools/Core/include/ProgressReporter.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/KinFitInterface.h"
#include "hh-bbtautau/Analysis/include/AnaTuple.h"

namespace analysis {

struct Arguments {
    REQ_ARG(std::string, input);
    REQ_ARG(std::string, input_reg);
    REQ_ARG(Channel, channel);
    REQ_ARG(std::string, output);
    REQ_ARG(std::string, dataIdName);
    OPT_ARG(unsigned, n_threads, 1);
    OPT_ARG(Long64_t, max_events, std::numeric_limits<Long64_t>::max());
};

struct JetEntry {
    float pt, eta, phi, m, resolution, bRegCorr, bRegRes, hadronFlavour, DeepFlavour;
};

class ExtractBRegression {
public:
    using DataId = EventAnalyzerDataId;
    using Hash = size_t;
    using DataIdBiMap = boost::bimap<DataId, Hash>;

    ExtractBRegression(const Arguments& _args) :
        args(_args), inputFile(root_ext::OpenRootFile(args.input())),
        inputFileReg(root_ext::OpenRootFile(args.input_reg())),
        outputFile(root_ext::CreateRootFile(args.output())), anaTuple(ToString(args.channel()), inputFile.get(), true)
    {
        ROOT::EnableThreadSafety();
        if(args.n_threads() > 1)
            ROOT::EnableImplicitMT(args.n_threads());

        bbtautau::AnaAuxTuple aux_tuple(inputFile.get(), true);
        aux_tuple.GetEntry(0);
        refDataId = FindRefDataId(aux_tuple(), args.dataIdName());
    }

    void Run()
    {
        const std::string channelName = ToString(args.channel());

        const Long64_t nEntries = std::min(anaTuple.GetEntries(), args.max_events());

        std::set<ULong64_t> evt_set;
        std::set<Long64_t> entry_set;
        for (Long64_t n = 0; n < nEntries; ++n) {
            anaTuple.GetEntry(n);
            const auto& event = anaTuple.data();
            evt_set.insert(event.evt);
            entry_set.insert(n);
        }


        LoadRegressionCache(evt_set);

        outputFile->cd();
        auto outputTree = std::make_unique<TTree>(channelName.c_str(), channelName.c_str());

        #define DEF_BR(name) \
            float name; \
            outputTree->Branch(#name, &name, #name "/F")

        #define BR(name, idx) DEF_BR(b##idx##_##name)

        #define B_BR(idx) \
            BR(pt, idx); \
            BR(eta, idx); \
            BR(phi, idx); \
            BR(m, idx); \
            BR(resolution, idx); \
            BR(hadronFlavour, idx); \
            BR(DeepFlavour, idx); \
            BR(bRegCorr, idx); \
            BR(bRegRes, idx)

        DEF_BR(weight);
        DEF_BR(m_sv);
        DEF_BR(tau1_pt);
        DEF_BR(tau1_eta);
        DEF_BR(tau1_phi);
        DEF_BR(tau1_m);
        DEF_BR(tau2_pt);
        DEF_BR(tau2_eta);
        DEF_BR(tau2_phi);
        DEF_BR(tau2_m);
        DEF_BR(MET_pt);
        DEF_BR(MET_phi);
        DEF_BR(MET_cov_00);
        DEF_BR(MET_cov_01);
        DEF_BR(MET_cov_11);
        DEF_BR(kinFit_m);
        DEF_BR(kinFit_chi2);
        DEF_BR(kinFit_m_reg);
        DEF_BR(kinFit_chi2_reg);
        B_BR(1);
        B_BR(2);

        #undef BR
        #undef DEF_BR

        // const Long64_t nEntries = std::min(anaTuple.GetEntries(), args.max_events());
        tools::ProgressReporter progressReporter(10, std::cout);
        progressReporter.SetTotalNumberOfEvents(entry_set.size());
        // for (Long64_t n = 0; n < nEntries; ++n) {
        size_t n_processed = 0;
        for (Long64_t n : entry_set) {
            anaTuple.GetEntry(n);
            const auto& event = anaTuple.data();

            bool id_found = false;
            for(size_t dataId : event.dataIds) {
                if(dataId == refDataId) {
                    id_found = true;
                    break;
                }
            }
            if(id_found) {

                #define BR(name, idx) b##idx##_##name = b##idx.name
                #define B_FILL(idx) \
                    const auto b##idx = FindJetEntry(event.evt, event.b##idx##_pt, event.b##idx##_eta, \
                                                     event.b##idx##_phi, event.b##idx##_m, event.b##idx##_resolution, \
                                                     event.b##idx##_hadronFlavour, event.b##idx##_DeepFlavour); \
                    B_BR(idx)

                B_FILL(1);
                B_FILL(2);

                #undef B_BR
                #undef B_FILL

                weight = event.all_weights.at(0);
                m_sv = event.SVfit_m;
                tau1_pt = event.tau1_pt;
                tau1_eta = event.tau1_eta;
                tau1_phi = event.tau1_phi;
                tau1_m = event.tau1_m;
                tau2_pt = event.tau2_pt;
                tau2_eta = event.tau2_eta;
                tau2_phi = event.tau2_phi;
                tau2_m = event.tau2_m;
                MET_pt = event.MET_pt;
                MET_phi = event.MET_phi;
                MET_cov_00 = event.MET_cov_00;
                MET_cov_01 = event.MET_cov_01;
                MET_cov_11 = event.MET_cov_11;
                kinFit_m = event.kinFit_m;
                kinFit_chi2 = event.kinFit_chi2;

                auto kinFit_results_reg = GetKinfit(event, b1, b2);


                kinFit_m_reg = kinFit_results_reg.HasValidMass() ? static_cast<float>(kinFit_results_reg.mass) : -1;
                kinFit_chi2_reg = kinFit_results_reg.HasValidMass() ? static_cast<float>(kinFit_results_reg.chi2) : -1;
                outputTree->Fill();
            }
            if (n_processed % 100 == 0)
                progressReporter.Report(++n_processed, false);
        }
        progressReporter.Report(entry_set.size(), true);

        std::cout << "Writing output... " << std::flush;
        outputFile->Write();
        std::cout << "done" << std::endl;
    }

    kin_fit::FitResults GetKinfit(const bbtautau::AnaEvent& event, const JetEntry& b1_corr,
                                  const JetEntry& b2_corr) const
    {
        LorentzVectorM tau1(event.tau1_pt, event.tau1_eta, event.tau1_phi, event.tau1_m);
        LorentzVectorM tau2(event.tau2_pt, event.tau2_eta, event.tau2_phi, event.tau2_m);
        LorentzVectorM b1(event.b1_pt, event.b1_eta, event.b1_phi, event.b1_m);
        LorentzVectorM b2(event.b2_pt, event.b2_eta, event.b2_phi, event.b2_m);
        LorentzVectorM met(event.MET_pt, 0, event.MET_phi, 0);
        double met_px = met.Px(), met_py = met.Py();
        TMatrixD met_cov(2, 2);
        met_cov(0, 0) = event.MET_cov_00;
        met_cov(0, 1) = event.MET_cov_01;
        met_cov(1, 0) = event.MET_cov_01;
        met_cov(1, 1) = event.MET_cov_11;


        LorentzVectorM b1_reg = b1 * b1_corr.bRegCorr;
        LorentzVectorM b2_reg = b2 * b2_corr.bRegCorr;

        met_px -= (b1_reg.Px() - b1.Px()) + (b2_reg.Px() - b2.Px());
        met_py -= (b1_reg.Py() - b1.Py()) + (b2_reg.Py() - b2.Py());


        return kin_fit::FitProducer::FitImpl(ConvertVector(tau1), ConvertVector(tau2), ConvertVector(b1_reg),
                                             ConvertVector(b2_reg), TVector2(met_px, met_py), met_cov,
                                             b1_corr.bRegRes * b1_reg.pt(), b2_corr.bRegRes * b2_reg.pt(), 0);
    }

    void LoadRegressionCache(const std::set<ULong64_t>& evt_set)
    {
        #define VALUE(type, name) TTreeReaderValue<type> name(reader, #name)
        #define ARRAY(type, name) TTreeReaderArray<type> name(reader, #name)

        TTreeReader reader("Events", inputFileReg.get());
        VALUE(ULong64_t, evt);
        // ARRAY(Float_t, jet_pt);
        ARRAY(Float_t, jet_eta);
        ARRAY(Float_t, jet_phi);
        ARRAY(Float_t, jet_bRegCorr);
        ARRAY(Float_t, jet_bRegRes);

        #undef VALUE
        #undef ARRAY

        while(reader.Next()) {
            if(regCache.count(*evt))
                throw exception("Duplicated event id = %1%.") % (*evt);
            if(!evt_set.count(*evt)) continue;
            for(size_t n = 0; n < jet_eta.GetSize(); ++n) {
                // if(!(Jet_pt[n] > 15 && std::abs(Jet_eta[n]) < 2.4)) continue;
                JetEntry jet;
                jet.eta = jet_eta[n];
                jet.phi = jet_phi[n];
                jet.bRegCorr = jet_bRegCorr[n];
                jet.bRegRes = jet_bRegRes[n];
                regCache[*evt].push_back(jet);
            }
        }
    }

    JetEntry FindJetEntry(ULong64_t evt, float pt, float eta, float phi, float m, float resolution, int hadronFlavour,
                          float DeepFlavour) const
    {
        static constexpr double dR_thr = 0.01;
        static constexpr double dR2_thr = dR_thr * dR_thr;
        const auto iter = regCache.find(evt);
        if(iter != regCache.end()) {
            const auto& jets = iter->second;
            std::set<size_t> matched;
            for(size_t n = 0; n < jets.size(); ++n) {
                const double dEta = jets.at(n).eta - eta;
                const double dPhi = ROOT::Math::VectorUtil::Phi_mpi_pi(jets.at(n).phi - phi);
                const double dR2 = std::pow(dEta, 2) + std::pow(dPhi, 2);
                if(dR2 < dR2_thr)
                    matched.insert(n);
            }
            if(matched.size() == 1) {
                JetEntry jet = jets.at(*matched.begin());
                jet.pt = pt;
                jet.m = m;
                jet.resolution = resolution;
                jet.hadronFlavour = hadronFlavour;
                jet.DeepFlavour = DeepFlavour;
                return jet;
            }
            if(matched.size() > 1)
                throw exception("Multiple matching. evt=%1%, eta=%2%, phi=%3%.") % evt % eta % phi;
        }
        throw exception("Jet entry not found. evt=%1%, eta=%2%, phi=%3%.") % evt % eta % phi;
    }

    static size_t FindRefDataId(const bbtautau::AnaAux& aux, const std::string& name)
    {
        const size_t N = aux.dataIds.size();
        if(aux.dataId_names.size() != N)
            throw exception("Inconsistent dataId info in AnaAux tuple.");
        for(size_t n = 0; n < N; ++n) {
            const auto hash = aux.dataIds.at(n);
            const auto& dataId_name = aux.dataId_names.at(n);
            if(dataId_name == name)
                return hash;
        }
        throw exception("Hash not found.");
    }

private:
    Arguments args;
    std::shared_ptr<TFile> inputFile, inputFileReg, outputFile;
    bbtautau::AnaTuple anaTuple;
    std::map<ULong64_t, std::vector<JetEntry>> regCache;
    size_t refDataId;
};

} // namespace analysis

PROGRAM_MAIN(analysis::ExtractBRegression, analysis::Arguments)
