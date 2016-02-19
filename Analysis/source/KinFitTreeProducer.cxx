/*! Study of kinematic fit performance.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "Analysis/include/LightBaseFlatTreeAnalyzer.h"
#include "h-tautau/Analysis/include/KinFitTree.h"

class KinFitTreeProducer : public analysis::LightBaseFlatTreeAnalyzer {
public:
    KinFitTreeProducer(const std::string& inputFileName, const std::string& outputFileName)
         : LightBaseFlatTreeAnalyzer(inputFileName, outputFileName)
    {
        recalc_kinfit = true;
        kinFitTree = std::shared_ptr<ntuple::KinFitTree>(new ntuple::KinFitTree("kinFitTree"));
    }

    virtual ~KinFitTreeProducer() { kinFitTree->Write(); }

protected:
    virtual void AnalyzeEvent(const analysis::FlatEventInfo& eventInfo, analysis::EventCategory category) override
    {
        using analysis::EventCategory;
        if (!analysis::TwoJetsEventCategories_MediumBjets.count(category)) return;
        const analysis::EventRegion eventRegion = DetermineEventRegion(eventInfo, category);
        if (eventRegion == analysis::EventRegion::Unknown) return;
        if (eventInfo.fitResults.convergence != 0) return;
        //std::cout << "event: " << eventInfo.event->evt << ", convergence= " << eventInfo.fitResults.fit_bb_tt.convergence << std::endl;
        kinFitTree->evt() = eventInfo.event->evt;
        kinFitTree->channel() = static_cast<int>(eventInfo.channel);
        kinFitTree->eventCategory() = static_cast<int>(category);
        kinFitTree->eventRegion() = static_cast<int>(eventRegion);
        //bjet1 bjet2
        TLorentzVector bjet1 = eventInfo.bjet_momentums.at(0);
        TLorentzVector bjet2 = eventInfo.bjet_momentums.at(1);
        kinFitTree->pt_b1() = bjet1.Pt();
        kinFitTree->eta_b1() = bjet1.Eta();
        kinFitTree->phi_b1() = bjet1.Phi();
        kinFitTree->energy_b1() = bjet1.E();
        kinFitTree->csv_b1() = eventInfo.event->csv_Bjets.at(0);
        kinFitTree->pt_b2() = bjet2.Pt();
        kinFitTree->eta_b2() = bjet2.Eta();
        kinFitTree->phi_b2() = bjet2.Phi();
        kinFitTree->energy_b2() = bjet2.E();
        kinFitTree->csv_b2() = eventInfo.event->csv_Bjets.at(1);
        //tau1 tau2
        TLorentzVector tau1 = eventInfo.lepton_momentums.at(0);
        TLorentzVector tau2 = eventInfo.lepton_momentums.at(1);
        kinFitTree->pt_1() = tau1.Pt();
        kinFitTree->eta_1() = tau1.Eta();
        kinFitTree->phi_1() = tau1.Phi();
        kinFitTree->energy_1() = tau1.E();
        kinFitTree->pt_2() = tau2.Pt();
        kinFitTree->eta_2() = tau2.Eta();
        kinFitTree->phi_2() = tau2.Phi();
        kinFitTree->energy_2() = tau2.E();
        //MET
        kinFitTree->mvamet() = eventInfo.MET.Pt();
        kinFitTree->mvametphi() = eventInfo.MET.Phi();
        kinFitTree->mvacov00() = eventInfo.event->mvacov00;
        kinFitTree->mvacov01() = eventInfo.event->mvacov01;
        kinFitTree->mvacov10() = eventInfo.event->mvacov10;
        kinFitTree->mvacov11() = eventInfo.event->mvacov11;
        //Htautau
        kinFitTree->pt_Htt() = eventInfo.event->pt_sv_MC;
        kinFitTree->m_Htt() = eventInfo.event->m_sv_MC;
        kinFitTree->eta_Htt() = eventInfo.Htt.Eta();
        kinFitTree->phi_Htt() = eventInfo.Htt.Phi();

        kinFitTree->Fill();
    }

private:
    std::shared_ptr<ntuple::KinFitTree> kinFitTree;
};

