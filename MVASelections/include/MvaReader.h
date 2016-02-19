/*!
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "h-tautau/Analysis/include/FlatTree.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/exception.h"
#include "h-tautau/Analysis/include/Particles.h"

#include <TLorentzVector.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "MVA_variables.h"

namespace MVA_Selections {

class MvaReader;
typedef std::shared_ptr<MvaReader> MvaReaderPtr;

class MvaReader {
public:
    MvaReader(const ParamId& paramId, const std::string& mvaXMLfile)
        : reader(new TMVA::Reader("Silent")), var_names_map(&Input_Variables_Map(paramId))
    {
        //        TMVA::Tools::Instance();
        std::ostringstream ss;
        ss << paramId.channel << "_" << paramId.category << "_" << paramId.method;
        label = ss.str();
        const str_vector& var_names = Input_Variables(paramId);
        vars.assign(var_names.size(), 0);
        for(size_t n = 0; n < var_names.size(); ++n)
            reader->AddVariable(var_names.at(n), &vars.at(n));

        reader->BookMVA(label.c_str(), mvaXMLfile.c_str());
    }

    double GetMva (const TLorentzVector& l1, const TLorentzVector& l2,
                   const TLorentzVector& b1, const TLorentzVector& b2, const TLorentzVector& MET)
    {
        const TLorentzVector BB = b1 + b2;
        const TLorentzVector TT = l1 + l2;
        const TLorentzVector H = BB + TT;
        const TLorentzVector TT_MET = TT + MET;

        vars.assign(vars.size(), 0);
        SetMvaInput("pt_mu", l1.Pt());
        SetMvaInput("pt_tau", l2.Pt());
        SetMvaInput("pt_b1", b1.Pt());
        SetMvaInput("pt_b2", b2.Pt());
        SetMvaInput("DR_bb", b1.DeltaR(b2));
        SetMvaInput("DPhi_BBMET", MET.DeltaPhi(BB));
        SetMvaInput("DR_ll", l1.DeltaR(l2));
        SetMvaInput("Pt_Htt", TT.Pt());
        SetMvaInput("DR_HBBHTT", TT.DeltaR(BB));
        SetMvaInput("Pt_Hbb", BB.Pt());
        SetMvaInput("DeltaPhi_METTT", MET.DeltaPhi(TT));
        SetMvaInput("PtH", H.Pt());
        SetMvaInput("mT2", analysis::Calculate_MT(l2, MET.Pt(), MET.Phi()));
        SetMvaInput("mT1", analysis::Calculate_MT(l1, MET.Pt(), MET.Phi()));
        SetMvaInput("Pt_Htt_MET", TT_MET.Pt());

        return reader->EvaluateMVA(label.c_str());
    }

private:
    void SetMvaInput(const std::string& name, float value)
    {
        if(var_names_map->count(name))
           vars.at(var_names_map->find(name)->second) = value;
    }

private:
    std::shared_ptr<TMVA::Reader> reader;
    std::vector<float> vars;
    const var_map* var_names_map;
    std::string label;

public:
    static MvaReaderPtr Get(const std::string& channel_name, const std::string& event_category_name, MvaMethod mva_method)
    {
        typedef std::map<ParamId, std::string> FileNameMap;
        typedef std::map<ParamId, MvaReaderPtr> MvaReaderMap;

        static const std::string path = "MVASelections/weights/";

        static FileNameMap file_names;
        if(!file_names.size()) {
            file_names[ParamId(MuTau, TwoJets_ZeroBtag, BDT)] = "TMVA_mutau_2jets0btag_BDT_8var_BDT.weights.xml";
            file_names[ParamId(MuTau, TwoJets_ZeroBtag, BDTD)] = "TMVA_mutau_2jets0btag_BDTD_8var_BDTD.weights.xml";
            file_names[ParamId(MuTau, TwoJets_ZeroBtag, BDTMitFisher)] = "TMVA_mutau_2jets0btag_BDTMitFisher_7var_BDTMitFisher.weights.xml";

            file_names[ParamId(MuTau, TwoJets_OneBtag, BDT)] = "TMVA_mutau_2jets1btag_BDT_8var_BDT.weights.xml";
            file_names[ParamId(MuTau, TwoJets_OneBtag, BDTD)] = "TMVA_mutau_2jets1btag_BDTD_8var_BDTD.weights.xml";
            file_names[ParamId(MuTau, TwoJets_OneBtag, BDTMitFisher)] = "TMVA_mutau_2jets1btag_BDTMitFisher_8var_BDTMitFisher.weights.xml";

            file_names[ParamId(MuTau, TwoJets_TwoBtag, BDT)] = "TMVA_mutau_2jets2btag_BDT_8var_BDT.weights.xml";
            file_names[ParamId(MuTau, TwoJets_TwoBtag, BDTD)] = "TMVA_mutau_2jets2btag_BDTD_7var_BDTD.weights.xml";
            file_names[ParamId(MuTau, TwoJets_TwoBtag, BDTMitFisher)] = "TMVA_mutau_2jets2btag_BDTMitFisher_8var_BDTMitFisher.weights.xml";
        }

        if(event_category_name == "Inclusive")
            return MvaReaderPtr();

        const ParamId key(ChannelFromString(channel_name), EventCategoryFromString(event_category_name), mva_method);
        if(!file_names.count(key))
            return MvaReaderPtr();

        static MvaReaderMap readers;
        if(!readers.count(key)) {
            const std::string full_file_name = path + file_names.at(key);
            readers[key] = MvaReaderPtr(new MvaReader(key, full_file_name));
        }

        return readers.at(key);
    }
};

} // namespace MVA_Selections
