/*!
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <vector>
#include <map>

namespace MVA_Selections {

enum EventCategory { UnknownCategory, OneJet_ZeroBtag, OneJet_OneBtag, TwoJets_ZeroBtag, TwoJets_OneBtag, TwoJets_TwoBtag };
inline EventCategory EventCategoryFromString(const std::string& category_name)
{
    static std::map<std::string, EventCategory> m;
    if(!m.size()) {
        m["1jet0btag"] = OneJet_ZeroBtag;
        m["1jet1btag"] = OneJet_OneBtag;
        m["2jets0btag"] = TwoJets_ZeroBtag;
        m["2jets1btag"] = TwoJets_OneBtag;
        m["2jets2btag"] = TwoJets_TwoBtag;
    }
    if(!m.count(category_name))
        throw std::runtime_error("Unknown category name");
    return m[category_name];
}

enum Channel { ETau, MuTau, TauTau };
inline Channel ChannelFromString(const std::string& channel_name)
{
    static std::map<std::string, Channel> m;
    if(!m.size()) {
        m["eTau"] = ETau;
        m["muTau"] = MuTau;
        m["tauTau"] = TauTau;
    }
    if(!m.count(channel_name))
        throw std::runtime_error("Unknown channel");
    return m[channel_name];
}

enum MvaMethod { BDT, BDTMitFisher, BDTD };
inline MvaMethod MvaMethodFromString(const std::string& method_name)
{
    static std::map<std::string, MvaMethod> m;
    if(!m.size()) {
        m["BDT"] = BDT;
        m["BDTMitFisher"] = BDTMitFisher;
        m["BDTD"] = BDTD;
    }
    if(!m.count(method_name))
        throw std::runtime_error("Unknown MVA method");
    return m[method_name];
}

struct ParamId {
    Channel channel;
    EventCategory category;
    MvaMethod method;

    ParamId(Channel _channel, EventCategory _category, MvaMethod _method)
        : channel(_channel), category(_category), method(_method) {}

    bool operator< (const ParamId& other) const
    {
        if(channel < other.channel) return true;
        if(channel > other.channel) return false;
        if(category < other.category) return true;
        if(category > other.category) return false;
        return method < other.method;
    }
};


typedef std::map<std::string, size_t> var_map;
typedef std::vector<std::string> str_vector;
typedef std::map<ParamId, str_vector> var_list;
typedef std::map<ParamId, var_map> var_map_list;

inline var_map MakeVarMap(const str_vector& vars)
{
    var_map result;
    for(size_t n = 0; n < vars.size(); ++n)
        result[vars[n]] = n;
    return result;
}

const str_vector& Input_Variables(const ParamId& key)
{
    static var_list l;
    if(!l.size()) {
        {
            str_vector& v = l[ParamId(MuTau, TwoJets_ZeroBtag, BDT)];
            v.push_back("pt_mu");
            v.push_back("pt_tau");
            //v.push_back("pt_b1");
            //v.push_back("pt_b2");
            //v.push_back("DR_bb");
            v.push_back("DPhi_BBMET");
            v.push_back("DR_ll");
            v.push_back("Pt_Htt");
            //v.push_back("DR_HBBHTT");
            //v.push_back("Pt_Hbb");
            //v.push_back("DeltaPhi_METTT");
            //v.push_back("PtH");
            v.push_back("mT2");
            v.push_back("mT1");
            v.push_back("Pt_Htt_MET");
        }
        {
            str_vector& v = l[ParamId(MuTau, TwoJets_ZeroBtag, BDTMitFisher)];
            v.push_back("pt_mu");
            v.push_back("pt_tau");
            //v.push_back("pt_b1");
            //v.push_back("pt_b2");
            //v.push_back("DR_bb");
            //v.push_back("DPhi_BBMET");
            v.push_back("DR_ll");
            v.push_back("Pt_Htt");
            //v.push_back("DR_HBBHTT");
            //v.push_back("Pt_Hbb");
            //v.push_back("DeltaPhi_METTT");
            //v.push_back("PtH");
            v.push_back("mT2");
            v.push_back("mT1");
            v.push_back("Pt_Htt_MET");
        }
        {
            str_vector& v = l[ParamId(MuTau, TwoJets_ZeroBtag, BDTD)];
            v.push_back("pt_mu");
            v.push_back("pt_tau");
            //v.push_back("pt_b1");
            //v.push_back("pt_b2");
            //v.push_back("DR_bb");
            v.push_back("DPhi_BBMET");
            v.push_back("DR_ll");
            v.push_back("Pt_Htt");
            //v.push_back("DR_HBBHTT");
            //v.push_back("Pt_Hbb");
            //v.push_back("DeltaPhi_METTT");
            //v.push_back("PtH");
            v.push_back("mT2");
            v.push_back("mT1");
            v.push_back("Pt_Htt_MET");
        }
        {
            str_vector& v = l[ParamId(MuTau, TwoJets_OneBtag, BDT)];
            v.push_back("pt_mu");
            v.push_back("pt_tau");
            //v.push_back("pt_b1");
            //v.push_back("pt_b2");
            v.push_back("DR_bb");
            //v.push_back("DPhi_BBMET");
            v.push_back("DR_ll");
            v.push_back("Pt_Htt");
            //v.push_back("DR_HBBHTT");
            //v.push_back("Pt_Hbb");
            //v.push_back("DeltaPhi_METTT");
            //v.push_back("PtH");
            v.push_back("mT2");
            v.push_back("mT1");
            v.push_back("Pt_Htt_MET");
        }
        {
            str_vector& v = l[ParamId(MuTau, TwoJets_OneBtag, BDTMitFisher)];
            v.push_back("pt_mu");
            v.push_back("pt_tau");
            //v.push_back("pt_b1");
            //v.push_back("pt_b2");
            v.push_back("DR_bb");
            //v.push_back("DPhi_BBMET");
            v.push_back("DR_ll");
            //v.push_back("Pt_Htt");
            //v.push_back("DR_HBBHTT");
            //v.push_back("Pt_Hbb");
            v.push_back("DeltaPhi_METTT");
            //v.push_back("PtH");
            v.push_back("mT2");
            v.push_back("mT1");
            v.push_back("Pt_Htt_MET");
        }
        {
            str_vector& v = l[ParamId(MuTau, TwoJets_OneBtag, BDTD)];
            v.push_back("pt_mu");
            v.push_back("pt_tau");
            //v.push_back("pt_b1");
            v.push_back("pt_b2");
            //v.push_back("DR_bb");
            //v.push_back("DPhi_BBMET");
            v.push_back("DR_ll");
            v.push_back("Pt_Htt");
            //v.push_back("DR_HBBHTT");
            //v.push_back("Pt_Hbb");
            //v.push_back("DeltaPhi_METTT");
            //v.push_back("PtH");
            v.push_back("mT2");
            v.push_back("mT1");
            v.push_back("Pt_Htt_MET");
        }
        {
            str_vector& v = l[ParamId(MuTau, TwoJets_TwoBtag, BDT)];
            v.push_back("pt_mu");
            v.push_back("pt_tau");
            v.push_back("pt_b1");
            v.push_back("pt_b2");
            v.push_back("DR_bb");
            //v.push_back("DPhi_BBMET");
            v.push_back("DR_ll");
            //v.push_back("Pt_Htt");
            //v.push_back("DR_HBBHTT");
            //v.push_back("Pt_Hbb");
            //v.push_back("DeltaPhi_METTT");
            //v.push_back("PtH");
            v.push_back("mT2");
            //v.push_back("mT1");
            v.push_back("Pt_Htt_MET");
        }
        {
            str_vector& v = l[ParamId(MuTau, TwoJets_TwoBtag, BDTMitFisher)];
            v.push_back("pt_mu");
            v.push_back("pt_tau");
            v.push_back("pt_b1");
            v.push_back("pt_b2");
            v.push_back("DR_bb");
            //v.push_back("DPhi_BBMET");
            v.push_back("DR_ll");
            //v.push_back("Pt_Htt");
            //v.push_back("DR_HBBHTT");
            //v.push_back("Pt_Hbb");
            //v.push_back("DeltaPhi_METTT");
            //v.push_back("PtH");
            v.push_back("mT2");
            //v.push_back("mT1");
            v.push_back("Pt_Htt_MET");
        }
        {
            str_vector& v = l[ParamId(MuTau, TwoJets_TwoBtag, BDTD)];
            v.push_back("pt_mu");
            v.push_back("pt_tau");
            v.push_back("pt_b1");
            v.push_back("pt_b2");
            v.push_back("DR_bb");
            //v.push_back("DPhi_BBMET");
            v.push_back("DR_ll");
            //v.push_back("Pt_Htt");
            //v.push_back("DR_HBBHTT");
            //v.push_back("Pt_Hbb");
            //v.push_back("DeltaPhi_METTT");
            //v.push_back("PtH");
            v.push_back("mT2");
            //v.push_back("mT1");
            //v.push_back("Pt_Htt_MET");
        }
    }
    if(!l.count(key))
        throw std::runtime_error("Unknown combination channel-category");
    return l[key];
}

const var_map& Input_Variables_Map(const ParamId& key)
{
    static var_map_list l;
    if(!l.count(key))
        l[key] = MakeVarMap(Input_Variables(key));
    return l[key];
}

} // namespace MVA_Selections
