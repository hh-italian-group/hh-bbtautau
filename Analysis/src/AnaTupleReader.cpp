/*! Read AnaTuples.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/AnaTupleReader.h"
#include "AnalysisTools/Core/include/RootExt.h"

namespace analysis {
namespace bbtautau {

namespace {

template <typename Result, typename ...Args>
std::function<Result(Args...)> make_function(Result (*fn)(Args...))
{
    return [=](auto&& ...args) { return (*fn)(std::forward<Args>(args)...); };
}

template <class Class, typename Result, typename ...Args>
std::function<Result(Args...)> bind_this(const Class& obj, Result (Class::*fn)(Args...) const)
{
    return [=](auto&& ...args) { return (obj.*fn)(std::forward<Args>(args)...); };
}

using UncMap = EventTagCreator::UncMap;
using RDF = AnaTupleReader::RDF;
void FillUncMap(UncMap& /*unc_map*/, const std::vector<UncertaintySource>& /*unc_sources*/, size_t /*index*/)
{
}

template<typename ...Args>
void FillUncMap(UncMap& unc_map, const std::vector<UncertaintySource>& unc_sources, size_t index,
                float up, float down, Args... args)
{
    static constexpr float def_val = std::numeric_limits<float>::lowest();
    if(index >= unc_sources.size())
        return;
    if(up != def_val)
        unc_map.insert({{unc_sources.at(index), UncertaintyScale::Up}, up});
    if(down == def_val && up != def_val)
        down = up;
    if(down != def_val)
        unc_map.insert({{unc_sources.at(index), UncertaintyScale::Down}, down});
    if(up == def_val && down != def_val)
        throw exception("Invalid uncertainty definition");

    FillUncMap(unc_map, unc_sources, index + 1, args...);
}

template<typename ...Args>
UncMap CreateUncMap(const std::vector<UncertaintySource>& unc_sources, Args... args)
{
    UncMap unc_map;
    if(!unc_sources.empty())
        FillUncMap(unc_map, unc_sources, 0, args...);
    return unc_map;
}

template<std::size_t N, typename = std::make_index_sequence<2*N>>
struct UncMapCreator;


template <typename Result, typename ...Args>
struct StdFunction {
    using Fn = typename std::function<Result(Args...)>;
};

template<std::size_t N, std::size_t... S>
struct UncMapCreator<N, std::index_sequence<S...>> {
    template<typename T, size_t>
    using Type = T;

    using Fn = typename StdFunction<UncMap, Type<float, S>...>::Fn;

    static RDF CallDefine(RDF& df, const std::vector<UncertaintySource>& unc_sources,
                          const std::vector<std::string>& column_names)
    {
        const Fn fn = [=](Type<float, S>... args) { return CreateUncMap(unc_sources, args...); };
        std::vector<std::string> column_names_cp = column_names;
        if(!column_names.size()) {
            return df.Define("unc_map", [](){ return UncMap(); }, {});
        }
        const std::string& default_name = column_names.back();
        column_names_cp.resize(2 * N, default_name);
        return df.Define("unc_map", fn, column_names_cp);
    }
};
}


std::vector<std::shared_ptr<TFile>> AnaTupleReader::OpenFiles(const std::string& file_name,
                                                              const std::vector<std::string>& input_friends)
{
    std::vector<std::shared_ptr<TFile>> files;
    files.emplace_back(root_ext::OpenRootFile(file_name));
    for(const std::string& friend_file_name : input_friends)
        files.emplace_back(root_ext::OpenRootFile(friend_file_name));
    return files;
}

std::vector<std::shared_ptr<TTree>> AnaTupleReader::ReadTrees(Channel channel,
                                                              const std::vector<std::shared_ptr<TFile>>& files)
{
    std::vector<std::shared_ptr<TTree>> trees;
    int i=0;
    for(const auto& file : files) {
        trees.emplace_back(root_ext::ReadObject<TTree>(*file, ToString(channel)));
        if(trees.size() > 1)
            {
                trees.front()->AddFriend(trees.back().get(),("friend_"+std::to_string(i)).c_str());
                i++;
            }
    }
    return trees;
}



AnaTupleReader::AnaTupleReader(const std::string& file_name, Channel channel, NameSet& active_var_names,
                               const std::vector<std::string>& input_friends, const EventTagCreator& event_tagger,
                               const std::string& _mdnn_version, std::set<UncertaintySource>& _norm_unc_sources, int& _hastune) :
        files(OpenFiles(file_name, input_friends)), trees(ReadTrees(channel, files)), dataFrame(*trees.front()),
        df(dataFrame), mdnn_version(_mdnn_version), norm_unc(_norm_unc_sources), hastune(_hastune)
{
    std::mutex mutex;
    auto extract_names = [&](const std::vector<unsigned>& hashes, const std::vector<std::string>& names,
                             DatasetBiMap& name_map) -> bool {
        std::lock_guard<std::mutex> lock(mutex);
        const size_t N = hashes.size();
        if(names.size() != N)
            throw exception("Inconsistent info in AnaAux tuple.");
        for(size_t n = 0; n < N; ++n) {
            const auto hash = hashes.at(n);
            const auto& name = names.at(n);
            if(name_map.right.count(hash) || name_map.left.count(name)) {
                if(name_map.right.count(hash) && name_map.right.at(hash) != name)
                    throw exception("Duplicated hash = %1% in AnaAux tuple for name = %2%.\n"
                                    "This hash is already assigned to %3%") % hash % name % name_map.right.at(hash);
                if(name_map.left.count(name) && name_map.left.at(name) != hash)
                    throw exception("Duplicated name = '%1%' in AnaAux tuple.") % name;
            } else {
                name_map.insert({name, hash});
            }
            //std::cout << names.at(n) << std::endl;
        }
        return true;
    };

    //using namespace std::placeholders;
    ROOT::RDataFrame aux_df("aux", files.front().get());
    aux_df.Foreach([&](const std::vector<unsigned>& hashes, const std::vector<std::string>& names) {
                   extract_names(hashes, names, known_datasets); }, { "dataset_hashes", "dataset_names" });
    DatasetBiMap known_regions_str;
    aux_df.Foreach([&](const std::vector<unsigned>& hashes, const std::vector<std::string>& names) {
                   extract_names(hashes, names, known_regions_str); }, { "region_hashes", "region_names" });

    for(const auto& [region_str, hash] : known_regions_str.left)
        known_regions.insert({Parse<EventRegion>(region_str), hash});

    for(const auto& column : df.GetColumnNames())
        branch_types[column] = df.GetColumnType(column);

    DefineBranches(active_var_names, active_var_names.empty(), event_tagger);
    if(active_var_names.empty()) {
        std::vector<std::vector<std::string>> names = {
            df.GetColumnNames(),
        };
        for(auto& other_df : skimmed_df)
            names.push_back(other_df.GetColumnNames());
        for(const auto& name_set : names) {
            active_var_names.insert(name_set.begin(), name_set.end());
        }
    } else {
        for(const auto& var_name : active_var_names) {
            if(var_name.back() == '+') {
                const std::string name = var_name.substr(0, var_name.size() - 1);
                for(const auto& column : df.GetColumnNames()) {
                    if(column.rfind(name, 0) == 0)
                        parametric_vars[name].insert(column);
                }
            }
        }
        for(const auto& [name, columns] : parametric_vars) {
            active_var_names.erase(name + "+");
            for(const auto& column : columns)
                active_var_names.insert(column);
        }
    }
}

void AnaTupleReader::DefineBranches(const NameSet& active_var_names, bool all, const EventTagCreator& event_tagger)
{
    const auto Define = [&](RDF& target_df, const std::string& var, auto expr,
                              const std::vector<std::string>& columns, bool force = false) {
        if(force || all || active_var_names.count(var))
            target_df = target_df.Define(var, expr, columns);
    };

    const auto Filter = [](RDF& target_df, const std::string& var) -> RDF {
        return target_df.Filter([](bool flag){ return flag; }, {var});
    };

    const auto FilterInt = [](RDF& target_df, const std::string& var) -> RDF {
        return target_df.Filter([](int flag) -> bool { return flag; }, {var});
    };
    const auto Find = [](RDF& target_df, std::string& name){
      std::vector<std::string> columns =  target_df.GetColumnNames();
      return std::find(columns.begin(), columns.end(), name) != columns.end() ;
    };
    const auto Sum = [](float a, float b) -> float { return a + b; };
    const auto Delta = [](float a, float b) -> float { return a - b; };

    const auto ReturnP4 = [](float pt, float eta, float phi, float mass) {
        return LorentzVectorM(pt, eta, phi, mass);
    };

    const auto ReturnMETP4 = [](float pt, float phi) {
        return LorentzVectorM(pt, 0, phi, 0);
    };

    const auto SumP4 = [](const LorentzVectorM& p4_1, const LorentzVectorM& p4_2) {
        return p4_1 + p4_2;
    };

    const auto GetPt = [](const LorentzVectorM& p4) -> float { return static_cast<float>(p4.pt()); };
    const auto GetEta = [](const LorentzVectorM& p4) -> float { return static_cast<float>(p4.eta()); };
    const auto GetPhi = [](const LorentzVectorM& p4) -> float { return static_cast<float>(p4.phi()); };
    const auto GetMass = [](const LorentzVectorM& p4) -> float { return static_cast<float>(p4.mass()); };

    const auto DeltaPhi = [](const LorentzVectorM& p4_1, const LorentzVectorM& p4_2) -> float {
        return static_cast<float>(ROOT::Math::VectorUtil::DeltaPhi(p4_1, p4_2));
    };
    const auto AbsDeltaPhi = [](const LorentzVectorM& p4_1, const LorentzVectorM& p4_2) -> float {
        return static_cast<float>(std::abs(ROOT::Math::VectorUtil::DeltaPhi(p4_1, p4_2)));
    };
    const auto DeltaEta = [](const LorentzVectorM& p4_1, const LorentzVectorM& p4_2) -> float {
        return static_cast<float>(p4_1.eta() - p4_2.eta());
    };
    const auto DeltaR = [](const LorentzVectorM& p4_1, const LorentzVectorM& p4_2) -> float {
        return static_cast<float>(ROOT::Math::VectorUtil::DeltaR(p4_1, p4_2));
    };

    const auto DefineP4 = [&](RDF& target_df, const std::string& prefix) {
        Define(target_df, prefix + "_p4", ReturnP4,
               { prefix + "_pt", prefix + "_eta", prefix + "_phi", prefix + "_m" }, true);
    };


    const auto _Calculate_MT = [](const LorentzVectorM& p4, const LorentzVectorM& MET_p4) -> float {
        return static_cast<float>(Calculate_MT(p4, MET_p4));
    };

    DefineP4(df, "tau1");
    DefineP4(df, "tau2");
    Define(df, "Htt_p4", SumP4, { "tau1_p4", "tau2_p4" }, true);
    Define(df, "MET_p4", ReturnMETP4, {"MET_pt", "MET_phi"}, true);
    Define(df, "HttMET_p4", SumP4, { "Htt_p4", "MET_p4" }, true);
    Define(df, "m_tt_vis", GetMass, {"Htt_p4"}, true);

    DefineP4(df, "b1");
    DefineP4(df, "b2");
    Define(df, "Hbb_p4", SumP4, { "b1_p4", "b2_p4" }, true);
    Define(df, "m_bb", GetMass, {"Hbb_p4"}, true);

    DefineP4(df, "SVfit");


    DefineP4(df, "VBF1");
    DefineP4(df, "VBF2");

    const auto convertToDataId = [&](unsigned dataset, unsigned event_region, int unc_source, int unc_scale) {
        return DataId(GetDatasetByHash(dataset), GetRegionByHash(event_region),
                      static_cast<UncertaintySource>(unc_source), static_cast<UncertaintyScale>(unc_scale));
    };
    Define(df, "dataId", convertToDataId, { "dataset", "event_region", "unc_source", "unc_scale" }, true);

    const auto Full_branch_name=[&](std::string branch_name){
        std::string full_name;
        std::vector<std::string> columns =  df.GetColumnNames();
        if(std::find(columns.begin(), columns.end(), branch_name) != columns.end() ) {
            full_name= branch_name;
        }
        else{
            for(unsigned int i=0; i<files.size(); i++){
              std::string name = "friend_"+std::to_string(i)+"."+branch_name;
                if(Find(df, name))
                        full_name="friend_"+std::to_string(i)+"."+branch_name;
            }
        }
        return full_name;
    };

    if(mdnn_version.empty()){
        auto fake_vbf_cat = [](){return std::make_pair(-1.f, analysis::VBF_Category::None);};
        Define(df, "vbf_cat", fake_vbf_cat, {}, true);
    }
    else{
        std::vector<std::string> mdnn_branches;
        static const std::vector<std::string> mdnn_score_names = {
            "_tt_dl", "_tt_sl", "_tt_lep", "_tt_fh", "_dy", "_hh_ggf", "_tth", "_tth_tautau", "_tth_bb", "_hh_vbf",
            "_hh_vbf_c2v", "_hh_vbf_sm"
        };
        for(const auto& score_name : mdnn_score_names) {
            //const std::string branch = "mdnn_" + mdnn_version + score_name;
            std::string branch = Full_branch_name("mdnn_" + mdnn_version + score_name);
            mdnn_branches.push_back(branch);
        }
        const auto vbf_cat = make_function(&EventTagCreator::FindVBFCategory);
        Define(df, "vbf_cat", vbf_cat, mdnn_branches, true);
    }
    Define(df, "mdnn_score", [](const std::pair<float, analysis::VBF_Category>& vbf_cat) { return vbf_cat.first; },
           {"vbf_cat"}, true);

    const auto GetTune = [&](unsigned dataset){
      bool is_TuneCP5=0;
      if(hastune==2){
        std::vector<std::string> datasets_tuneCP5 = {"TTTo2L2Nu", "TTToSemiLeptonic", "TTToHadronic", "ST_tW_antitop", "ST_tW_top", "ST_t-channel_antitop", "ST_t-channel_top"};
        std::vector<unsigned> datasets_tuneCP5_ids;
        for(auto& dataset : datasets_tuneCP5){
            datasets_tuneCP5_ids.push_back(known_datasets.left.at(dataset));
        }
        int count_dataset=0;
        for(auto& k: datasets_tuneCP5_ids){
            if(k==dataset)
            count_dataset+=1;
        }
        is_TuneCP5 = (count_dataset>0) ? 1 : 0 ;
      }
      else if(hastune==0){
        is_TuneCP5 = 0;
      }
      return is_TuneCP5;
    };
    std::vector<std::string> columns =  df.GetColumnNames();
    if(std::find(columns.begin(), columns.end(), "is_TuneCP5")!=columns.end()) {
        hastune=1;
    }
    if (hastune==0 || hastune == 2 )
      Define(df, "is_TuneCP5", GetTune, {"dataset"}, true);
    /*
    if(hastune==2){ // 2016, but does not have tuneCP5
            std::vector<std::string> datasets_tuneCP5 = {"TTTo2L2Nu", "TTToSemiLeptonic", "TTToHadronic", "ST_tW_antitop", "ST_tW_top", "ST_t-channel_antitop", "ST_t-channel_top"};
            std::vector<unsigned> datasets_tuneCP5_ids;
            for(auto& dataset : datasets_tuneCP5){
                datasets_tuneCP5_ids.push_back(known_datasets.left.at(dataset));
            }
            const auto GetTune = [&](unsigned dataset){
                int count_dataset=0;
                for(auto& k: datasets_tuneCP5_ids){
                    if(k==dataset)
                        count_dataset+=1;
                }
                int is_TuneCP5 = (count_dataset>0) ? 1 : 0 ;
                return is_TuneCP5;
            };
    }
    else if(hastune==0) {
        auto fake_is_TuneCP5 = [](){return 0;};
        Define(df, "is_TuneCP5", fake_is_TuneCP5, {}, true);
    }*/

    const std::vector<UncertaintySource> norm_unc_sources(norm_unc.begin(),norm_unc.end());

    std::vector<std::string> norm_unc_names;
    for (auto i:norm_unc_sources){
        norm_unc_names.push_back("unc_"+analysis::ToString(i)+"_Up") ;
        norm_unc_names.push_back("unc_"+analysis::ToString(i)+"_Down") ;
    }
    df = UncMapCreator<40>::CallDefine(df, norm_unc_sources, norm_unc_names);

    const auto create_event_tags = bind_this(event_tagger, &EventTagCreator::CreateEventTags);
    Define(df, "event_tags", create_event_tags,
           { "dataId", "weight", "is_data", "weight_btag_Loose", "weight_btag_Medium", "weight_btag_Tight",
             "weight_btag_IterativeFit", "num_central_jets", "has_b_pair", "num_btag_Loose", "num_btag_Medium",
             "num_btag_Tight", "is_vbf", "is_boosted", "is_TuneCP5", "vbf_cat", "SVfit_p4", "unc_map", "m_bb", "m_tt_vis",
             "kinFit_convergence", "SVfit_valid" }, true);


    auto df_bb = Filter(df, "has_b_pair");


    auto df_vbf = Filter(df_bb, "has_VBF_pair");
    //DefineP4(df_vbf, "VBF1");
    //DefineP4(df_vbf, "VBF2");

    DefineP4(df, "central_jet1");
    DefineP4(df, "central_jet2");
    DefineP4(df, "central_jet3");
    DefineP4(df, "central_jet4");
    DefineP4(df, "central_jet5");

    auto df_sv = FilterInt(df, "SVfit_valid");


    auto df_bb_sv = FilterInt(df_bb, "SVfit_valid");
    //DefineP4(df_bb_sv, "SVfit");


    Define(df, "pt_H_tt", GetPt, {"Htt_p4"});
    Define(df, "eta_H_tt", GetEta, {"Htt_p4"});
    Define(df, "phi_H_tt", GetPhi, {"Htt_p4"});
    Define(df, "mt_1", _Calculate_MT, {"tau1_p4", "MET_p4"});
    Define(df, "mt_2", _Calculate_MT, {"tau2_p4", "MET_p4"});
    Define(df_sv, "MT_htautau", _Calculate_MT, {"SVfit_p4", "MET_p4"});

    Define(df, "dR_l1l2", DeltaR, {"tau1_p4", "tau2_p4"});
    Define(df, "abs_dphi_l1MET", AbsDeltaPhi, {"tau1_p4", "MET_p4"});
    Define(df_sv, "dphi_htautauMET", DeltaPhi, {"SVfit_p4", "MET_p4"});
    Define(df, "dR_l1l2MET", DeltaR, {"Htt_p4", "MET_p4"});
    Define(df_sv, "dR_l1l2Pt_htautau", [](const LorentzVectorM& p4_1, const LorentzVectorM& p4_2, float pt)
        { return ROOT::Math::VectorUtil::DeltaR(p4_1, p4_2) * pt; }, {"tau1_p4", "tau2_p4", "SVfit_pt"});
    Define(df, "mass_l1l2MET", GetMass, {"HttMET_p4"});
    Define(df, "pt_l1l2MET", GetPt, {"HttMET_p4"});
    Define(df, "p_zeta", [](const LorentzVectorM& tau1_p4, const LorentzVectorM& tau2_p4, const LorentzVectorM& MET_p4)
        { return Calculate_Pzeta(tau1_p4, tau2_p4, MET_p4); }, {"tau1_p4", "tau2_p4", "MET_p4"});
    Define(df, "p_zetavisible", [](const LorentzVectorM& tau1_p4, const LorentzVectorM& tau2_p4)
        { return Calculate_visiblePzeta(tau1_p4, tau2_p4); }, {"tau1_p4", "tau2_p4"});
    Define(df, "mt_tot", [](const LorentzVectorM& tau1_p4, const LorentzVectorM& tau2_p4, const LorentzVectorM& MET_p4)
        { return Calculate_TotalMT(tau1_p4, tau2_p4, MET_p4); }, {"tau1_p4", "tau2_p4", "MET_p4"});


    Define(df_bb, "pt_H_bb", GetPt, {"Hbb_p4"});
    Define(df_bb_sv, "dphi_hbbhtautau", DeltaPhi, {"Hbb_p4", "SVfit_p4"});
    Define(df_bb_sv, "deta_hbbhtautau", DeltaEta, {"Hbb_p4", "SVfit_p4"});
    Define(df_bb, "costheta_METhbb", [](const LorentzVectorM& MET_p4, const LorentzVectorM& Hbb_p4)
        { return four_bodies::Calculate_cosTheta_2bodies(MET_p4, Hbb_p4); }, {"MET_p4", "Hbb_p4"});
    Define(df_bb, "dR_b1b2", DeltaR, {"b1_p4", "b2_p4"});
    Define(df_bb, "p_zeta", [](const LorentzVectorM& b1_p4, const LorentzVectorM& b2_p4, const LorentzVectorM& Hbb_p4)
        { return four_bodies::Calculate_dR_boosted(b1_p4, b2_p4, Hbb_p4); }, {"b1_p4", "b2_p4", "Hbb_p4"});
    Define(df_bb, "dR_lj", [](const LorentzVectorM& tau1_p4, const LorentzVectorM& tau2_p4,
                              const LorentzVectorM& b1_p4, const LorentzVectorM& b2_p4)
        { return four_bodies::Calculate_min_dR_lj(tau1_p4, tau2_p4, b1_p4, b2_p4); },
        { "tau1_p4", "tau2_p4", "b1_p4", "b2_p4"});

    Define(df_bb, "mass_top_pair", [](const LorentzVectorM& tau1_p4, const LorentzVectorM& tau2_p4,
                                      const LorentzVectorM& b1_p4, const LorentzVectorM& b2_p4,
                                      const LorentzVectorM& MET_p4)
        { return four_bodies::Calculate_topPairMasses(tau1_p4, tau2_p4, b1_p4, b2_p4, MET_p4); },
        { "tau1_p4", "tau2_p4", "b1_p4", "b2_p4", "MET_p4"});

    Define(df_bb, "mass_top1", [](const std::pair<double, double>& m_top) { return m_top.first; },
        {"mass_top_pair"});
    Define(df_bb, "mass_top2", [](const std::pair<double, double>& m_top) { return m_top.second; },
        {"mass_top_pair"});
    Define(df_bb, "hh_btag_b1b2", Sum, {"b1_HHbtag", "b2_HHbtag"});
    Define(df_bb, "hh_btag_b1_minus_b2", Delta, {"b1_HHbtag", "b2_HHbtag"});
    Define(df_vbf, "hh_btag_VBF1VBF2", Sum, {"VBF1_HHbtag", "VBF2_HHbtag"});


    skimmed_df.push_back(df_bb);
    skimmed_df.push_back(df_vbf);
    skimmed_df.push_back(df_sv);
    skimmed_df.push_back(df_bb_sv);
}

const std::string& AnaTupleReader::GetDatasetByHash(unsigned hash) const
{
    const auto iter = known_datasets.right.find(hash);
    if(iter == known_datasets.right.end())
        throw exception("Dataset not found for hash = %1%") % hash;
    return iter->second;
}

const EventRegion& AnaTupleReader::GetRegionByHash(unsigned hash) const
{
    const auto iter = known_regions.right.find(hash);
    if(iter == known_regions.right.end())
        throw exception("Event region not found for hash = %1%") % hash;
    return iter->second;
}
size_t AnaTupleReader::GetNumberOfEntries() const { return static_cast<size_t>(trees.front()->GetEntries()); }
const AnaTupleReader::RDF& AnaTupleReader::GetDataFrame() const { return df; }
const std::list<AnaTupleReader::RDF>& AnaTupleReader::GetSkimmedDataFrames() const { return skimmed_df; }

const std::map<std::string, std::set<std::string>>& AnaTupleReader::GetParametricVariables() const
{
    return parametric_vars;
}

boost::optional<std::string> AnaTupleReader::TryGetBranchType(const std::string& branch_name) const
{
    boost::optional<std::string> branch_type;
    auto iter = branch_types.find(branch_name);
    if(iter != branch_types.end())
        branch_type = iter->second;
    return branch_type;
}

std::string HyperPoint::ToString()
{
    std::vector<std::string> points;
    if(spin)
        points.push_back("s" + ::analysis::ToString(*spin));
    if(mass)
        points.push_back("m" + ::analysis::ToString(*mass));
    if(kl)
        points.push_back("kl" + ::analysis::ToString(*kl));
    std::ostringstream ss;
    if(!points.size())
        throw exception("Empty hyperparameter.");
    ss << points.at(0);
    for(size_t n = 1; n < points.size(); ++n)
        ss << "_" << points.at(n);
    return ss.str();
}

} // namespace bbtautau
} // namespace analysis
