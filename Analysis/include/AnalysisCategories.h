/*! Definition of data and event categories used in HH->bbTauTau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <boost/optional/optional.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include "AnalysisTools/Core/include/Tools.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "h-tautau/Core/include/AnalysisTypes.h"

namespace analysis {

enum class SampleType { Data, MC, DY, QCD, TT, NonResHH };
ENUM_NAMES(SampleType) = {
    { SampleType::Data, "Data" }, { SampleType::MC, "MC" }, { SampleType::DY, "DY" }, { SampleType::QCD, "QCD" },
    { SampleType::TT, "TT" }, { SampleType::NonResHH, "NonResHH" }
};

enum class DYFitModel { None = 0, NbjetBins = 1, NbjetBins_htBins = 2 , NbjetBins_NjetBins = 3, NbjetBins_ptBins = 4};
ENUM_NAMES(DYFitModel) = {
    { DYFitModel::None, "None" } , { DYFitModel::NbjetBins, "NbjetBins" } ,
    { DYFitModel::NbjetBins_htBins, "NbjetBins_htBins"}, { DYFitModel::NbjetBins_NjetBins, "NbjetBins_NjetBins"},
    { DYFitModel::NbjetBins_ptBins, "NbjetBins_ptBins"}
};

struct EventRegion {
    static const EventRegion& Unknown();
    static const EventRegion& OS_Isolated();
    static const EventRegion& OS_AntiIsolated();
    static const EventRegion& SS_Isolated();
    static const EventRegion& SS_LooseIsolated();
    static const EventRegion& SS_AntiIsolated();
    static const EventRegion& SignalRegion();

    EventRegion() {}
    EventRegion& SetCharge(bool _os);
    EventRegion& SetLowerIso(DiscriminatorWP wp);
    EventRegion& SetUpperIso(DiscriminatorWP wp);
    bool HasCharge() const;
    bool HasLowerIso() const;
    bool HasUpperIso() const;

    DiscriminatorWP GetLowerIso() const;
    DiscriminatorWP GetUpperIso() const;
    bool GetCharge() const;

    bool Implies(const EventRegion& other) const;
    bool operator ==(const EventRegion& er) const;
    bool operator !=(const EventRegion& er) const;
    bool operator <(const EventRegion& er) const;

    std::string ToString() const;
    static EventRegion Parse(const std::string& str);

private:
    boost::optional<bool> os;
    boost::optional<DiscriminatorWP> iso_lower, iso_upper;
};

std::istream& operator>>(std::istream& s, EventRegion& eventRegion);
std::ostream& operator<<(std::ostream& os, const EventRegion& eventRegion);


struct EventCategory {
    static const EventCategory& Inclusive();

    EventCategory() {}
    explicit EventCategory(size_t _n_jets);
    EventCategory(size_t _n_jets, size_t _n_btag, bool _strict_n_btag, DiscriminatorWP _btag_wp);
    EventCategory(size_t _n_jets, size_t _n_btag, bool _strict_n_btag, DiscriminatorWP _btag_wp, bool _boosted);
    EventCategory(size_t _n_jets, size_t _n_btag, bool _strict_n_btag, DiscriminatorWP _btag_wp,
                  boost::optional<bool> _boosted, bool _is_VBF);

    bool HasJetConstraint() const;
    size_t N_jets() const;

    bool HasBtagConstraint() const;
    size_t N_btag() const;
    DiscriminatorWP BtagWP() const;
    bool IsN_btagStrict() const;

    bool HasBoostConstraint() const;
    bool IsBoosted() const;

    bool HasVBFConstraint() const;
    bool isVBF() const;
    bool operator ==(const EventCategory& ec) const;
    bool operator !=(const EventCategory& ec) const;
    bool operator <(const EventCategory& ec) const;

    std::string ToString() const;
    static EventCategory Parse(const std::string& str);

    bool Contains(size_t num_jets, const std::map<DiscriminatorWP, size_t>& num_btag, bool is_vbf,
                  bool is_boosted) const;
private:
    boost::optional<size_t> n_jets, n_btag ;
    boost::optional<bool> strict_n_btag;
    boost::optional<DiscriminatorWP> btag_wp;
    boost::optional<bool> boosted, is_VBF;
};

std::ostream& operator<<(std::ostream& os, const EventCategory& eventCategory);
std::istream& operator>>(std::istream& is, EventCategory& eventCategory);

#define DECL_MVA_SEL(z, n, first) MVA##n = n + first,
#define MVA_CUT_LIST(first, count) BOOST_PP_REPEAT(count, DECL_MVA_SEL, first)

enum class SelectionCut { mh = 0, mhVis = 1, mhMET = 2, KinematicFitConverged = 3, lowMET = 4, lowHT = 5, medHT = 6, highHT = 7,
                          vlowPt = 8, lowPt = 9, medPt = 10, highPt = 11, mtt = 12, MVA_CUT_LIST(13, 100) MVA_first = MVA0, MVA_last = MVA99 };

#undef MVA_CUT_LIST
#undef DECL_MVA_SEL

namespace detail {
inline std::map<SelectionCut, std::string> CreateSelectionCutNames()
{
    std::map<SelectionCut, std::string> names;
    names[SelectionCut::mh] = "mh";
    names[SelectionCut::mhVis] = "mhVis";
    names[SelectionCut::mhMET] = "mhMET";
    names[SelectionCut::KinematicFitConverged] = "KinematicFitConverged";
    names[SelectionCut::lowMET] = "lowMET";
    names[SelectionCut::lowHT] = "lowHT";
    names[SelectionCut::medHT] = "medHT";
    names[SelectionCut::highHT] = "highHT";
    names[SelectionCut::vlowPt] = "vlowPt";
    names[SelectionCut::lowPt] = "lowPt";
    names[SelectionCut::medPt] = "medPt";
    names[SelectionCut::highPt] = "highPt";
    names[SelectionCut::mtt] = "mtt";
    const size_t MVA_first_index = static_cast<size_t>(SelectionCut::MVA_first);
    const size_t MVA_last_index = static_cast<size_t>(SelectionCut::MVA_last);
    const size_t n_mva_cuts = MVA_last_index - MVA_first_index + 1;
    for(size_t n = 0; n < n_mva_cuts; ++n) {
        const SelectionCut cut = static_cast<SelectionCut>(n + MVA_first_index);
        std::ostringstream ss;
        ss << "MVA" << n;
        names[cut] = ss.str();
    }
    return names;
}
} // namespace detail

ENUM_NAMES(SelectionCut) = detail::CreateSelectionCutNames();

struct EventSubCategory {
    using BitsContainer = boost::multiprecision::uint128_t;
    static constexpr size_t MaxNumberOfCuts = std::numeric_limits<BitsContainer>::digits;

    static const EventSubCategory& NoCuts();

    EventSubCategory();

    bool HasCut(SelectionCut cut) const;
    bool Passed(SelectionCut cut) const;
    bool Failed(SelectionCut cut) const;

    EventSubCategory& SetCutResult(SelectionCut cut, bool result);

    BitsContainer GetPresenceBits() const;
    BitsContainer GetResultBits() const;

    bool operator ==(const EventSubCategory& sc) const;
    bool operator !=(const EventSubCategory& sc) const;
    bool operator <(const EventSubCategory& sc) const;

    bool Implies(const EventSubCategory& sc) const;
    bool TryGetLastMvaCut(SelectionCut& cut) const;

    std::string ToString(const std::map<SelectionCut, std::string>& sel_aliases = {}) const;
    static EventSubCategory Parse(const std::string& str);

private:
    static size_t GetIndex(SelectionCut cut);

private:
    BitsContainer presence, results;
    boost::optional<SelectionCut> last_mva_cut;
};

std::ostream& operator<<(std::ostream& os, const EventSubCategory& eventSubCategory);
std::istream& operator>>(std::istream& is, EventSubCategory& eventSubCategory);

using EventRegionSet = std::set<EventRegion>;
using EventCategorySet = std::set<EventCategory>;
using EventSubCategorySet = std::set<EventSubCategory>;
using Dataset = std::string;

} // namespace analysis
