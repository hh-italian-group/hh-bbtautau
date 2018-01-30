/*! Definition of data and event categories used in HH->bbTauTau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <boost/optional/optional.hpp>
#include <boost/bimap.hpp>
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/Tools.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include <boost/multiprecision/cpp_int.hpp>

namespace analysis {

enum class SampleType { Data, MC, DY, QCD, TT, NonResHH };
ENUM_NAMES(SampleType) = {
    { SampleType::Data, "Data" }, { SampleType::MC, "MC" }, { SampleType::DY, "DY" }, { SampleType::QCD, "QCD" },
    { SampleType::TT, "TT" }, { SampleType::NonResHH, "NonResHH" }
};

struct EventRegion {
    using EventRegionMapString = boost::bimap<std::string, EventRegion>;

    static const EventRegion& Unknown() { static const EventRegion er; return er; }
    static const EventRegion& OS_Isolated()
    {
        static const EventRegion er = EventRegion().SetCharge(true).SetLowerIso(DiscriminatorWP::Medium);
        return er;
    }

    static const EventRegion& OS_AntiIsolated()
    {
        static const EventRegion er =
                EventRegion().SetCharge(true).SetLowerIso(DiscriminatorWP::VLoose).SetUpperIso(DiscriminatorWP::Medium);
        return er;
    }

    static const EventRegion& SS_Isolated()
    {
        static const EventRegion er = EventRegion().SetCharge(false).SetLowerIso(DiscriminatorWP::Medium);
        return er;
    }

    static const EventRegion& SS_LooseIsolated()
    {
        static const EventRegion er = EventRegion().SetCharge(false).SetLowerIso(DiscriminatorWP::Loose);
        return er;
    }

    static const EventRegion& SS_AntiIsolated()
    {
        static const EventRegion er =
                EventRegion().SetCharge(false).SetLowerIso(DiscriminatorWP::VLoose).SetUpperIso(DiscriminatorWP::Medium);
        return er;
    }

    static const EventRegion& SignalRegion() { return OS_Isolated(); }

    EventRegion() {}

    EventRegion& SetCharge(bool _os) { os = _os; return *this; }
    EventRegion& SetLowerIso(DiscriminatorWP wp)
    {
        if(HasUpperIso() && iso_upper <= wp)
            throw exception("HasUpperIso - Iso Upper limit is not greater than Iso lower limit");
        iso_lower = wp;
        return *this;
    }
    EventRegion& SetUpperIso(DiscriminatorWP wp)
    {
        if(HasLowerIso() && wp <= iso_lower)
            throw exception("HasLowerIso - Iso Upper limit is not greater than Iso lower limit");
        iso_upper = wp;
        return *this;
    }

    bool HasCharge() const { return os.is_initialized(); }
    bool HasLowerIso() const { return iso_lower.is_initialized(); }
    bool HasUpperIso() const { return iso_upper.is_initialized(); }

    DiscriminatorWP GetLowerIso() const
    {
        if(!HasLowerIso())
            throw exception("Lower isolation bound not set.");
        return *iso_lower;
    }

    DiscriminatorWP GetUpperIso() const
    {
        if(!HasUpperIso())
            throw exception("Upper isolation bound not set.");
        return *iso_upper;
    }

    bool GetCharge() const
    {
        if(!HasCharge())
            throw exception("Charge info not set.");
        return *os;
    }


    bool Implies(const EventRegion& other) const
    {
        if(other.HasCharge() && (!HasCharge() || GetCharge() != other.GetCharge())) return false;
        if(other.HasLowerIso() && (!HasLowerIso() || GetLowerIso() < other.GetLowerIso())) return false;
        return !other.HasUpperIso() || (HasUpperIso() && GetUpperIso() <= other.GetUpperIso());
    }

    bool operator ==(const EventRegion& er) const { return os == er.os && iso_lower == er.iso_lower && iso_upper == er.iso_upper; }
    bool operator !=(const EventRegion& er) const { return !(*this == er); }
    bool operator <(const EventRegion& er) const
    {
        if(os != er.os) return os < er.os;
        if(iso_lower != er.iso_lower) return iso_lower < er.iso_lower;
        return iso_upper < er.iso_upper;
    }



    std::string ToString() const
    {
        if(*this == Unknown()) return "Unknown";
        if(!EventRegionMapToString().right.count(*this))
            throw exception("Unknown EventRegion. No conversion to String");
        return EventRegionMapToString().right.at(*this);
    }

    static EventRegion Parse(const std::string& str)
    {
        if(!EventRegionMapToString().left.count(str))
            throw exception("Unknown EventRegion = '%1%'.") % str;
        return EventRegionMapToString().left.at(str);
    }



private:
    boost::optional<bool> os;
    boost::optional<DiscriminatorWP> iso_lower, iso_upper;

    static const EventRegionMapString& EventRegionMapToString()
    {
        static EventRegionMapString predefined_regions;
        if(!predefined_regions.left.size()){
            predefined_regions.insert({ "Unknown", Unknown() });
            predefined_regions.insert({ "OS_Isolated", OS_Isolated() });
            predefined_regions.insert({ "OS_AntiIsolated", OS_AntiIsolated() });
            predefined_regions.insert({ "SS_LooseIsolated", SS_LooseIsolated() });
            predefined_regions.insert({ "SS_Isolated", SS_Isolated() });
            predefined_regions.insert({ "SS_AntiIsolated", SS_AntiIsolated() });
            predefined_regions.insert({ "SignalRegion", SignalRegion() });
        }
        return predefined_regions;
    }


};

inline std::istream& operator>>(std::istream& s, EventRegion& eventRegion)
{
    std::string str;
    s >> str;
    eventRegion = EventRegion::Parse(str);
    return s;
}

inline std::ostream& operator<<(std::ostream& os, const EventRegion& eventRegion)
{
    os << eventRegion.ToString();
    return os;
}

#define DEF_ES(name, n_jets, ...) \
    static const EventCategory& name() { static const EventCategory ec(n_jets, ##__VA_ARGS__); return ec; } \
    /**/

struct EventCategory {
    static const EventCategory& Inclusive() { static const EventCategory ec; return ec; }
    DEF_ES(TwoJets_Inclusive, 2)
    DEF_ES(TwoJets_ZeroBtag, 2, 0, DiscriminatorWP::Medium)
    DEF_ES(TwoJets_OneBtag, 2, 1, DiscriminatorWP::Medium)
    DEF_ES(TwoJets_TwoBtag, 2, 2, DiscriminatorWP::Medium)
    DEF_ES(TwoJets_ZeroLooseBtag, 2, 0, DiscriminatorWP::Loose)
    DEF_ES(TwoJets_OneLooseBtag, 2, 1, DiscriminatorWP::Loose)
    DEF_ES(TwoJets_TwoLooseBtag, 2, 2, DiscriminatorWP::Loose)
    DEF_ES(TwoJets_ZeroBtag_Resolved, 2, 0, DiscriminatorWP::Medium, false)
    DEF_ES(TwoJets_OneBtag_Resolved, 2, 1, DiscriminatorWP::Medium, false)
    DEF_ES(TwoJets_TwoBtag_Resolved, 2, 2, DiscriminatorWP::Medium, false)
    DEF_ES(TwoJets_TwoLooseBtag_Boosted, 2, 2, DiscriminatorWP::Loose, true)

    EventCategory() {}
    explicit EventCategory(size_t _n_jets) : n_jets(_n_jets) {}
    EventCategory(size_t _n_jets, size_t _n_btag, DiscriminatorWP _btag_wp) :
        n_jets(_n_jets), n_btag(_n_btag), btag_wp(_btag_wp)
    {
        if(n_btag > n_jets)
            throw exception("Number of btag can't be greater than number of jets");
    }

    EventCategory(size_t _n_jets, size_t _n_btag, DiscriminatorWP _btag_wp, bool _boosted) :
        n_jets(_n_jets), n_btag(_n_btag), btag_wp(_btag_wp), boosted(_boosted)
    {
        if(n_btag > n_jets)
            throw exception("Number of btag can't be greater than number of jets");
    }

    bool HasJetConstraint() const { return n_jets.is_initialized(); }
    size_t N_jets() const
    {
        if(!HasJetConstraint())
            throw exception("Jet constraint is not defined.");
        return *n_jets;
    }

    bool HasBtagConstraint() const { return n_btag.is_initialized(); }
    size_t N_btag() const
    {
        if(!HasBtagConstraint())
            throw exception("Btag constraint is not defined.");
        return *n_btag;
    }
    DiscriminatorWP BtagWP() const
    {
        if(!HasBtagConstraint())
            throw exception("Btag constraint is not defined.");
        return *btag_wp;
    }

    bool HasBoostConstraint() const { return boosted.is_initialized(); }
    bool IsBoosted() const
    {
        if(!HasBoostConstraint())
            throw exception("Boost constraint is not defined.");
        return *boosted;
    }

    bool operator ==(const EventCategory& ec) const
    {
        return n_jets == ec.n_jets && n_btag == ec.n_btag && btag_wp == ec.btag_wp && boosted == ec.boosted;
    }
    bool operator !=(const EventCategory& ec) const { return !(*this == ec); }
    bool operator <(const EventCategory& ec) const
    {
        if(n_jets != ec.n_jets) return n_jets < ec.n_jets;
        if(n_btag != ec.n_btag) return n_btag < ec.n_btag;
        if(btag_wp != ec.btag_wp) return btag_wp < ec.btag_wp;
        return boosted < ec.boosted;
    }

    std::string ToString() const
    {
        if(*this == Inclusive()) return "Inclusive";
        std::ostringstream s;
        s << *n_jets << "jets";
        if(HasBtagConstraint()) {
            s << *n_btag;
            if(*btag_wp != DiscriminatorWP::Medium)
                s << __DiscriminatorWP_short_names.EnumToString(*btag_wp);
            s << "btag";
        }
        if(HasBoostConstraint()) {
            const std::string boosted_str = IsBoosted() ? "B" : "R";
            s << boosted_str;
        }
        return s.str();
    }

    static EventCategory Parse(const std::string& str)
    {
        static const std::string numbers = "0123456789";
        static const std::string jets_suffix = "jets", btag_suffix = "btag";
        static const std::map<char, bool> boosted_suffix = { { 'R', false }, { 'B', true } };

        if(str == "Inclusive") return Inclusive();
        try {
            const size_t jet_pos = str.find_first_not_of(numbers);
            if(jet_pos == std::string::npos && str.substr(jet_pos, jets_suffix.size()) != "jets")
                throw exception("");
            const size_t n_jets = ::analysis::Parse<size_t>(str.substr(0, jet_pos));
            const size_t btag_str_pos = jet_pos + jets_suffix.size();
            if(str.size() == btag_str_pos)
                return EventCategory(n_jets);
            const size_t btag_wp_pos = str.find_first_not_of(numbers, btag_str_pos);
            if(btag_wp_pos == std::string::npos)
                throw exception("");
            const size_t n_btag = ::analysis::Parse<size_t>(str.substr(btag_str_pos, btag_wp_pos - btag_str_pos));
            const size_t btag_pos = str.find(btag_suffix, btag_wp_pos);
            if(btag_pos == std::string::npos)
                throw exception("");
            const DiscriminatorWP btag_wp = btag_wp_pos == btag_pos ? DiscriminatorWP::Medium
                    : __DiscriminatorWP_short_names.Parse(str.substr(btag_wp_pos, btag_pos - btag_wp_pos));
            const size_t boosted_pos = btag_pos + btag_suffix.size();
            if(str.size() == boosted_pos)
                return EventCategory(n_jets, n_btag, btag_wp);
            const char boosted_flag = str.at(boosted_pos);
            if(str.size() != boosted_pos + 1 || !boosted_suffix.count(boosted_flag))
                throw exception("");
            const bool is_boosted = boosted_suffix.at(boosted_flag);
            return EventCategory(n_jets, n_btag, btag_wp, is_boosted);
        }catch(exception& e) {
            throw exception("Invalid EventCategory '%1%'. %2%") % str % e.message();
        }
    }

private:
    boost::optional<size_t> n_jets, n_btag;
    boost::optional<DiscriminatorWP> btag_wp;
    boost::optional<bool> boosted;
};

#undef DEF_ES

std::ostream& operator<<(std::ostream& os, const EventCategory& eventCategory)
{
    os << eventCategory.ToString();
    return os;
}

std::istream& operator>>(std::istream& is, EventCategory& eventCategory)
{
    std::string str;
    is >> str;
    eventCategory = EventCategory::Parse(str);
    return is;
}

#define DECL_MVA_SEL(z, n, first) MVA##n = n + first,
#define MVA_CUT_LIST(first, count) BOOST_PP_REPEAT(count, DECL_MVA_SEL, first)

enum class SelectionCut { mh = 0, mhVis = 1, mhMET = 2, KinematicFitConverged = 3, lowMET = 4,
                          MVA_CUT_LIST(5, 100) MVA_first = MVA0, MVA_last = MVA99 };

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

    static const EventSubCategory& NoCuts() { static const EventSubCategory esc; return esc; }

    EventSubCategory() : presence(0), results(0) {}

    bool HasCut(SelectionCut cut) const
    {
        const BitsContainer mask = BitsContainer(1) << GetIndex(cut);
        return (presence & mask) != BitsContainer(0);
    }

    bool Passed(SelectionCut cut) const
    {
        if(!HasCut(cut))
            throw exception("Cut '%1%' is not defined.") % cut;
        const BitsContainer mask = BitsContainer(1) << GetIndex(cut);
        return (results & mask) != BitsContainer(0);
    }
    bool Failed(SelectionCut cut) const { return !Passed(cut); }

    EventSubCategory& SetCutResult(SelectionCut cut, bool result)
    {
        if(HasCut(cut))
            throw exception("Cut '%1%' is already defined.") % cut;
        const size_t index = GetIndex(cut);
        const BitsContainer mask = BitsContainer(1) << index;
        results = (results & ~mask) | (BitsContainer(result) << index);
        presence |= mask;
        if(cut >= SelectionCut::MVA_first && cut <= SelectionCut::MVA_last)
            last_mva_cut = cut;
        return *this;
    }

    BitsContainer GetPresenceBits() const { return presence; }
    BitsContainer GetResultBits() const { return results; }

    bool operator ==(const EventSubCategory& sc) const
    {
        return GetPresenceBits() == sc.GetPresenceBits() && GetResultBits() == sc.GetResultBits();
    }
    bool operator !=(const EventSubCategory& sc) const { return !(*this == sc); }
    bool operator <(const EventSubCategory& sc) const
    {
        if(GetPresenceBits() != sc.GetPresenceBits()) return GetPresenceBits() < sc.GetPresenceBits();
        return GetResultBits() < sc.GetResultBits();
    }

    bool Implies(const EventSubCategory& sc) const
    {
        const BitsContainer pres_a = GetPresenceBits(), pres_b = sc.GetPresenceBits();
        if(((pres_a ^ pres_b) & pres_b) != BitsContainer(0)) return false;
        const BitsContainer res_a = GetResultBits(), res_b = sc.GetResultBits();
        return (res_a & pres_b) == res_b;
    }

    bool TryGetLastMvaCut(SelectionCut& cut) const
    {
        if(!last_mva_cut.is_initialized()) return false;
        cut = *last_mva_cut;
        return true;
    }

    std::string ToString(const std::map<SelectionCut, std::string>& sel_aliases = {}) const
    {
        if(*this == NoCuts())
            return "NoCuts";
        std::ostringstream s;
        for(size_t n = 0; n < MaxNumberOfCuts; ++n) {
            const BitsContainer mask = BitsContainer(1) << n;
            if((presence & mask) == BitsContainer(0)) continue;
            if((results & mask) == BitsContainer(0)) s << "not";
            const SelectionCut cut = static_cast<SelectionCut>(n);
            if(sel_aliases.count(cut))
                s << sel_aliases.at(cut);
            else
                s << cut;
            s << "_";
        }
        std::string str = s.str();
        if(str.size())
            str.erase(str.size() - 1);
        return str;
    }

    static EventSubCategory Parse(const std::string& str)
    {
        if(str == "NoCuts")
            return NoCuts();
        EventSubCategory sub;
        const std::vector<std::string> cut_strings = SplitValueList(str, false, "_");
        try {
            for(auto cut_str : cut_strings) {
                bool result = true;
                if(cut_str.substr(0, 3) == "not") {
                    cut_str.erase(0, 3);
                    result = false;
                }
                const auto cut = ::analysis::Parse<SelectionCut>(cut_str);
                sub.SetCutResult(cut, result);
            }
        } catch(exception& e) {
            throw exception("Invalid event sub-category '%1%'. %2%") % str % e.message();
        }

        return sub;
    }

private:
    static size_t GetIndex(SelectionCut cut)
    {
        size_t index = static_cast<size_t>(cut);
        if(index >= MaxNumberOfCuts)
            throw exception("Cut index is out of range for cut '%1%'.") % cut;
        return index;
    }

private:
    BitsContainer presence, results;
    boost::optional<SelectionCut> last_mva_cut;
};

std::ostream& operator<<(std::ostream& os, const EventSubCategory& eventSubCategory)
{
    os << eventSubCategory.ToString();
    return os;
}

std::istream& operator>>(std::istream& is, EventSubCategory& eventSubCategory)
{
    std::string str;
    is >> str;
    eventSubCategory = EventSubCategory::Parse(str);
    return is;
}

using EventRegionSet = std::set<EventRegion>;
using EventCategorySet = std::set<EventCategory>;
using EventSubCategorySet = std::set<EventSubCategory>;
using Dataset = std::string;

} // namespace analysis
