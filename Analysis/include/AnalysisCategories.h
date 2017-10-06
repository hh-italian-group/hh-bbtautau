/*! Definition of data and event categories used in HH->bbTauTau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <boost/optional/optional.hpp>
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/Tools.h"
#include "AnalysisTools/Core/include/TextIO.h"

namespace analysis {

enum class SampleType { Data, MC, DY, QCD };
ENUM_NAMES(SampleType) = {
    { SampleType::Data, "Data" }, { SampleType::MC, "MC" }, { SampleType::DY, "DY" }, { SampleType::QCD, "QCD" }
};

struct EventRegion {
    static const EventRegion& Unknown() { static const EventRegion er; return er; }
    static const EventRegion& OS_Isolated() { static const EventRegion er(true, true); return er; }
    static const EventRegion& OS_AntiIsolated() { static const EventRegion er(true, false); return er; }
    static const EventRegion& SS_Isolated() { static const EventRegion er(false, true); return er; }
    static const EventRegion& SS_AntiIsolated() { static const EventRegion er(false, true); return er; }
    static const EventRegion& SignalRegion() { return OS_Isolated(); }

    EventRegion() {}
    EventRegion(bool _os) : os(_os) {}
    EventRegion(bool _os, bool _iso) : os(_os), iso(_iso) {}

    bool OS() const { return os.is_initialized() && *os; }
    bool SS() const { return os.is_initialized() && !*os; }
    bool Iso() const { return iso.is_initialized() && *iso; }
    bool AntiIso() const { return iso.is_initialized() && !*iso; }

    bool operator ==(const EventRegion& er) const { return os == er.os && iso == er.iso; }
    bool operator !=(const EventRegion& er) const { return !(*this == er); }
    bool operator <(const EventRegion& er) const
    {
        if(os != er.os) return os < er.os;
        return iso < er.iso;
    }

    std::string ToString() const
    {
        if(*this == Unknown()) return "Unknown";
        std::ostringstream s;
        s << SignPairStr(*os);
        if(iso.is_initialized())
            s << "_" << IsoStr(*iso);
        return s.str();
    }

private:
    boost::optional<bool> os, iso;

    static std::string SignPairStr(bool sign_pair) { return sign_pair ? "OS" : "SS"; }
    static std::string IsoStr(bool iso) { return iso ? "Isolated" : "AntiIsolated"; }
};

std::ostream& operator<<(std::ostream& os, const EventRegion& eventRegion)
{
    os << eventRegion.ToString();
    return os;
}

struct EventCategory {
    static const EventCategory& Inclusive() { static const EventCategory ec; return ec; }
    static const EventCategory& TwoJets_Inclusive() { static const EventCategory ec(2); return ec; }
    static const EventCategory& TwoJets_ZeroBtag()
        { static const EventCategory ec(2, 0, DiscriminatorWP::Medium); return ec; }
    static const EventCategory& TwoJets_OneBtag()
        { static const EventCategory ec(2, 1, DiscriminatorWP::Medium); return ec; }
    static const EventCategory& TwoJets_TwoBtag()
        { static const EventCategory ec(2, 2, DiscriminatorWP::Medium); return ec; }
    static const EventCategory& TwoJets_ZeroLooseBtag()
        { static const EventCategory ec(2, 0, DiscriminatorWP::Loose); return ec; }
    static const EventCategory& TwoJets_OneLooseBtag()
        { static const EventCategory ec(2, 1, DiscriminatorWP::Loose); return ec; }
    static const EventCategory& TwoJets_TwoLooseBtag()
        { static const EventCategory ec(2, 2, DiscriminatorWP::Loose); return ec; }

    EventCategory() {}
    explicit EventCategory(size_t _n_jets) : n_jets(_n_jets) {}
    EventCategory(size_t _n_jets, size_t _n_btag, DiscriminatorWP _btag_wp) :
        n_jets(_n_jets), n_btag(_n_btag), btag_wp(_btag_wp)
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

    bool operator ==(const EventCategory& ec) const
    {
        return n_jets == ec.n_jets && n_btag == ec.n_btag && btag_wp == ec.btag_wp;
    }
    bool operator !=(const EventCategory& ec) const { return !(*this == ec); }
    bool operator <(const EventCategory& ec) const
    {
        if(n_jets != ec.n_jets) return n_jets < ec.n_jets;
        if(n_btag != ec.n_btag) return n_btag < ec.n_btag;
        return btag_wp < ec.btag_wp;
    }

    std::string ToString() const
    {
        if(*this == Inclusive()) return "Inclusive";
        std::ostringstream s;
        s << *n_jets << "jets";
        if(HasBtagConstraint()) {
            s << *n_btag;
            if(*btag_wp != DiscriminatorWP::Medium)
                s << *btag_wp;
            s << "btag";
        }
        return s.str();
    }

    static EventCategory Parse(const std::string& str)
    {
        static const std::string numbers = "0123456789";
        static const std::string jets_suffix = "jets", btag_suffix = "btag";

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
                    : ::analysis::Parse<DiscriminatorWP>(str.substr(btag_wp_pos, btag_pos - btag_wp_pos));
            return EventCategory(n_jets, n_btag, btag_wp);
        }catch(exception& e) {
            throw exception("Invalid EventCategory '%1%'. %2%") % str % e.message();
        }
    }

private:
    boost::optional<size_t> n_jets, n_btag;
    boost::optional<DiscriminatorWP> btag_wp;
};

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

enum class SelectionCut { InsideMassWindow = 0, KinematicFitConverged = 1,
                          MVA_CUT_LIST(2, 100) MVA_first = MVA0, MVA_last = MVA99 };

#undef MVA_CUT_LIST
#undef DECL_MVA_SEL

namespace detail {
inline std::map<SelectionCut, std::string> CreateSelectionCutNames()
{
    std::map<SelectionCut, std::string> names;
    names[SelectionCut::InsideMassWindow] = "InsideMassWindow";
    names[SelectionCut::KinematicFitConverged] = "KinematicFitConverged";
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
    using BitsContainer = unsigned long long;
    static constexpr size_t MaxNumberOfCuts = std::numeric_limits<BitsContainer>::digits;
    using Bits = std::bitset<MaxNumberOfCuts>;

    static const EventSubCategory& NoCuts() { static const EventSubCategory esc; return esc; }

    EventSubCategory() {}

    bool HasCut(SelectionCut cut) const { return presence[GetIndex(cut)]; }
    bool Passed(SelectionCut cut) const
    {
        if(!HasCut(cut))
            throw exception("Cut '%1%' is not defined.") % cut;
        return results[GetIndex(cut)];
    }
    bool Failed(SelectionCut cut) const { return !Passed(cut); }
    EventSubCategory& SetCutResult(SelectionCut cut, bool result)
    {
        if(HasCut(cut))
            throw exception("Cut '%1%' is aready defined.") % cut;
        results[GetIndex(cut)] = result;
        presence[GetIndex(cut)] = true;
        return *this;
    }

    BitsContainer GetPresenceBits() const { return presence.to_ulong(); }
    BitsContainer GetResultBits() const { return results.to_ulong(); }

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
        if((pres_a ^ pres_b) & pres_b) return false;
        const BitsContainer res_a = GetResultBits(), res_b = sc.GetResultBits();
        return (res_a & pres_b) == res_b;
    }

    std::string ToString() const
    {
        std::ostringstream s;
        for(size_t n = 0; n < MaxNumberOfCuts; ++n) {
            if(!presence[n]) continue;
            if(!results[n]) s << "not";
            s << static_cast<SelectionCut>(n) << "_";
        }
        std::string str = s.str();
        if(str.size())
            str.erase(str.size() - 1);
        return str;
    }

    static EventSubCategory Parse(const std::string& str)
    {
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
    Bits presence, results;
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
