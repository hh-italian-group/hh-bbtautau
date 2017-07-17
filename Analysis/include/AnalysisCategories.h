/*! Definition of data and event categories used in HH->bbTauTau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <boost/optional/optional.hpp>
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/Tools.h"

namespace analysis {

enum class DataCategoryType { Signal, Signal_SM, Background, DataDrivenBkg, Data };
ENUM_NAMES(DataCategoryType) = {
    { DataCategoryType::Signal, "Signal" }, { DataCategoryType::Signal_SM, "Signal_SM" },
    { DataCategoryType::Background, "Background" }, { DataCategoryType::DataDrivenBkg, "DataDrivenBkg" },
    { DataCategoryType::Data, "Data" },
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

private:
    boost::optional<size_t> n_jets, n_btag;
    boost::optional<DiscriminatorWP> btag_wp;
};

std::ostream& operator<<(std::ostream& os, const EventCategory& eventCategory)
{
    os << eventCategory.ToString();
    return os;
}

enum class SelectionCut { InsideMassWindow = 0, MVA = 1, KinematicFitConverged = 2 };
ENUM_NAMES(SelectionCut) = {
    { SelectionCut::InsideMassWindow, "InsideMassWindow" }, { SelectionCut::MVA, "MVA" },
    { SelectionCut::KinematicFitConverged, "KinematicFitConverged" }
};

struct EventSubCategory {
    using BitsContainer = unsigned long;
    static constexpr size_t MaxNumberOfCuts = std::numeric_limits<BitsContainer>::digits;
    using Bits = std::bitset<MaxNumberOfCuts>;

    static const EventSubCategory& NoCuts() { static const EventSubCategory esc; return esc; }

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

    std::string ToString() const
    {
        bool is_first = true;
        std::ostringstream s;
        for(size_t n = 0; n < MaxNumberOfCuts; ++n) {
            if(!presence[n]) continue;
            if(!is_first) {
                s << "_";
            }
            is_first = false;
            if(!results[n])
                s << "not";
            s << static_cast<SelectionCut>(n);
        }
        return s.str();
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

using EventRegionSet = std::set<EventRegion>;
using EventCategorySet = std::set<EventCategory>;
using EventSubCategorySet = std::set<EventSubCategory>;
using Dataset = std::string;

} // namespace analysis
