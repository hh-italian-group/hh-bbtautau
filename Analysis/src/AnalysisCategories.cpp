/*! Definition of data and event categories used in HH->bbTauTau analysis.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "hh-bbtautau/Analysis/include/AnalysisCategories.h"
#include <boost/bimap.hpp>

namespace analysis {
const std::unique_ptr<boost::bimap<std::string, EventRegion>> EventRegion::predefined_regions =
    std::make_unique<boost::bimap<std::string, EventRegion>>();
const std::unique_ptr<boost::optional<EventRegion>> EventRegion::_OS_Isolated =
    std::make_unique<boost::optional<EventRegion>>();
const std::unique_ptr<boost::optional<EventRegion>> EventRegion::_OS_AntiIsolated =
    std::make_unique<boost::optional<EventRegion>>();
const std::unique_ptr<boost::optional<EventRegion>> EventRegion::_SS_Isolated =
    std::make_unique<boost::optional<EventRegion>>();
const std::unique_ptr<boost::optional<EventRegion>> EventRegion::_SS_LooseIsolated =
    std::make_unique<boost::optional<EventRegion>>();
const std::unique_ptr<boost::optional<EventRegion>> EventRegion::_SS_AntiIsolated =
    std::make_unique<boost::optional<EventRegion>>();

void EventRegion::Initialize(DiscriminatorWP iso_lower, DiscriminatorWP anti_iso_lower, DiscriminatorWP anti_iso_upper)
{
    *_OS_Isolated = EventRegion().SetCharge(true).SetLowerIso(iso_lower);
    *_OS_AntiIsolated = EventRegion().SetCharge(true).SetLowerIso(anti_iso_lower).SetUpperIso(anti_iso_upper);
    *_SS_Isolated = EventRegion().SetCharge(false).SetLowerIso(iso_lower);
    *_SS_LooseIsolated = EventRegion().SetCharge(false).SetLowerIso(DiscriminatorWP::Loose);
    *_SS_AntiIsolated = EventRegion().SetCharge(false).SetLowerIso(anti_iso_lower).SetUpperIso(anti_iso_upper);

    predefined_regions->clear();
    predefined_regions->insert({ "Unknown", EventRegion::Unknown() });
    predefined_regions->insert({ "OS_Isolated", EventRegion::OS_Isolated() });
    predefined_regions->insert({ "OS_AntiIsolated", EventRegion::OS_AntiIsolated() });
    predefined_regions->insert({ "SS_LooseIsolated", EventRegion::SS_LooseIsolated() });
    predefined_regions->insert({ "SS_Isolated", EventRegion::SS_Isolated() });
    predefined_regions->insert({ "SS_AntiIsolated", EventRegion::SS_AntiIsolated() });
    predefined_regions->insert({ "SignalRegion", EventRegion::SignalRegion() });

}

const boost::bimap<std::string, analysis::EventRegion>& EventRegion::EventRegionMapToString()
{
    if(predefined_regions->size() == 0)
        throw exception("predefined regions maps is not initialized");
    return *predefined_regions;
}

const EventRegion& EventRegion::Unknown()
{
    static const EventRegion er;
    return er;
}

const EventRegion& EventRegion::OS_Isolated()
{
    if(!_OS_Isolated->is_initialized())
        throw exception("Event regions (OS_Isolated) are not initialized");
    return **_OS_Isolated;
}

const EventRegion& EventRegion::OS_AntiIsolated()
{
    if(!_OS_AntiIsolated->is_initialized())
        throw exception("Event regions (OS_AntiIsolated) are not initialized");
    return **_OS_AntiIsolated;
}

const EventRegion& EventRegion::SS_Isolated()
{
    if(!_SS_Isolated->is_initialized())
        throw exception("Event regions (SS_Isolated) are not initialized");
    return **_SS_Isolated;
}

const EventRegion& EventRegion::SS_LooseIsolated()
{
    if(!_SS_LooseIsolated->is_initialized())
        throw exception("Event regions (SS_LooseIsolated) not initialized");
    return **_SS_LooseIsolated;
}

const EventRegion& EventRegion::SS_AntiIsolated()
{
    if(!_SS_AntiIsolated->is_initialized())
        throw exception("Event regions (SS_AntiIsolated) are not initialized");
    return **_SS_AntiIsolated;
}

const EventRegion& EventRegion::SignalRegion()
{
    return OS_Isolated();
}

EventRegion& EventRegion::SetCharge(bool _os)
{
    os = _os; return *this;
}

EventRegion& EventRegion::SetLowerIso(DiscriminatorWP wp)
{
    if(HasUpperIso() && iso_upper <= wp)
        throw exception("HasUpperIso - Iso Upper limit is not greater than Iso lower limit");
    iso_lower = wp;
    return *this;
}

EventRegion& EventRegion::SetUpperIso(DiscriminatorWP wp)
{
    if(HasLowerIso() && wp <= iso_lower)
        throw exception("HasLowerIso - Iso Upper limit is not greater than Iso lower limit");
    iso_upper = wp;
    return *this;
}

bool EventRegion::HasCharge() const { return os.is_initialized(); }
bool EventRegion::HasLowerIso() const { return iso_lower.is_initialized(); }
bool EventRegion::HasUpperIso() const { return iso_upper.is_initialized(); }

DiscriminatorWP EventRegion::GetLowerIso() const
{
    if(!HasLowerIso())
        throw exception("Lower isolation bound not set.");
    return *iso_lower;
}

DiscriminatorWP EventRegion::GetUpperIso() const
{
    if(!HasUpperIso())
        throw exception("Upper isolation bound not set.");
    return *iso_upper;
}

bool EventRegion::GetCharge() const
{
    if(!HasCharge())
        throw exception("Charge info not set.");
    return *os;
}

bool EventRegion::Implies(const EventRegion& other) const
{
    if(other.HasCharge() && (!HasCharge() || GetCharge() != other.GetCharge())) return false;
    if(other.HasLowerIso() && (!HasLowerIso() || GetLowerIso() < other.GetLowerIso())) return false;
    return !other.HasUpperIso() || (HasUpperIso() && GetUpperIso() <= other.GetUpperIso());
}

bool EventRegion::operator ==(const EventRegion& er) const
{
    return os == er.os && iso_lower == er.iso_lower && iso_upper == er.iso_upper;
}

bool EventRegion::operator !=(const EventRegion& er) const { return !(*this == er); }

bool EventRegion::operator <(const EventRegion& er) const
{
    if(os != er.os) return os < er.os;
    if(iso_lower != er.iso_lower) return iso_lower < er.iso_lower;
    return iso_upper < er.iso_upper;
}

std::string EventRegion::ToString() const
{
    if(*this == Unknown()) return "Unknown";
    if(!EventRegionMapToString().right.count(*this))
        throw exception("Unknown EventRegion. No conversion to String");
    return EventRegionMapToString().right.at(*this);
}

EventRegion EventRegion::Parse(const std::string& str)
{
    if(!EventRegionMapToString().left.count(str))
        throw exception("Unknown EventRegion = '%1%'.") % str;
    return EventRegionMapToString().left.at(str);
}

std::istream& operator>>(std::istream& s, EventRegion& eventRegion)
{
    std::string str;
    s >> str;
    eventRegion = EventRegion::Parse(str);
    return s;
}

std::ostream& operator<<(std::ostream& os, const EventRegion& eventRegion)
{
    os << eventRegion.ToString();
    return os;
}

const EventCategory& EventCategory::Inclusive() { static const EventCategory ec; return ec; }

EventCategory::EventCategory(size_t _n_jets) : n_jets(_n_jets) {}

EventCategory::EventCategory(size_t _n_jets, size_t _n_btag, bool _strict_n_btag, DiscriminatorWP _btag_wp) :
        n_jets(_n_jets), n_btag(_n_btag), strict_n_btag(_strict_n_btag), btag_wp(_btag_wp)
{
    if(n_btag > n_jets)
        throw exception("Number of btag can't be greater than number of jets");
}

EventCategory::EventCategory(size_t _n_jets, size_t _n_btag, bool _strict_n_btag, DiscriminatorWP _btag_wp,
                             bool _boosted) :
    n_jets(_n_jets), n_btag(_n_btag), strict_n_btag(_strict_n_btag), btag_wp(_btag_wp), boosted(_boosted)
{
    if(n_btag > n_jets)
        throw exception("Number of btag can't be greater than number of jets");
}

EventCategory::EventCategory(size_t _n_jets, size_t _n_btag, bool _strict_n_btag, DiscriminatorWP _btag_wp,
                             boost::optional<bool> _boosted, bool _is_VBF):
    n_jets(_n_jets), n_btag(_n_btag), strict_n_btag(_strict_n_btag), btag_wp(_btag_wp), boosted(_boosted), is_VBF(_is_VBF)
{
    if(n_btag > n_jets)
        throw exception("Number of btag can't be greater than number of jets");
}

bool EventCategory::HasJetConstraint() const { return n_jets.is_initialized(); }
size_t EventCategory::N_jets() const
{
    if(!HasJetConstraint())
        throw exception("Jet constraint is not defined.");
    return *n_jets;
}

bool EventCategory::HasBtagConstraint() const { return n_btag.is_initialized(); }
size_t EventCategory::N_btag() const
{
    if(!HasBtagConstraint())
        throw exception("Btag constraint is not defined.");
    return *n_btag;
}
DiscriminatorWP EventCategory::BtagWP() const
{
    if(!HasBtagConstraint())
        throw exception("Btag constraint is not defined.");
    return *btag_wp;
}

bool EventCategory::IsN_btagStrict() const
{
    if(!strict_n_btag.is_initialized())
        throw exception("Strict Btag constraint is not defined.");
    return *strict_n_btag;
}

bool EventCategory::HasBoostConstraint() const { return boosted.is_initialized(); }
bool EventCategory::IsBoosted() const
{
    if(!HasBoostConstraint())
        throw exception("Boost constraint is not defined.");
    return *boosted;
}

bool EventCategory::HasVBFConstraint() const { return is_VBF.is_initialized(); }
bool EventCategory::isVBF() const
{
    if(!HasVBFConstraint())
        throw exception("VBF constraint is not defined.");
    return *is_VBF;
}

bool EventCategory::operator ==(const EventCategory& ec) const
{
    return n_jets == ec.n_jets && n_btag == ec.n_btag && strict_n_btag == ec.strict_n_btag && btag_wp == ec.btag_wp && boosted == ec.boosted && is_VBF == ec.is_VBF;
}
bool EventCategory::operator !=(const EventCategory& ec) const { return !(*this == ec); }
bool EventCategory::operator <(const EventCategory& ec) const
{
    if(n_jets != ec.n_jets) return n_jets < ec.n_jets;
    if(n_btag != ec.n_btag) return n_btag < ec.n_btag;
    if (strict_n_btag != ec.strict_n_btag) return strict_n_btag < ec.strict_n_btag;
    if(btag_wp != ec.btag_wp) return btag_wp < ec.btag_wp;
    if (boosted != ec.boosted) return boosted < ec.boosted;
    return is_VBF < ec.is_VBF;
}


std::string EventCategory::ToString() const
{
    if(*this == Inclusive()) return "Inclusive";
    std::ostringstream s;
    s << *n_jets << "j";
    if(HasBtagConstraint()) {
        s << *n_btag;
        if(*btag_wp != DiscriminatorWP::Medium)
            s << __DiscriminatorWP_short_names.EnumToString(*btag_wp);
        s << "b";
        const std::string stricbtag_str = !IsN_btagStrict() ? "+" : "";
        s << stricbtag_str;
    }
    if(HasBoostConstraint()) {
        const std::string boosted_str = IsBoosted() ? "B" : "R";
        s << boosted_str;
    }
    if(HasVBFConstraint()) {
        const std::string VBF_str = isVBF() ? "_VBF" : "_noVBF";
        s << VBF_str;
    }
    return s.str();
}

EventCategory EventCategory::Parse(const std::string& str)
{
    static const std::string numbers = "0123456789";
    static const std::string jets_suffix = "j", btag_suffix = "b";
    static const std::map<char, bool> boosted_suffix = { { 'R', false }, { 'B', true } };
    static const std::map<std::string, bool> VBF_suffix = { { "noVBF", false }, { "VBF", true } };

    if(str == "Inclusive") return Inclusive();
    try {
        const size_t jet_pos = str.find_first_not_of(numbers);
        if(jet_pos == std::string::npos && str.substr(jet_pos, jets_suffix.size()) != "j")
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
        const size_t strictbtag_pos = btag_pos + btag_suffix.size();
        const char strictbtag_flag = strictbtag_pos == str.size() ? ' ' : str.at(strictbtag_pos);
        const bool is_strictbtag = strictbtag_flag != '+'?  true : false;
        const size_t boosted_pos = strictbtag_flag == '+' ? strictbtag_pos + 1 : strictbtag_pos;
        if(str.size() == boosted_pos)
            return EventCategory(n_jets, n_btag, is_strictbtag, btag_wp);
        const char boosted_flag = str.at(boosted_pos);
        if(!boosted_suffix.count(boosted_flag) && boosted_flag!='_')
            throw exception("");
        boost::optional<bool> is_boosted = boost::make_optional(false, false);
        bool boosted = false;
        if (boosted_suffix.count(boosted_flag)){
            is_boosted = boosted_suffix.at(boosted_flag);
            boosted = boosted_suffix.at(boosted_flag);
        }
        if(str.size() == boosted_pos+1)
            return EventCategory(n_jets, n_btag, is_strictbtag, btag_wp, boosted);

        const size_t isVBF_pos = str.find("_")+1;
        const std::string isVBF_flag = str.substr(isVBF_pos);
        if(!VBF_suffix.count(isVBF_flag))
            throw exception("");
        const bool is_VBF = VBF_suffix.at(isVBF_flag);
        return EventCategory(n_jets, n_btag, is_strictbtag, btag_wp, is_boosted, is_VBF);
    }catch(exception& e) {
        throw exception("Invalid EventCategory '%1%'. %2%") % str % e.message();
    }
}

bool EventCategory::Contains(size_t num_jets, const std::map<DiscriminatorWP, size_t>& num_btag, bool is_vbf,
                             bool is_boosted) const
{
    if(btag_wp && !num_btag.count(*btag_wp))
        throw exception("The btag_wp, is not defined") ;

    return (!n_jets || num_jets >= *n_jets) && (!n_btag
                        || (*strict_n_btag ? (num_btag.at(*btag_wp) == n_btag) : (num_btag.at(*btag_wp) >= *n_btag)))
                        && (!is_VBF || is_vbf == *is_VBF)
                        && (!boosted || is_boosted == *boosted);
}

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


const EventSubCategory& EventSubCategory::NoCuts() { static const EventSubCategory esc; return esc; }

EventSubCategory::EventSubCategory() : presence(0), results(0) {}

bool EventSubCategory::HasCut(SelectionCut cut) const
{
    const BitsContainer mask = BitsContainer(1) << GetIndex(cut);
    return (presence & mask) != BitsContainer(0);
}

bool EventSubCategory::Passed(SelectionCut cut) const
{
    if(!HasCut(cut))
        throw exception("Cut '%1%' is not defined.") % cut;
    const BitsContainer mask = BitsContainer(1) << GetIndex(cut);
    return (results & mask) != BitsContainer(0);
}
bool EventSubCategory::Failed(SelectionCut cut) const { return !Passed(cut); }

EventSubCategory& EventSubCategory::SetCutResult(SelectionCut cut, bool result)
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

EventSubCategory::BitsContainer EventSubCategory::GetPresenceBits() const { return presence; }
EventSubCategory::BitsContainer EventSubCategory::GetResultBits() const { return results; }

bool EventSubCategory::operator ==(const EventSubCategory& sc) const
{
    return GetPresenceBits() == sc.GetPresenceBits() && GetResultBits() == sc.GetResultBits();
}
bool EventSubCategory::operator !=(const EventSubCategory& sc) const { return !(*this == sc); }
bool EventSubCategory::operator <(const EventSubCategory& sc) const
{
    if(GetPresenceBits() != sc.GetPresenceBits()) return GetPresenceBits() < sc.GetPresenceBits();
    return GetResultBits() < sc.GetResultBits();
}

bool EventSubCategory::Implies(const EventSubCategory& sc) const
{
    const BitsContainer pres_a = GetPresenceBits(), pres_b = sc.GetPresenceBits();
    if(((pres_a ^ pres_b) & pres_b) != BitsContainer(0)) return false;
    const BitsContainer res_a = GetResultBits(), res_b = sc.GetResultBits();
    return (res_a & pres_b) == res_b;
}

bool EventSubCategory::TryGetLastMvaCut(SelectionCut& cut) const
{
    if(!last_mva_cut.is_initialized()) return false;
    cut = *last_mva_cut;
    return true;
}

std::string EventSubCategory::ToString(const std::map<SelectionCut, std::string>& sel_aliases) const
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

EventSubCategory EventSubCategory::Parse(const std::string& str)
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

size_t EventSubCategory::GetIndex(SelectionCut cut)
{
    size_t index = static_cast<size_t>(cut);
    if(index >= MaxNumberOfCuts)
        throw exception("Cut index is out of range for cut '%1%'.") % cut;
    return index;
}

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

} // namespace analysis
