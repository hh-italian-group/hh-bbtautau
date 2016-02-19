/*! Classes to define uncertainties and read these definitions from a configuration file.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"

namespace analysis {
namespace limits {

enum class UncertaintyType { Normalization, Shape };

enum class UncertaintyRange { Global, Channel, Category };

enum class SampleCategory { Signal, Background, Data, StandardModel };

namespace detail {
const std::map<UncertaintyType, std::string> UncertaintyTypeNames = {
    { UncertaintyType::Normalization, "lnN" },
    { UncertaintyType::Shape, "shape" }
};

const std::map<UncertaintyRange, std::string> UncertaintyRangeNames = {
    { UncertaintyRange::Global, "global" },
    { UncertaintyRange::Channel, "channel" },
    { UncertaintyRange::Category, "category" }
};

const std::map<SampleCategory, std::string> SampleCategoryNames = {
    { SampleCategory::Signal, "signal" },
    { SampleCategory::Background, "background" },
    { SampleCategory::Data, "data" },
    { SampleCategory::StandardModel, "SM" }
};

} // namespace detail

std::ostream& operator<< (std::ostream& s, UncertaintyType t)
{
    s << detail::UncertaintyTypeNames.at(t);
    return s;
}

std::istream& operator>> (std::istream& s, UncertaintyType& t) {
    std::string name;
    s >> name;
    for(const auto& map_entry : detail::UncertaintyTypeNames) {
        if(map_entry.second == name) {
            t = map_entry.first;
            return s;
        }
    }
    throw exception("Unknown uncertainty type '") << name << "'.";
}

std::ostream& operator<< (std::ostream& s, UncertaintyRange r)
{
    s << detail::UncertaintyRangeNames.at(r);
    return s;
}

std::istream& operator>> (std::istream& s, UncertaintyRange& r) {
    std::string name;
    s >> name;
    for(const auto& map_entry : detail::UncertaintyRangeNames) {
        if(map_entry.second == name) {
            r = map_entry.first;
            return s;
        }
    }
    throw exception("Unknown uncertainty range '") << name << "'.";
}

std::ostream& operator<< (std::ostream& s, SampleCategory c)
{
    s << detail::SampleCategoryNames.at(c);
    return s;
}

std::istream& operator>> (std::istream& s, SampleCategory& c) {
    std::string name;
    s >> name;
    for(const auto& map_entry : detail::SampleCategoryNames) {
        if(map_entry.second == name) {
            c = map_entry.first;
            return s;
        }
    }
    throw exception("Unknown uncertainty range '") << name << "'.";
}

typedef std::set<std::string> SampleNameSet;
typedef std::map<SampleCategory, SampleNameSet> SampleCategoryMap;
typedef std::map<std::string, SampleNameSet> SampleSuffixMap;

struct UncertaintyDescriptor {
    std::string name;
    std::string description, source;
    std::string name_prefix, name_suffix;
    std::set<std::string> samples;
    std::map<std::string, double> sample_values;
    UncertaintyType type;
    UncertaintyRange range;
    bool calculate_value;
    double value;
    double threshold;

    UncertaintyDescriptor()
        : type(UncertaintyType::Normalization), range(UncertaintyRange::Global), calculate_value(false), value(1.),
          threshold(0.)
    {}

    std::string FullName(const std::string& channel_name, const std::string& category_name) const
    {
        static const std::string separator = "_";

        std::ostringstream s_full;
        if(name_prefix.size())
            s_full << name_prefix << separator;
        s_full << name;
        if(range == UncertaintyRange::Channel || range == UncertaintyRange::Category) {
            s_full << separator << channel_name;
            if(range == UncertaintyRange::Category)
                s_full << separator << category_name;
        }
        if(name_suffix.size())
            s_full << separator << name_suffix;
        return s_full.str();
    }

    std::string SampleList(const SampleNameSet& category_samples) const
    {
        static const std::string separator = ",";

        std::ostringstream s_list;
        bool first = true;
        for(const std::string& sample : samples) {
            if(!category_samples.count(sample)) continue;
            if(!first)
                s_list << separator;
            s_list << sample;
            first = false;
        }
        return s_list.str();
    }
};

typedef std::vector<const UncertaintyDescriptor*> UncertaintyDescriptorPtrVector;
typedef std::map<std::string, UncertaintyDescriptor> UncertaintyDescriptorMap;

class UncertaintyDescriptorCollection {
public:
    void Add(const UncertaintyDescriptor& descriptor)
    {
        if(Contains(descriptor.name))
            throw exception("Uncertainty descriptor with name '") << descriptor.name << "' already exists.";

        uncertainties[descriptor.name] = descriptor;
        orderedUncertainties.push_back(&uncertainties[descriptor.name]);
    }

    bool Contains(const std::string& descriptorName) const { return uncertainties.count(descriptorName); }
    const UncertaintyDescriptor& Get(const std::string& descriptorName) const
    {
        if(!Contains(descriptorName))
            throw exception("Uncertainty descriptor '") << descriptorName << "' not found.";
        return uncertainties.at(descriptorName);
    }

    const UncertaintyDescriptorPtrVector& GetOrderedCollection() const { return orderedUncertainties; }

private:
    UncertaintyDescriptorMap uncertainties;
    UncertaintyDescriptorPtrVector orderedUncertainties;
};


class SampleCategoryCollection {
public:
    void AddSample(SampleCategory sample_category, const std::string& sample_name)
    {
        if(all_samples.count(sample_name))
            throw exception("Sample with name '") << sample_name << "' already exist in the collection.";
        all_samples.insert(sample_name);
        all_samples.insert(detail::SampleCategoryNames.at(sample_category));
        samples[sample_category].insert(sample_name);
    }

    void AddSampleSuffix(const std::string& sample_name, const std::string& suffix)
    {
        if(!all_samples.count(sample_name))
            throw exception("Sample with name '") << sample_name << "' does not exist.";
        suffixes[sample_name].insert(suffix);
    }

    void Include(const SampleCategoryCollection& other)
    {
        all_samples.insert(other.all_samples.begin(), other.all_samples.end());
        for(const auto& sample_entry : other.samples)
            samples[sample_entry.first].insert(sample_entry.second.begin(), sample_entry.second.end());
        for(const auto& suffix_entry : other.suffixes)
            suffixes[suffix_entry.first].insert(suffix_entry.second.begin(), suffix_entry.second.end());
    }

    SampleNameSet GenerateSampleListToProcess(const std::string& base_sample_name) const
    {
        SampleNameSet result;
        SampleCategory sample_category;
        if(TryConvertToSampleCategory(base_sample_name, sample_category)) {
            const auto& sub_samples = samples.at(sample_category);
            for(const std::string& sub_sample : sub_samples)
                AddSampleWithSuffixes(sub_sample, result);
        } else {
            AddSampleWithSuffixes(base_sample_name, result);
        }
        return result;
    }

    const std::string& GetName() const { return name; }
    void SetName(const std::string& _name) { name = _name; }
    const SampleNameSet& GetAllSamples() const { return all_samples; }
    const SampleCategoryMap& GetCategorisedSamples() const { return samples; }
    const SampleSuffixMap& GetSuffixes() const { return suffixes; }

    static bool TryConvertToSampleCategory(const std::string& sample_name, SampleCategory& sample_category)
    {
        try {
            std::istringstream s_name(sample_name);
            s_name >> sample_category;
            return true;
        } catch(exception&) {}
        return false;
    }

private:
    void AddSampleWithSuffixes(const std::string& sample_name, SampleNameSet& name_set) const
    {
        if(!suffixes.count(sample_name) || !suffixes.at(sample_name).size())
            name_set.insert(sample_name);
        else {
            for(const auto& suffix : suffixes.at(sample_name))
                name_set.insert(sample_name + suffix);
        }
    }

private:
    std::string name;
    SampleNameSet all_samples;
    SampleCategoryMap samples;
    SampleSuffixMap suffixes;
};

typedef std::map<std::string, SampleCategoryCollection> SampleCategoryCollectionMap;

struct CategoryDescriptor {
    std::string name;
    std::string category_name;
    std::string channel_name;
    std::string index;
    SampleCategoryCollection samples;
};

typedef std::map<std::string, CategoryDescriptor> CategoryDescriptorMap;

class SampleCategoryCollectionReader : public ConfigEntryReader {
public:
    SampleCategoryCollectionReader(SampleCategoryCollectionMap& _output) : output(&_output) {}

    virtual void StartEntry(const std::string& name) override
    {
        if(output->count(name))
            throw exception("Samples collection with name '") << name << "' already exists.";
        current = SampleCategoryCollection();
        current.SetName(name);
    }

    virtual void EndEntry() override
    {
        (*output)[current.GetName()] = current;
    }

    virtual void ReadParameter(const std::string& param_name, const std::string& param_value) override
    {
        ReadSamples(param_name, param_value, current, *output);
    }

    static void ReadSamples(const std::string& param_name, const std::string& param_value,
                            SampleCategoryCollection& collection, const SampleCategoryCollectionMap& all_collections)
    {
        std::istringstream ss(param_value);
        ss >> std::boolalpha;
        if(param_name == "samples") {
            SampleCategory category;
            std::string s_names;
            ss >> category >> s_names;
            const auto names = ConfigReader::ParseParameterList(s_names);
            for(const std::string& name : names)
                collection.AddSample(category, name);
        } else if(param_name == "include_samples") {
            if(!all_collections.count(param_value))
                throw exception("Can't include samples from '") << param_value << "'. Samples collection not found.";
            collection.Include(all_collections.at(param_value));
        } else if(param_name == "sample_suffixes") {
            std::string sample_name, suffix_name_list;
            ss >> sample_name >> suffix_name_list;
            const auto suffixes = ConfigReader::ParseParameterList(suffix_name_list);
            for(const auto& suffix : suffixes)
                collection.AddSampleSuffix(sample_name, suffix);
        } else
            throw exception("Unsupported parameter '") << param_name << "'.";
    }

private:
    SampleCategoryCollectionMap* output;
    SampleCategoryCollection current;
};

class CategoryDescriptorReader : public ConfigEntryReader {
public:
    CategoryDescriptorReader(const SampleCategoryCollectionMap& _sampleCategoryCollections,
                             CategoryDescriptorMap& _output)
        : sampleCategoryCollections(&_sampleCategoryCollections), output(&_output) {}

    virtual void StartEntry(const std::string& name) override
    {
        if(output->count(name))
            throw exception("Category descriptor with name '") << name << "' already exists.";

        current = CategoryDescriptor();
        current.name = name;
    }

    virtual void EndEntry() override
    {
        (*output)[current.name] = current;
    }

    virtual void ReadParameter(const std::string& param_name, const std::string& param_value) override
    {
        std::istringstream ss(param_value);
        ss >> std::boolalpha;
        if(param_name == "samples" || param_name == "include_samples") {
            SampleCategoryCollectionReader::ReadSamples(param_name, param_value, current.samples,
                                                        *sampleCategoryCollections);
        } else if(param_name == "category_name") {
            ss >> current.category_name;
        } else if(param_name == "channel_name") {
            ss >> current.channel_name;
        } else if(param_name == "index") {
            ss >> current.index;
        } else
            throw exception("Unsupported parameter '") << param_name << "'.";
    }

private:
    const SampleCategoryCollectionMap* sampleCategoryCollections;
    CategoryDescriptorMap* output;
    CategoryDescriptor current;
};

class UncertaintyDescriptorReader : public ConfigEntryReader {
public:
    UncertaintyDescriptorReader(UncertaintyDescriptorCollection& _output) : output(&_output) {}

    virtual void StartEntry(const std::string& name) override
    {
        if(output->Contains(name))
            throw exception("Uncertainty descriptor with name '") << name << "' already exists.";

        current = UncertaintyDescriptor();
        current.name = name;
    }

    virtual void EndEntry() override
    {
        output->Add(current);
    }

    virtual void ReadParameter(const std::string& param_name, const std::string& param_value) override
    {
        std::istringstream ss(param_value);
        ss >> std::boolalpha;
        if(param_name == "description") {
            std::string description;
            ss >> description;
            current.description += "\n" + description;
        } else if(param_name == "source") {
            std::string source;
            ss >> current.source;
            current.source += "\n" + source;
        } else if(param_name == "name_prefix") {
            ss >> current.name_prefix;
        } else if(param_name == "name_suffix") {
            ss >> current.name_suffix;
        } else if(param_name == "samples") {
            std::string s_names;
            ss >> s_names;
            const auto names = ConfigReader::ParseParameterList(s_names);
            current.samples.insert(names.begin(), names.end());
        } else if(param_name == "type") {
            ss >> current.type;
        } else if(param_name == "range") {
            ss >> current.range;
        } else if(param_name == "calc_value") {
            ss >> current.calculate_value;
        } else if(param_name == "value") {
            ss >> current.value;
        } else if(param_name == "threshold") {
            ss >> current.threshold;
        } else if(param_name == "sample_value") {
            std::string sample_name;
            double value;
            ss >> sample_name >> value;
            if(!current.samples.count(sample_name))
                throw exception("Sample '") << sample_name << "' not listed in the samples list.";
            current.sample_values[sample_name] = value;
        } else
            throw exception("Unsupported parameter '") << param_name << "'.";
    }

private:
    UncertaintyDescriptorCollection* output;
    UncertaintyDescriptor current;
};

struct UncertaintyInterval;
std::ostream& operator<< (std::ostream& s, const UncertaintyInterval& i);

struct UncertaintyInterval {
private:
    typedef PhysicalValue (UncertaintyInterval::*ValuePtr);

public:
    PhysicalValue down, up;

    UncertaintyInterval() : down(PhysicalValue::One), up(PhysicalValue::One) {}
    explicit UncertaintyInterval(const PhysicalValue& unc)
        : down(PhysicalValue::One - unc), up(PhysicalValue::One + unc) {}
    UncertaintyInterval(const PhysicalValue& _down, const PhysicalValue& _up) : down(_down), up(_up) {}

    static UncertaintyInterval WeightedAverage(const std::vector<UncertaintyInterval>& intervals)
    {
        const auto down_collection = CollectValues(intervals, &UncertaintyInterval::down);
        const auto up_collection = CollectValues(intervals, &UncertaintyInterval::up);
        const PhysicalValue new_down = PhysicalValue::WeightedAverage(down_collection);
        const PhysicalValue new_up = PhysicalValue::WeightedAverage(up_collection);
        return UncertaintyInterval(new_down, new_up);
    }

private:
    static std::vector<PhysicalValue> CollectValues(const std::vector<UncertaintyInterval>& intervals,
                                                    ValuePtr value_ptr)
    {
        std::vector<PhysicalValue> result;
        for(const UncertaintyInterval& unc : intervals)
            result.push_back(unc.*value_ptr);
        return result;
    }
};

std::ostream& operator<< (std::ostream& s, const UncertaintyInterval& i)
{
    s << "(" << i.down << ", " << i.up << ")";
    return s;
}

} // namespace limits
} // namespace analysis
