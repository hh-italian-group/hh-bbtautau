/*! Parameters, benchmarks and kinematic parameterizations used in HH effective-field theory.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/TextIO.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"

namespace analysis {

namespace NonResHH_EFT {

struct Point {
    double kl{1}, kt{1}, c2{0}, cg{0}, c2g{0};

    static const std::vector<std::string> ParameterNames()
    {
        const std::vector<std::string> names = { "kl", "kt", "c2", "cg", "c2g" };
        return names;
    }

    static const Point& SM_Point() { static const Point p(1, 1, 0, 0, 0); return p; }

    Point() {}
    Point(double _kl, double _kt, double _c2, double _cg, double _c2g) :
        kl(_kl), kt(_kt), c2(_c2), cg(_cg), c2g(_c2g) {}
    Point(const std::vector<double>& values)
    {
        if(values.size() != 5)
            throw exception("Invalid numer of dimensions for non-resonat HH EFT parameters");
        kl = values[0];
        kt = values[1];
        c2 = values[2];
        cg = values[3];
        c2g = values[4];
    }

    std::vector<double> ToVector() const { return { kl, kt, c2, cg, c2g }; }
};

inline std::ostream& operator<<(std::ostream& s, const Point& p)
{
    s << p.kl << " " << p.kt << " " << p.c2 << " " << p.cg << " " << p.c2g;
    return s;
}

inline std::istream& operator>>(std::istream& s, Point& p)
{
    try {
        auto values = ReadValueList(s, 5, true, " /t", true);
        p.kl = Parse<double>(values.at(0));
        p.kt = Parse<double>(values.at(1));
        p.c2 = Parse<double>(values.at(2));
        p.cg = Parse<double>(values.at(3));
        p.c2g = Parse<double>(values.at(4));
    } catch(std::exception& e) {
        throw exception("Invalid EFT point. %1%") % e.what();
    }
    return s;
}

struct Benchmark {
    Point point;
    int id;
    std::string name;

    Benchmark() {}
    Benchmark(const Point _point, int _id, const std::string& _name = "") :
        point(_point), id(_id), name(_name)
    {
        if(name.empty())
            name = std::to_string(id);
    }
};

class BenchmarkCollection {
public:
    using BenchmarkMap = std::map<int, Benchmark>;
    using const_iterator = BenchmarkMap::const_iterator;

    static const BenchmarkCollection& Benchmarks()
    {
        static BenchmarkCollection col;
        if(!col.size()) {
            col.insert({ 7.5, 1., -1., 0., 0. }, 1);
            col.insert({ 1., 1., 0.5, -0.8, 0.6 }, 2);
            col.insert({ 1., 1., -1.5, 0., -0.8 }, 3);
            col.insert({ -3.5, 1.5, -3., 0., 0. }, 4);
            col.insert({ 1., 1., 0., 0.8, -1. }, 5);
            col.insert({ 2.4, 1., 0., 0.2, -0.2 }, 6);
            col.insert({ 5., 1., 0., 0.2, -0.2 }, 7);
            col.insert({ 15., 1., 0., -1., 1. }, 8);
            col.insert({ 1., 1., 1., -0.6, 0.6 }, 9);
            col.insert({ 10, 1.5, -1.0, 0., 0. }, 10);
            col.insert({ 2.4, 1., 0., 1., -1. }, 11);
            col.insert({ 15., 1., 1., 0., 0. }, 12);
            col.insert({ 0., 1., 0., 0., 0. }, 13, "box");
            col.insert(Point::SM_Point(), 14, "SM");
        }
        return col;
    }

    static BenchmarkCollection CreateScanBenchmarks(const std::vector<std::vector<double>>& values)
    {
        if(values.size() != 5)
            throw exception("Invalid numer of dimensions for non-resonat HH EFT parameters grid.");
        std::vector<size_t> limits;
        for(const auto& v : values)
            limits.push_back(v.size());
        Grid_ND grid(limits);
        BenchmarkCollection col;
        int id = 0;
        for(const auto& grid_point : grid) {
            std::vector<double> p_vector(values.size());
            for(size_t n = 0; n < values.size(); ++n)
                p_vector[n] = values.at(n).at(grid_point.at(n));
            const Point point(p_vector);
            std::ostringstream p_name;
            for(size_t n = 0; n < limits.size(); ++n) {
                if(limits.at(n) > 1) {
                    p_name << Point::ParameterNames().at(n) << p_vector.at(n);
                }
            }
            col.insert(point, id++, p_name.str());
        }
        return col;
    }

    void insert(const Benchmark& benchmark)
    {
        if(benchmarks.count(benchmark.id))
            throw exception("Benchmark with id = %1% already present in the collection.") % benchmark.id;
        if(names.count(benchmark.name))
            throw exception("Benchmark with name = '%1%' already present in the collection.") % benchmark.name;
        benchmarks[benchmark.id] = benchmark;
        names[benchmark.name] = benchmark.id;
        ids.push_back(benchmark.id);
    }

    void insert(const Point& point, int id, const std::string& name = "")
    {
        insert(Benchmark(point, id, name));
    }

    const Benchmark& at(int id) const
    {
        auto iter = benchmarks.find(id);
        if(iter == benchmarks.end())
            throw exception("Benchmark with id = %1% does not present in the collection.") % id;
        return iter->second;
    }

    const Benchmark& at(const std::string& name) const
    {
        auto iter = names.find(name);
        if(iter == names.end())
            throw exception("Benchmark with name = '%1%' does not present in the collection.") % name;
        return at(iter->second);
    }

    const Benchmark& at_index(size_t index) const
    {
        if(index >= ids.size())
            throw exception("Benchark index is out of range.");
        return at(ids.at(index));
    }

    const_iterator begin() const { return benchmarks.begin(); }
    const_iterator end() const { return benchmarks.end(); }
    size_t size() const { return benchmarks.size(); }

private:
    std::map<int, Benchmark> benchmarks;
    std::map<std::string, int> names;
    std::vector<int> ids;
};

// Definition taken from here:
// https://github.com/cms-hh/HHStatAnalysis/blob/master/AnalyticalModels/python/NonResonantModel.py
struct GF_Parameterization {
    static constexpr size_t n_params = 15;
    using Array = std::array<double, n_params>;

    static const GF_Parameterization& Denominator_13TeV()
    {
        static const GF_Parameterization p({
            2.09078, 10.1517, 0.282307, 0.101205, 1.33191, -8.51168, -1.37309, 2.82636, 1.45767, -4.91761,
            -0.675197, 1.86189, 0.321422, -0.836276, -0.568156
        });
        return p;
    }

    GF_Parameterization() { a_err.fill(0); }
    GF_Parameterization(const Array& _a) : a(_a) { a_err.fill(0); }

    double Evaluate(const Point& p) const
    {
        const double kl = p.kl, kl_2 = std::pow(kl, 2);
        const double kt = p.kt, kt_2 = std::pow(kt, 2), kt_4 = std::pow(kt, 4);
        const double c2 = p.c2, c2_2 = std::pow(c2, 2);
        const double cg = p.cg, cg_2 = std::pow(cg, 2);
        const double c2g = p.c2g, c2g_2 = std::pow(c2g, 2);
        return a[0] * kt_4 + a[1] * c2_2 + (a[2] * kt_2 + a[3] * cg_2 ) * kl_2 + a[4] * c2g_2
                + (a[5] * c2 + a[6] * kt * kl) * kt_2 + (a[7] * kt * kl + a[8] * cg * kl) * c2 + a[9] * c2 * c2g
                + (a[10] * cg * kl + a[11] * c2g) * kt_2 + (a[12] * kl * cg + a[13] * c2g) * kt * kl
                + a[14] * cg * c2g * kl;
    }

    Array a, a_err;
};

inline std::ostream& operator<<(std::ostream& s, const GF_Parameterization& p)
{
    s << p.a[0];
    for(size_t n = 1; n < GF_Parameterization::n_params; ++n)
        s << " " << p.a[n];
    for(size_t n = 0; n < GF_Parameterization::n_params; ++n)
        s << " " << p.a_err[n];
    return s;
}

inline std::istream& operator>>(std::istream& s, GF_Parameterization& p)
{
    try {
        auto values = ReadValueList(s, GF_Parameterization::n_params * 2, true, " /t", true);
        for(size_t n = 0; n < GF_Parameterization::n_params; ++n) {
            p.a[n] = Parse<double>(values.at(n));
            p.a_err[n] = Parse<double>(values.at(n + GF_Parameterization::n_params));
        }
    } catch(std::exception& e) {
        throw exception("Invalid EFT parameterization definition. %1%") % e.what();
    }
    return s;
}

} // namespace NonResHH_EFT
} // namespace analysis

