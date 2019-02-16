/*! Parameters, benchmarks and kinematic parameterizations used in HH effective-field theory.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <array>
#include <vector>
#include <map>
#include <istream>
#include <ostream>

namespace analysis {
namespace NonResHH_EFT {

struct Point {
    double kl{1}, kt{1}, c2{0}, cg{0}, c2g{0};

    static const std::vector<std::string> ParameterNames();
    static const Point& SM_Point();

    Point() {}
    Point(double _kl, double _kt, double _c2, double _cg, double _c2g);
    Point(const std::vector<double>& values);
    std::vector<double> ToVector() const;
};

std::ostream& operator<<(std::ostream& s, const Point& p);
std::istream& operator>>(std::istream& s, Point& p);

struct Benchmark {
    Point point;
    int id;
    std::string name;

    Benchmark() {}
    Benchmark(const Point _point, int _id, const std::string& _name = "");
};

class BenchmarkCollection {
public:
    using BenchmarkMap = std::map<int, Benchmark>;
    using const_iterator = BenchmarkMap::const_iterator;

    static const BenchmarkCollection& Benchmarks();
    static BenchmarkCollection CreateScanBenchmarks(const std::vector<std::vector<double>>& values);

    void insert(const Benchmark& benchmark);
    void insert(const Point& point, int id, const std::string& name = "");

    const Benchmark& at(int id) const;
    const Benchmark& at(const std::string& name) const;
    const Benchmark& at_index(size_t index) const;

    const_iterator begin() const;
    const_iterator end() const;
    size_t size() const;

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

    static const GF_Parameterization& Denominator_13TeV();

    GF_Parameterization();
    GF_Parameterization(const Array& _a);

    double Evaluate(const Point& p) const;

    Array a, a_err;
};

std::ostream& operator<<(std::ostream& s, const GF_Parameterization& p);
std::istream& operator>>(std::istream& s, GF_Parameterization& p);

} // namespace NonResHH_EFT
} // namespace analysis
