/*! Parameters, benchmarks and kinematic parameterizations used in HH effective-field theory.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include "AnalysisTools/Core/include/TextIO.h"

namespace analysis {

namespace NonResHH_EFT {

struct Point {
    double kl{1}, kt{1}, c2{0}, cg{0}, c2g{0};

    Point() {}
    Point(double _kl, double _kt, double _c2, double _cg, double _c2g) :
        kl(_kl), kt(_kt), c2(_c2), cg(_cg), c2g(_c2g) {}
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

    GF_Parameterization() {}
    GF_Parameterization(const Array& _a) : a(_a) {}

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

    Array a;
};

inline std::ostream& operator<<(std::ostream& s, const GF_Parameterization& p)
{
    s << p.a[0];
    for(size_t n = 1; n < GF_Parameterization::n_params; ++n)
        s << " " << p.a[n];
    return s;
}

inline std::istream& operator>>(std::istream& s, GF_Parameterization& p)
{
    try {
        auto values = ReadValueList(s, GF_Parameterization::n_params, true, " /t", true);
        for(size_t n = 0; n < GF_Parameterization::n_params; ++n)
            p.a[n] = Parse<double>(values.at(n));
    } catch(std::exception& e) {
        throw exception("Invalid EFT parameterization definition. %1%") % e.what();
    }
    return s;
}

} // namespace NonResHH_EFT
} // namespace analysis

