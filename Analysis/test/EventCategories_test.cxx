/*! Test PhysicalValue class.
This file is part of https://github.com/hh-italian-group/AnalysisTools. */

#include <iostream>
#include "hh-bbtautau/Analysis/include/AnalysisCategories.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"

#define BOOST_TEST_MODULE EventCategories_t
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

using EC = analysis::EventCategory;
using WP = analysis::DiscriminatorWP;

BOOST_AUTO_TEST_CASE(ec_to_string)
{
    {
        EC e(2,2,true, analysis::DiscriminatorWP::Medium, boost::optional<bool>(), true);
        BOOST_TEST(EC::Parse("2j2b_VBF") == e);
    }

    {
        EC e(2,2,false, analysis::DiscriminatorWP::Medium, boost::optional<bool>(), true);
        BOOST_TEST(EC::Parse("2j2b+_VBF") == e);
    }

    {
        EC e(2,2,false, analysis::DiscriminatorWP::Medium, boost::optional<bool>(), false);
        BOOST_TEST(EC::Parse("2j2b+_noVBF") == e);
    }

    {
        EC e(2,2,true, analysis::DiscriminatorWP::Loose, boost::optional<bool>(), true);
        BOOST_TEST(EC::Parse("2j2Lb_VBF") == e);
    }

    {
        EC e(2,2,true, analysis::DiscriminatorWP::Medium, false, true);
        BOOST_TEST(EC::Parse("2j2bR_VBF") == e);
    }

    {
        EC e(2,2,false, analysis::DiscriminatorWP::Medium, false, true);
        BOOST_TEST(EC::Parse("2j2b+R_VBF") == e);
    }

    {
        EC e(2,2,false, analysis::DiscriminatorWP::Medium, false, false);
        BOOST_TEST(EC::Parse("2j2b+R_noVBF") == e);
    }

    {
        EC e(2,2,true, analysis::DiscriminatorWP::Loose, false, true);
        BOOST_TEST(EC::Parse("2j2LbR_VBF") == e);
    }

    {
        EC e(2,2,true, analysis::DiscriminatorWP::Loose, true, true);
        BOOST_TEST(EC::Parse("2j2LbB_VBF") == e);
    }

    {
        EC e(2,2,true, analysis::DiscriminatorWP::Medium, true, true);
        BOOST_TEST(EC::Parse("2j2bB_VBF") == e);
    }

    {
        EC e(2,2,false, analysis::DiscriminatorWP::Medium, true, true);
        BOOST_TEST(EC::Parse("2j2b+B_VBF") == e);
    }

    {
        EC e(2,2,false, analysis::DiscriminatorWP::Medium, true, false);
        BOOST_TEST(EC::Parse("2j2b+B_noVBF") == e);
    }

    {
        EC e(2,2,true, analysis::DiscriminatorWP::Medium, boost::optional<bool>(), true);
        BOOST_TEST("2j2b_VBF" == e.ToString());
    }

    {
        EC e(2,2,false, analysis::DiscriminatorWP::Medium, boost::optional<bool>(), true);
        BOOST_TEST("2j2b+_VBF" ==e.ToString());
    }

    {
        EC e(2,2,false, analysis::DiscriminatorWP::Medium, boost::optional<bool>(), false);
        BOOST_TEST("2j2b+_noVBF" == e.ToString());
    }

    {
        EC e(2,2,true, analysis::DiscriminatorWP::Loose, boost::optional<bool>(), true);
        BOOST_TEST("2j2Lb_VBF" ==e.ToString());
    }

    {
        EC e(2,2,false, analysis::DiscriminatorWP::Loose, boost::optional<bool>(), false);
        BOOST_TEST("2j2Lb+_noVBF" ==e.ToString());
    }

    {
        EC e(2,2,true, analysis::DiscriminatorWP::Medium, false, true);
        BOOST_TEST("2j2bR_VBF" == e.ToString());
    }

    {
        EC e(2,2,false, analysis::DiscriminatorWP::Medium, false, true);
        BOOST_TEST("2j2b+R_VBF" ==e.ToString());
    }

    {
        EC e(2,2,false, analysis::DiscriminatorWP::Medium, false, false);
        BOOST_TEST("2j2b+R_noVBF" == e.ToString());
    }

    {
        EC e(2,2,true, analysis::DiscriminatorWP::Loose, false, true);
        BOOST_TEST("2j2LbR_VBF" == e.ToString());
    }

    {
        EC e(2,2,true, analysis::DiscriminatorWP::Loose, true, true);
        BOOST_TEST("2j2LbB_VBF" == e.ToString());
    }

    {
        EC e(2,2,true, analysis::DiscriminatorWP::Medium, true, true);
        BOOST_TEST("2j2bB_VBF" == e.ToString());
    }

    {
        EC e(2,2,false, analysis::DiscriminatorWP::Medium, true, true);
        BOOST_TEST("2j2b+B_VBF" == e.ToString());
    }

    {
        EC e(2,2,false, analysis::DiscriminatorWP::Medium, true, false);
        BOOST_TEST("2j2b+B_noVBF" == e.ToString());
    }
}

BOOST_AUTO_TEST_CASE(determine_ec_1)
{
    const std::map<WP, size_t> btags = { {WP::Medium, 0} };
    BOOST_TEST(EC::Parse("Inclusive").Contains(1, btags, false, false));
    BOOST_TEST(EC::Parse("0j").Contains(1, btags, false, false));
    BOOST_TEST(EC::Parse("1j").Contains(1, btags, false, false));
    BOOST_TEST(EC::Parse("1j0b").Contains(2, btags, false, false));
    BOOST_TEST(EC::Parse("1j0b+").Contains(1, btags, false, false));
    BOOST_TEST(EC::Parse("2j").Contains(2, btags, false, false));
    // BOOST_CHECK_THROW(!EC::Parse("2j1LbR_VBF").Contains(2, btags, true, false), analysis::exception);
}

BOOST_AUTO_TEST_CASE(determine_ec_2)
{
    const std::map<WP, size_t> btags = { {WP::Loose, 3}, {WP::Medium, 2} };
    BOOST_TEST(EC::Parse("1j1b+").Contains(3, btags, false, false));
    BOOST_TEST(EC::Parse("Inclusive").Contains(3, btags, false, false));
    BOOST_TEST(EC::Parse("0j").Contains(3, btags, false, false));
    BOOST_TEST(EC::Parse("3j3LbR_noVBF").Contains(3, btags, false, false));
    BOOST_TEST(EC::Parse("2j2b+").Contains(3, btags, false, false));
    BOOST_TEST(EC::Parse("2j2b").Contains(3, btags, false, false));
    BOOST_TEST(EC::Parse("2j1b+R").Contains(3, btags, false, false));
    BOOST_TEST(EC::Parse("2j2b+R_noVBF").Contains(3, btags, false, false));
    BOOST_TEST(!EC::Parse("2j0b").Contains(3, btags, false, false));
    BOOST_TEST(EC::Parse("2j2b+_noVBF").Contains(3, btags, false, false));
    BOOST_TEST(!EC::Parse("2j2b_VBF").Contains(3, btags, false, false));
    BOOST_TEST(EC::Parse("2j2b").Contains(3, btags, false, false));
}
