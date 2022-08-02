/**
 * @file tests/src/SphericalPeakCount_test.cpp
 * @date 06/18/21
 * @author user
 *
 * @copyright (C) 2012-2020 Euclid Science Ground Segment
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 3.0 of the License, or (at your option)
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */

#include "LE3_2D_MASS_WL_PEAK_COUNT/SphericalPeakCount.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <random>

using LE3_2D_MASS_WL_PEAK_COUNT::SphericalPeakCount;
using namespace LE3_2D_MASS_WL_UTILITIES;

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE (SphericalPeakCount_test)
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( cartesianPeakCount_simple_test ) {
    Healpix_Map<double> map(4, RING, SET_NSIDE);
    map.fill(0.);
    map[0] = 10;

    unsigned int scale = 0;
    SphericalPeakCount pc(map, scale);
    pc.findPeaksAtTheta(map, scale);

    auto peakCat = pc.getPeakCatalog();

    std::cout << "Number of peak: " << peakCat.getNentries() << std::endl;
    BOOST_CHECK_EQUAL(peakCat.getNentries(), 1);

    double ra, dec, zmin, zmax, theta, snr;
    peakCat.getEntry(0, ra, dec, zmin, zmax, theta, snr);
    BOOST_CHECK_EQUAL(snr, 10);
    BOOST_CHECK_EQUAL(zmin, 0);
    BOOST_CHECK_EQUAL(zmax, 0);
    BOOST_CHECK_EQUAL(theta, 0);
    printf("%f %f %f %f %f %f\n", ra, dec, zmin, zmax, theta, snr);
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( cartesianPeakCount_kappa_test ) {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> normal_dist{0,1};

    Healpix_Map<double> map(4, RING, SET_NSIDE);

    // fill map with random values from N(0,1)
    for(int i=0; i<map.Npix(); i++)
    {
        map[i] = normal_dist(gen);
    }
    auto radec = getIndex2RaDec(map, 0);
    printf("Pix0 is ra/dec: %f %f\n", radec.first, radec.second);
    map[0] = 42;

    unsigned int scale = 2;
    SphericalPeakCount pc(map, scale);
    pc.findPeaks_hp();

    auto peakCat = pc.getPeakCatalog();

    std::cout << "Number of peak: " << peakCat.getNentries() << std::endl;

    double ra, dec, zmin, zmax, theta, snr;
    for(unsigned int i=0; i<peakCat.getNentries(); i++)
    {
        peakCat.getEntry(i, ra, dec, zmin, zmax, theta, snr);
        printf("%d %f %f %f %f %f %f\n", i, ra, dec, zmin, zmax, theta, snr);
    }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
