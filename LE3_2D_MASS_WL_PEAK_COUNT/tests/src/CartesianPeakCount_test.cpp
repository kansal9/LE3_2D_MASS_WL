/**
 * @file tests/src/CartesianPeakCount_test.cpp
 * @date 06/30/22
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

#include "LE3_2D_MASS_WL_PEAK_COUNT/CartesianPeakCount.h"
#include "LE3_2D_MASS_WL_UTILITIES/PatchDef.h"

#include <boost/test/unit_test.hpp>
#include <random>

using LE3_2D_MASS_WL_PEAK_COUNT::CartesianPeakCount;
using LE3_2D_MASS_WL_UTILITIES::PatchDef;
using LE3_2D_MASS_WL_PEAK_COUNT::estimator;

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE (CartesianPeakCount_test)
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( cartesianPeakCount_simple_test ) {

    PatchDef patch(0, 0, 10, 1);
    CartesianPeakCount pc(patch, 0, 10, estimator::APERTURE_MASS);

    Matrix map(10, 10);
    map.clear();
    map(1,1) = 10;
    map(5,5) = 11;

    pc.findPeaksAtTheta(map, 0);

    auto peakCat = pc.getPeakCatalog();

    std::cout << "Number of peak: " << peakCat.getNentries() << std::endl;

    BOOST_CHECK_EQUAL(peakCat.getNentries(), 2);

    double ra, dec, zmin, zmax, theta, snr;
    for(unsigned int i=0; i<peakCat.getNentries(); i++)
    {
        peakCat.getEntry(i, ra, dec, zmin, zmax, theta, snr);
        BOOST_CHECK_EQUAL(snr, 10+i);
        BOOST_CHECK_EQUAL(zmin, 0);
        BOOST_CHECK_EQUAL(zmax, 10);
        BOOST_CHECK_EQUAL(theta, 0);
        printf("%d %f %f %f %f %f %f\n", i, ra, dec, zmin, zmax, theta, snr);
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( cartesianPeakCount_wavelet_test ) {
    PatchDef patch(0, 0, 10, 1);
    CartesianPeakCount pc(patch, 0, 10, estimator::WAVELET);
    ConvergenceMap convMap(10, 10, 1);
    convMap[0](5, 5) = 10;
    pc.findPeaksFromConvergence(convMap, 2);
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
BOOST_AUTO_TEST_CASE( cartesianPeakCount_aperture_test ) {

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> normal_dist{0,1};

    PatchDef patch(0, 0, 10, 1);
    CartesianPeakCount pc(patch, 0, 10, estimator::APERTURE_MASS);

    ShearMap shearMap(10, 10, 2);
    ShearMap noisyShearMap(10, 10, 2);

    // fill maps with random values from N(0,1)
    for(size_t i=0; i<shearMap.getXdim(); i++)
    {
        for(size_t j=0; j<shearMap.getYdim(); j++)
        {
            shearMap[0](i,j) = normal_dist(gen);
            shearMap[1](i,j) = normal_dist(gen);
            noisyShearMap[0](i,j) = normal_dist(gen);
            noisyShearMap[1](i,j) = normal_dist(gen);
        }
    }

    shearMap[0](5, 5) = 10;

    std::vector<double> radius = {1, 2};
    pc.findPeaksFromShear(shearMap, noisyShearMap, radius);
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


