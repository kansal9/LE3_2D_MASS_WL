/**
 * @file tests/src/SphericalUtils_test.cpp
 * @date 06/17/21
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

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "ElementsKernel/Logging.h"
#include "ElementsKernel/Temporary.h"
#include "ElementsServices/DataSync.h"

#include "LE3_2D_MASS_WL_SPHERICAL/SphericalUtils.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"

#include <iostream>
#include "../../LE3_2D_MASS_WL_SPHERICAL/SphMapMaker.h"

using LE3_2D_MASS_WL_SPHERICAL::SphMapMaker;
using LE3_2D_MASS_WL_SPHERICAL::SphericalParam;
using Euclid::WeakLensing::TwoDMass::Spherical::SphericalIO;

using namespace Euclid::WeakLensing::TwoDMass::Spherical;
using namespace Euclid::WeakLensing::TwoDMass;
using namespace ElementsServices::DataSync;

using boost::filesystem::path;
boost::filesystem::path test_path;

static Elements::Logging logger = Elements::Logging::getLogger(
        "SphericalUtils_test");

//-----------------------------------------------------------------------------
struct SphericalUtilsDataSyncFixture
{
    DataSync sync;
    path testParamFile, fitsFilePath;
    SphericalParam params;
    SphericalUtilsDataSyncFixture() :
            sync("LE3_2D_MASS_WL_SPHERICAL/datasync_webdav.conf",
                 "LE3_2D_MASS_WL_SPHERICAL/test_file_list.txt"),
            testParamFile(sync.absolutePath("Sparam_test.xml")),
            fitsFilePath(sync.absolutePath("data/Gamma_NSide16_Sphere.fits"))
    {
        sync.download();
        CatalogData dummy;
        readSphericalParameterFile(testParamFile, params, dummy);
    }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (SphericalUtils_test, SphericalUtilsDataSyncFixture)
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( readSphericalParameterFile_test )
{
    std::cout << "-- SphericalUtils: readSphericalParameterFile_Test"
              << std::endl;
    int nside = params.getNside();
    BOOST_CHECK(nside == 16);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( applyGaussianFilter_hp_test )
{
    std::cout << "-- SphericalUtils: applyGaussianFilter_hp_Test" << std::endl;
    double sigma = 0.003;
    Healpix_Map<double> map;
    read_Healpix_map_from_fits((fitsFilePath).native(), map);
    applyGaussianFilter_hp(map, sigma);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( spline_test )
{
    logger.info() << "-- SphericalUtils: spline_test";
    Healpix_Map<double> map(7, RING); //nside = 128
    map.fill(1.);
    // Perform an arbitrary wavelet decomposition
    std::vector<Healpix_Map<double> > myBand = transformBspline_hp(map, 5);

    // Reconstruct the image
    Healpix_Map<double> mapBack = reconsBspline_hp(myBand);

    int m_npix = map.Npix();
    // Check the values of the original and back image are close
    for (int i = 0; i < m_npix; i++)
    {
        BOOST_CHECK_CLOSE(map[i], mapBack[i], 0.1);
    }

}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( applyThreshold_test )
{
    logger.info() << "-- SphericalUtils: applyThreshold_test";

    Alm<xcomplex<double> > map_lm( 10, 10);
    map_lm.SetToZero();
    for (int l = 0; l <= 10; ++l)
    {
        for (int m = 0; m <= l; ++m)
        {
            map_lm(l, m) = double(l + m);
            ;
        }
    }
    applyThreshold(map_lm, 10, 10);

    for (int l = 0; l <= 10; ++l)
    {
        for (int m = 0; m <= l; ++m)
        {
            //  if ((l+m) < 10) {
            //   BOOST_CHECK_CLOSE(map_lm(l,m), 0, 1);
            //  } else {
            //   BOOST_CHECK_CLOSE(map_lm(l,m), l+m, 1);
            //  }
        }
    }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( max_abs_alm_test )
{
    logger.info() << "-- SphericalUtils: max_abs_alm_test";

    Alm<xcomplex<double> > map_lm( 10, 10);
    map_lm.SetToZero();
    for (int l = 0; l <= 10; ++l)
    {
        for (int m = 0; m <= l; ++m)
        {
            map_lm(l, m) = double(l + m);
            ;
        }
    }
    double max = max_abs_alm(map_lm, 10);

    BOOST_CHECK_CLOSE(max, 20, 1);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( getStd_test )
{
    logger.info() << "-- SphericalUtils: getStd_test";

    Healpix_Map<double> map(7, RING); //nside = 128
    map.fill(1.);
    int m_npix = map.Npix();

    for (int i = 0; i < m_npix; i++)
    {
        map[i] = double(i + i);
    }
    double stdev = getStd(map);
    logger.info() << "stdev: " << stdev;
    BOOST_CHECK_CLOSE(stdev, 113511.68, 1);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
