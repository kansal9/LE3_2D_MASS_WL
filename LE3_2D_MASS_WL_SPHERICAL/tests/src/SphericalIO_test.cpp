/**
 * @file tests/src/SphericalIO_test.cpp
 * @date 12/24/20
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

#include "ElementsServices/DataSync.h"
#include "ElementsKernel/Temporary.h"

#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalUtils.h"

#include <iostream>
#include "../../LE3_2D_MASS_WL_SPHERICAL/SphMapMaker.h"

using LE3_2D_MASS_WL_SPHERICAL::SphMapMaker;
using LE3_2D_MASS_WL_SPHERICAL::SphericalParam;
using Euclid::WeakLensing::TwoDMass::Spherical::SphericalIO;

using namespace Euclid::WeakLensing::TwoDMass::Spherical;
using namespace Euclid::WeakLensing::TwoDMass;
using namespace ElementsServices::DataSync;

boost::filesystem::path test_path;
//-----------------------------------------------------------------------------
struct SphericalIOFixture
{
    DataSync sync;
    path testParamFile;
    SphericalParam params;
    CatalogData cat;
    SphericalIOFixture() :
        sync("LE3_2D_MASS_WL_SPHERICAL/datasync_webdav.conf",
             "LE3_2D_MASS_WL_SPHERICAL/test_file_list.txt"),
        testParamFile(sync.absolutePath("Sparam_test.xml"))

    {
        sync.download();
        cat.fillTest(100);
        readSphericalParameterFile(testParamFile, params, cat);
    }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (SphericalIO_test, SphericalIOFixture)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( wrtieMap_test )
{
    std::cout << "-- SphericalIO: writeMap_Test" << std::endl;

    using Elements::TempDir;
    TempDir one;
    test_path = one.path();

    SphMapMaker mapMaker(params);
    Healpix_Map<double> g1_hmap;
    Healpix_Map<double> g2_hmap;
    Healpix_Map<double> ngal_hp;
    mapMaker.create_ShearMap(cat, g1_hmap, g2_hmap, ngal_hp);

    SphericalIO SphericalIO(params);
    boost::filesystem::path filename = test_path / "GammaMap.fits";
    std::cout << "-- SphericalIO: writing Map1" << std::endl;
    SphericalIO.write_Map(filename.native(), g1_hmap, "GAMMA1");
    std::cout << "-- SphericalIO: writing Map2" << std::endl;
    SphericalIO.write_Map(filename.native(), g2_hmap, "GAMMA2");
    std::cout << "-- SphericalIO: writing GalCount Map" << std::endl;
    boost::filesystem::path galfilename = test_path / "GalCountMap.fits";
    SphericalIO.write_Map(galfilename.native(), ngal_hp, "GALCOUNT");
    std::cout << "-- SphericalIO: Complete test" << std::endl;
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( writeSphericalXMLfile_test )
{
    std::cout << "-- SphericalIO: writeSphericalXMLfile _Test" << std::endl;

    using Elements::TempDir;
    TempDir one;
    test_path = one.path();

    SphericalIO SphericalIO(params);

    boost::filesystem::path outputKappa = boost::filesystem::path(
            "SphereConvergence.json");
    std::ofstream outfile;
    outfile.open((test_path / outputKappa).string(), std::ios_base::app);
    outfile
            << "[NoisySphereConvergence.fits, DenoisedSphereConvergence.fits, MCSphereConvergence1.fits, MCSphereConvergence2.fits]";
    outfile.close();

    boost::filesystem::path outputGalCount = boost::filesystem::path(
            "SphereGalCount.json");
    std::ofstream galOutfile;
    galOutfile.open((test_path / outputGalCount).string(), std::ios_base::app);
    galOutfile << "[abc.fits]";
    galOutfile.close();
    boost::filesystem::path outXMLConvergenceMap = boost::filesystem::path(
            "outXMLConvergenceMap.xml");
    std::cout << "writeSphericalXMLfile" << std::endl;
    SphericalIO.writeSphericalXMLfile(test_path, outputKappa, outputGalCount,
            outXMLConvergenceMap);
    std::cout << "writeSphericalMCXMLfile" << std::endl;
    SphericalIO.writeSphericalMCXMLfile(test_path, outputKappa,
            outXMLConvergenceMap);
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
