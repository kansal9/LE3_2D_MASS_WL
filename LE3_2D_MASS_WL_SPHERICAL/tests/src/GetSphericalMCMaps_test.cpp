/**
 * @file tests/src/GetSphericalMCMaps_test.cpp
 * @date 10/13/20
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

#include "ElementsServices/DataSync.h"
#include "ElementsKernel/Temporary.h"
#include "ElementsKernel/Logging.h"

#include "LE3_2D_MASS_WL_SPHERICAL/GetSphericalMCMaps.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalUtils.h"

#include <ios>
#include <sstream>
#include <iostream>
#include "../../LE3_2D_MASS_WL_SPHERICAL/SphMapMaker.h"
#include "../../LE3_2D_MASS_WL_SPHERICAL/SphMassMapping.h"

using LE3_2D_MASS_WL_SPHERICAL::GetSphericalMCMaps;
using LE3_2D_MASS_WL_SPHERICAL::SphericalParam;

using namespace Euclid::WeakLensing::TwoDMass::Spherical;
using namespace Euclid::WeakLensing::TwoDMass;
using namespace ElementsServices::DataSync;

// handle on created path names
boost::filesystem::path test_path;
using boost::filesystem::path;
static Elements::Logging logger = Elements::Logging::getLogger(
        "SphericalMCMaps_Test");

struct GetSphericalMCMapsTestEnv
{
    DataSync sync;
    path testParamFile;
    std::vector<fs::path> MCShearMaps;
    SphericalParam params;
    CatalogData cat;
    GetSphericalMCMapsTestEnv() :
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
BOOST_FIXTURE_TEST_SUITE (GetSphericalMCMaps_test, GetSphericalMCMapsTestEnv)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( getDeNoisedShearMap_test )
{
    logger.info() << "-- GetSphericalMCMaps: getDeNoisedShearMap_test";
    GetSphericalMCMaps smc(cat, params);
    std::pair<Healpix_Map<double>, Healpix_Map<double> > DenoisedShearMapPair =
            smc.getDeNoisedShearMap();
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( getNoisedShearMaps_test )
{
    logger.info() << "-- GetSphericalMCMaps: getNoisedShearMaps_test";
    GetSphericalMCMaps smc(cat, params);
    std::vector<std::pair<Healpix_Map<double>, Healpix_Map<double> > > ShearMapList;
    ShearMapList = smc.getNoisedShearMaps();
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( performAddition_test )
{
    logger.info() << "-- GetSphericalMCMaps: performAddition_test";
    GetSphericalMCMaps smc(cat, params);

    std::pair<Healpix_Map<double>, Healpix_Map<double> > DenoisedShearMapPair =
            smc.getDeNoisedShearMap();

    std::vector<std::pair<Healpix_Map<double>, Healpix_Map<double> > > NoisedShearMapList;
    NoisedShearMapList = smc.getNoisedShearMaps();

    Elements::TempDir one;
    test_path = one.path();

    for (size_t i = 0; i < NoisedShearMapList.size(); i++)
    {
        boost::filesystem::path shearMapName = (fs::path(
                "EUC_LE3_WL_ShearMap_0" + std::to_string(i) + ".fits"));
        MCShearMaps.push_back(test_path / shearMapName);
    }

    for (size_t i = 0; i < NoisedShearMapList.size(); i++)
    {
        smc.performAddition(DenoisedShearMapPair, NoisedShearMapList[i],
                MCShearMaps[i].native());
    }

    for (size_t i = 0; i < MCShearMaps.size(); i++)
    {
        BOOST_CHECK(boost::filesystem::exists(MCShearMaps[i]));
    }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
