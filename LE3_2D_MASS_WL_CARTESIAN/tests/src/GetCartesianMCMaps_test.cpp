/**
 * @file tests/src/GetCartesianMCMaps_test.cpp
 * @date 10/13/20
 * @author vkansal
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

#include "ElementsKernel/Temporary.h"
#include "ElementsServices/DataSync.h"

#include "LE3_2D_MASS_WL_UTILITIES/PatchDef.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GetCartesianMCMaps.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include <ios>
#include <sstream>
#include <iostream>

using namespace ElementsServices::DataSync;
using namespace LE3_2D_MASS_WL_CARTESIAN;
using namespace LE3_2D_MASS_WL_UTILITIES;
namespace fs = boost::filesystem;

// handle on created path names
fs::path test_path;

struct GetCartesianMCMapsTestEnv
{
    DataSync sync;
    path testParamFile;
    std::vector<fs::path> MCShearMaps;
    CartesianParam params;
    CatalogData cat;
    double rmin, rmax, dmin, dmax, zmin, zmax;

    GetCartesianMCMapsTestEnv() :
            sync("LE3_2D_MASS_WL_CARTESIAN/datasync_webdav.conf",
                 "LE3_2D_MASS_WL_CARTESIAN/test_file_list.txt"),
            testParamFile(sync.absolutePath("Cparam_test_GetCartesianMCMaps.xml"))
    {
        sync.download();
        CatalogData dummy;
        if (checkFileType(testParamFile, signXML))
        {
            if (fileHasField(testParamFile,"DpdTwoDMassParamsConvergencePatch"))
            {
                params.readConvPatchXMLFile(testParamFile.native(), dummy);
            }
        }
        auto& patch = params.getPatches()[0];
        rmin = patch.getRaMin();
        rmax = patch.getRaMax();
        dmin = patch.getDecMin();
        dmax = patch.getDecMax();
        cat.fillTest(100, "LENSMC", rmin*rad2deg, rmax*rad2deg,
                     dmin*rad2deg, dmax*rad2deg);
        vecMinMax(cat["z"], zmin, zmax);
    }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (GetCartesianMCMaps_test, GetCartesianMCMapsTestEnv)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( getDeNoisedShearMap_test )
{
    std::cout << "-- GetCartesianMCMaps: getDeNoisedShearMap_test" << std::endl;
    // PatchDef CB(rmin, rmax, dmin, dmax, zmin, zmax);
    GetCartesianMCMaps mc(cat, params);
    ShearMap DeNoisedShearMap;
    mc.getDeNoisedShearMap(DeNoisedShearMap);
    // TODO : implement a relevant test or delete
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getNoisedShearMaps_test )
{
    std::cout << "-- GetCartesianMCMaps: getNoisedShearMaps_test" << std::endl;
    // PatchDef CB(rmin, params.getRaMax(rmin), dmin, params.getDecMax(dmin), zmin, zmax);
    GetCartesianMCMaps mc(cat, params);
    std::vector<ShearMap> NoisedShearMapList;
    mc.getNoisedShearMaps(NoisedShearMapList);
    BOOST_CHECK(NoisedShearMapList.size() > 0);
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( performAddition_test )
{
    std::cout << "-- GetCartesianMCMaps: performAddition_test" << std::endl;
    // PatchDef CB(rmin, params.getRaMax(rmin), dmin, params.getDecMax(dmin), zmin, zmax);
    GetCartesianMCMaps mc(cat, params);
    ShearMap DeNoisedShearMap;
    mc.getDeNoisedShearMap(DeNoisedShearMap);
    std::vector<ShearMap> NoisedShearMapList;
    mc.getNoisedShearMaps(NoisedShearMapList);
    Elements::TempDir one;
    test_path = one.path();
    for (size_t i = 0; i < NoisedShearMapList.size(); i++)
    {
        fs::path shearMapName(
                "EUC_LE3_WL_ShearMap_0" + std::to_string(i) + ".fits");
        MCShearMaps.push_back(test_path / shearMapName);
    }

    for (size_t i = 0; i < NoisedShearMapList.size(); i++)
    {
        mc.performAddition(DeNoisedShearMap, NoisedShearMapList[i],
                           MCShearMaps[i].native());
    }
    for (size_t i = 0; i < MCShearMaps.size(); i++)
    {
        BOOST_CHECK(boost::filesystem::exists(MCShearMaps[i]));
    }
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
