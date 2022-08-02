/**
 * @file tests/src/MapMaker_test.cpp
 * @date 05/13/19
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

#include "ElementsKernel/Logging.h"
#include "ElementsServices/DataSync.h"

#include "LE3_2D_MASS_WL_CARTESIAN/MapMaker.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GenericMap.h"
#include "LE3_2D_MASS_WL_UTILITIES/PatchDef.h"
#include <ios>
#include <sstream>
#include <iostream>

using namespace LE3_2D_MASS_WL_CARTESIAN;
using namespace LE3_2D_MASS_WL_UTILITIES;
using namespace ElementsServices::DataSync;

using LE3_2D_MASS_WL_CARTESIAN::CartesianParam;

using boost::filesystem::path;
static Elements::Logging logger = Elements::Logging::getLogger("MapMaker_test");

//-----------------------------------------------------------------------------
struct MapMakerTestEnv
{
    DataSync sync;
    path testParamFile;
    CartesianParam params;
    ShearMap m_ShearMap;
    CatalogData m_cat;
    double rmin, rmax, dmin, dmax;
    MapMakerTestEnv() :
        sync("LE3_2D_MASS_WL_CARTESIAN/datasync_webdav.conf",
             "LE3_2D_MASS_WL_CARTESIAN/test_file_list.txt"),
        testParamFile(sync.absolutePath("Cparam_test.xml"))

    {
        sync.download();

        if (fileHasField(testParamFile, "DpdTwoDMassParamsConvergencePatch"))
        {
            logger.info() << "Parameter file is for Convergence Patch..";
            params.readConvPatchXMLFile(testParamFile.native(), m_cat);
        }

        auto& patch = params.getPatches()[0];
        m_ShearMap = ShearMap(patch.getXbin(), patch.getYbin(), 3);
        rmin = patch.getRaMin();
        rmax = patch.getRaMax();
        dmin = patch.getDecMin();
        dmax = patch.getDecMax();
        m_cat.fillTest(100, "LENSMC", rmin*rad2deg, rmax*rad2deg,
                                      dmin*rad2deg, dmax*rad2deg);
    }
};
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE (MapMaker_test)
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE(get_shearMap_test, MapMakerTestEnv)
{
    logger.info() << "get_shearMap_test";
    MapMaker map(m_cat);
    std::cout << "rmin: " << rmin << std::endl;
    std::cout << "rmax: " << rmax << std::endl;
    std::cout << "dmin: " << rmin << std::endl;
    std::cout << "dmax: " << dmax << std::endl;
    auto& patch = params.getPatches()[0];
    map.getShearMap(patch, m_ShearMap, 0, 10);
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
