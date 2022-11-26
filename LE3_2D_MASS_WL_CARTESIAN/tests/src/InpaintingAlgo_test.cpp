/**
 * @file tests/src/InpaintingAlgo_test.cpp
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

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MapMaker.h"
#include "LE3_2D_MASS_WL_CARTESIAN/InpaintingAlgo.h"

#include <ios>
#include <sstream>
#include <iostream>

using namespace LE3_2D_MASS_WL_CARTESIAN;
using namespace ElementsServices::DataSync;
using namespace LE3_2D_MASS_WL_UTILITIES;
using std::string;

static Elements::Logging logger = Elements::Logging::getLogger(
        "Inpainting_test");

//-----------------------------------------------------------------------------

struct InpaintingAlgoDataSyncFixture
{
    DataSync sync; // This is the synchronizer
    // These are just shortcuts
    path testParamFile, testParamFile2, testshearMap, testconvMap;
    InpaintingAlgoDataSyncFixture() :
            sync(// Here is the connection configuration file
                 "LE3_2D_MASS_WL_CARTESIAN/datasync_webdav.conf",
                 // Here is the dependency configuration file
                 "LE3_2D_MASS_WL_CARTESIAN/test_file_list.txt"),
              testParamFile(sync.absolutePath("Cparam_test.xml")),
              testParamFile2(sync.absolutePath("Cparam_Inp_test.xml")),
              testshearMap(sync.absolutePath("data/shearMap_test.fits")),
              testconvMap(sync.absolutePath("data/convMap_test.fits"))
    {
        sync.download();
    }
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (InpaintingAlgo_test, InpaintingAlgoDataSyncFixture)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( inpaint_test )
{
    logger.info() << "Inpaint_test";
    CartesianParam params;
    CatalogData dummy;
    // Case 1: check file is of XML format
    if (checkFileType(testParamFile, signXML))
    {
        // read XML, get parameter file type
        // check if it's to get a patch
        if (fileHasField(testParamFile, "DpdTwoDMassParamsConvergencePatch"))
        {
            params.readConvPatchXMLFile(testParamFile.native(), dummy);
        }
    }

    // Import shear and convergence maps
    ShearMap myShearMap(testshearMap.native());
    ConvergenceMap myConvMap(testconvMap.native());
    GenericMap myMask(myShearMap.getXdim(), myShearMap.getYdim(), 1);

    // Then try to perform inpainting
    InpaintingAlgo myInPainting(myShearMap, myMask, params);
    myInPainting.performInPaintingAlgo(myConvMap);

    // TODO: implement a test

}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( inpaint_test2 )
{
    logger.info() << "Inpaint_test apply wavelet boundary";
    CartesianParam params;
    CatalogData dummy;
    // Case 1: check file is of XML format
    if (checkFileType(testParamFile, signXML))
    {
        // read XML, get parameter file type
        // check if it's to get a patch
        if (fileHasField(testParamFile2, "DpdTwoDMassParamsConvergencePatch"))
        {
            params.readConvPatchXMLFile(testParamFile2.native(), dummy);
        }
    }

    // Import shear and convergence maps
    ShearMap myShearMap(testshearMap.native());
    ConvergenceMap myConvMap(testconvMap.native());
    GenericMap myMask(myShearMap.getXdim(), myShearMap.getYdim(), 1);

    // Then try to perform inpainting
    InpaintingAlgo myInPainting(myShearMap, myMask, params);
    myInPainting.performInPaintingAlgo(myConvMap);

    // TODO: implement a test

}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
