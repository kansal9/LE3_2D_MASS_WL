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
//#define private public
//#define protected public
#include "ElementsKernel/Logging.h"
#include <boost/filesystem.hpp>
#include "ElementsKernel/Auxiliary.h"
#include "ElementsServices/DataSync.h"

#include <boost/test/unit_test.hpp>
#include "LE3_2D_MASS_WL_CARTESIAN/MassMapping.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MapMaker.h"
#include "LE3_2D_MASS_WL_CARTESIAN/InpaintingAlgo.h"

#include <ios>
#include <sstream>
#include <iostream>

using namespace LE3_2D_MASS_WL_CARTESIAN;
using namespace Euclid::WeakLensing::TwoDMass;
using LE3_2D_MASS_WL_CARTESIAN::CartesianParam;
using namespace ElementsServices::DataSync;
using std::string;

using LE3_2D_MASS_WL_CARTESIAN::InpaintingAlgo;
static Elements::Logging logger = Elements::Logging::getLogger("Inpainting_test");

//-----------------------------------------------------------------------------

struct InpaintingAlgoDataSyncFixture {
  DataSync sync; // This is the synchronizer
  // These are just shortcuts
  path testParamFile, testParamFile2, testshearMap, testconvMap;
  InpaintingAlgoDataSyncFixture () :
      sync(
          // Here is the connection configuration file
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

//BOOST_AUTO_TEST_SUITE (InpaintingAlgo_test)
BOOST_FIXTURE_TEST_SUITE (InpaintingAlgo_test, InpaintingAlgoDataSyncFixture)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( copyConstructor_test ) {
 logger.info() << "copyConstructor_test";
  CartesianParam params;
   // Case 1: check file is of XML format
   if (true == checkFileType(testParamFile, Euclid::WeakLensing::TwoDMass::signXML)) {
    // read XML, get parameter file type
    // check if it's to get a patch
    if (true == fileHasField(testParamFile, "DpdTwoDMassParamsConvergencePatch")) {
      params.ReadConvPatchXMLFile (testParamFile.native());
    }
   }

  // Import shear and convergence maps
  ShearMap myShearMap(testshearMap.native());
  ConvergenceMap myConvMap(testconvMap.native());

  // Create the InPaintingAlgo object with the maps
  InpaintingAlgo *myInPainting = new InpaintingAlgo(myShearMap, myConvMap, params);

  // Create an InPaintingAlgo using copy constructor
  InpaintingAlgo *myCopy = new InpaintingAlgo(*myInPainting);

  // Delete the original object
  delete myInPainting;
  myInPainting = nullptr;

  // Then try to perform inpainting (should not work if copy constructor did not work)
  ConvergenceMap* myInPaintedMap = myCopy->performInPaintingAlgo();
  // Check a convergence map is well returned
  BOOST_CHECK(myInPaintedMap!=nullptr);

  delete myInPaintedMap;
  myInPaintedMap = nullptr;

  delete myCopy;
  myCopy = nullptr;
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( Inpaint_test ) {
 logger.info() << "Inpaint_test";
  CartesianParam params;
   // Case 1: check file is of XML format
   if (true == checkFileType(testParamFile, Euclid::WeakLensing::TwoDMass::signXML)) {
    // read XML, get parameter file type
    // check if it's to get a patch
    if (true == fileHasField(testParamFile, "DpdTwoDMassParamsConvergencePatch")) {
      params.ReadConvPatchXMLFile (testParamFile.native());
    }
   }

  // Import shear and convergence maps
  ShearMap myShearMap(testshearMap.native());
  ConvergenceMap myConvMap(testconvMap.native());

  InpaintingAlgo myInPainting(myShearMap, myConvMap, params);

  // Then try to perform inpainting
  LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap *myInPaintedMap = myInPainting.performInPaintingAlgo();
  // Check a convergence map is well returned
  BOOST_CHECK(myInPaintedMap!=nullptr);
  //BOOST_CHECK(true);

  delete myInPaintedMap;
  myInPaintedMap = nullptr;

}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( Inpaint_test2 ) {
 logger.info() << "Inpaint_test apply wavelet boundary";
  CartesianParam params;
   // Case 1: check file is of XML format
   if (true == checkFileType(testParamFile, Euclid::WeakLensing::TwoDMass::signXML)) {
    // read XML, get parameter file type
    // check if it's to get a patch
    if (true == fileHasField(testParamFile2, "DpdTwoDMassParamsConvergencePatch")) {
      params.ReadConvPatchXMLFile (testParamFile2.native());
    }
   }

  // Import shear and convergence maps
  ShearMap myShearMap(testshearMap.native());
  ConvergenceMap myConvMap(testconvMap.native());

  InpaintingAlgo myInPainting(myShearMap, myConvMap, params);

  // Then try to perform inpainting
  LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap *myInPaintedMap = myInPainting.performInPaintingAlgo();
  // Check a convergence map is well returned
  BOOST_CHECK(myInPaintedMap!=nullptr);
  //BOOST_CHECK(true);

  delete myInPaintedMap;
  myInPaintedMap = nullptr;

}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( Inpaint_Block_test ) {
 logger.info() << "Inpaint_Block_test ";
  CartesianParam params;
   // Case 1: check file is of XML format
   if (true == checkFileType(testParamFile, Euclid::WeakLensing::TwoDMass::signXML)) {
    // read XML, get parameter file type
    // check if it's to get a patch
    if (true == fileHasField(testParamFile, "DpdTwoDMassParamsConvergencePatch")) {
      params.ReadConvPatchXMLFile (testParamFile.native());
    }
   }

  // Import shear and convergence maps
  ShearMap myShearMap(testshearMap.native());
  ConvergenceMap myConvMap(testconvMap.native());

  InpaintingAlgo myInPainting(myShearMap, myConvMap, params);

  // Then try to perform inpainting
  LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap *myInPaintedMap = myInPainting.performInPaintingAlgo(64, 64);
  // Check a convergence map is well returned
  BOOST_CHECK(myInPaintedMap!=nullptr);
  //BOOST_CHECK(true);

  delete myInPaintedMap;
  myInPaintedMap = nullptr;

}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
