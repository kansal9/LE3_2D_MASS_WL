/*
 * Copyright (C) 2012-2020 Euclid Science Ground Segment
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
 */

/**
 * @file tests/src/CartesianAlgoKS_test.cpp
 * @date 10/21/19
 * @author user
 */
//#define private public
//#define protected public
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "ElementsKernel/Auxiliary.h"
#include "ElementsKernel/Temporary.h"
#include "ElementsServices/DataSync.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CoordinateBound.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include <ios>
#include <sstream>
#include <iostream>
#include <exception>

using namespace LE3_2D_MASS_WL_CARTESIAN;
using namespace Euclid::WeakLensing::TwoDMass;
using LE3_2D_MASS_WL_CARTESIAN::CartesianParam;
using namespace ElementsServices::DataSync;
using std::string;
namespace fs = boost::filesystem;
 // handle on created path names
 boost::filesystem::path test_path;

//-----------------------------------------------------------------------------

struct CartesianAlgoKSDataSyncFixture {
  DataSync sync; // This is the synchronizer
  // These are just shortcuts
  path testParamFile, testCatFile;
  LE3_2D_MASS_WL_CARTESIAN::ShearMap *m_ShearMap = nullptr;
  LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap *m_ConvergenceMap = nullptr;
  CartesianAlgoKSDataSyncFixture () :
      sync(
          // Here is the connection configuration file
          "LE3_2D_MASS_WL_CARTESIAN/datasync_webdav.conf",
           // Here is the dependency configuration file
          "LE3_2D_MASS_WL_CARTESIAN/test_file_list.txt"),
         testParamFile(sync.absolutePath("Cparam_test.xml")),
         testCatFile(sync.absolutePath("data/lensmcCat.fits"))
 {
    sync.download();
  }
};

//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_SUITE (CartesianAlgoKS_test)
BOOST_FIXTURE_TEST_SUITE (CartesianAlgoKS_test, CartesianAlgoKSDataSyncFixture)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( OverallCartesianAlgo_test ) {
 std::cout <<"OverallCartesianAlgo_test" << std::endl;
 CatalogData read;
 std::vector<std::vector<double> > testCatData;
 read.readCatalog(testCatFile.native(), testCatData);

 double zMin, zMax;
 vecMinMax(testCatData[5], &zMin, &zMax);

 Elements::TempDir one;
 test_path = one.path();

 auto datadir = test_path / "data";
 fs::create_directories (datadir);

  CartesianParam params;
   // Case 1: check file is of XML format
   if (true == Euclid::WeakLensing::TwoDMass::checkFileType(testParamFile, Euclid::WeakLensing::TwoDMass::signXML)) {
    // read XML, get parameter file type
    // check if it's to get a patch
    if (true == Euclid::WeakLensing::TwoDMass::fileHasField(testParamFile, "DpdTwoDMassParamsConvergencePatch")) {
      params.ReadConvPatchXMLFile (testParamFile.native());
    }
   }
 std::vector<double> centerX = params.getMapCenterX();
 std::vector<double> centerY = params.getMapCenterY();
 for (size_t it=0; it<centerX.size(); it++){
  Euclid::WeakLensing::TwoDMass::CartesianKS::CartesianAlgoKS CartesainAlgo(params);
  double ramin, decmin;
  ramin = params.getRaMin(centerX[it]);
  decmin = params.getDecMin(centerY[it]);
  LE3_2D_MASS_WL_CARTESIAN::CoordinateBound m_CB(ramin, params.getRaMax(ramin), decmin,
                                                 params.getDecMax(decmin), zMin, zMax);
  bool status;
  std::cout <<"ExtractShearMap_test" << std::endl;
  boost::filesystem::path shearMap= "ExtractedShearMapFile.fits";
  status = CartesainAlgo.extractShearMap ((datadir/shearMap).native(), testCatData, m_CB);
  if (boost::filesystem::exists( (datadir/shearMap).native() )) {
    BOOST_CHECK(true);
  }
  BOOST_CHECK(status);

  std::cout <<"PerformKSMassMapping_test" << std::endl;
  boost::filesystem::path outConvergenceMap= "ConvMapKSFile.fits";
  boost::filesystem::path outConvergenceMapXML= "ConvMapKSFile.xml";
  /*status = CartesainAlgo.performKSMassMapping((test_path/shearMap).native(), (test_path/outConvergenceMap).native());
  if (boost::filesystem::exists( (test_path/outConvergenceMap).native() )) {
    BOOST_CHECK(true);
  }
  BOOST_CHECK(status);*/
  std::cout <<"getNoisedConvergenceMap_test" << std::endl;
  status = CartesainAlgo.getNoisyConvergenceMap(test_path, shearMap, outConvergenceMap);
  //if (boost::filesystem::exists( (test_path/outConvergenceMapXML).native() )) {
    BOOST_CHECK(true);
  //}
  BOOST_CHECK(status);

  if (params.getSigmaGauss() != 0) {
    std::cout <<"getDenoisedConvergenceMap_test" << std::endl;
    status = CartesainAlgo.getDenoisedConvergenceMap(test_path, shearMap, outConvergenceMap);
    //if (boost::filesystem::exists( (test_path/outConvergenceMapXML).native() )) {
      BOOST_CHECK(true);
   // }
    BOOST_CHECK(status);
  }

  std::cout <<"PerformInverseKSMassMapping_test" << std::endl;
  boost::filesystem::path outShearMap= "ShearMapKSFile.fits";
  status = CartesainAlgo.performInverseKSMassMapping((test_path/outConvergenceMap).native(),
                                       (datadir/outShearMap).native());
  if (boost::filesystem::exists((datadir/outShearMap).native() )) {
    BOOST_CHECK(true);
  }
  BOOST_CHECK(status);

  auto testRSParamFile = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/CparamRS_test.xml");
  CartesianParam paramsRS;
  readParameterFile(testRSParamFile, paramsRS);
  Euclid::WeakLensing::TwoDMass::CartesianKS::CartesianAlgoKS CartesainAlgoRS(paramsRS);

  std::cout <<"PerformReducedShearComputation_test" << std::endl;
  if (params.getNItReducedShear() != 0) {
   status = CartesainAlgoRS.performReducedShear(shearMap, outConvergenceMap, test_path, outConvergenceMapXML);
  }
  if (boost::filesystem::exists((test_path/outConvergenceMapXML).native() )) {
    BOOST_CHECK(true);
  }
  BOOST_CHECK(status);

 }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
