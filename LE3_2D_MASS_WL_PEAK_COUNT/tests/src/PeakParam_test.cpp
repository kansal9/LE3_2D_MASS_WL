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
 * @file tests/src/PeakParam_test.cpp
 * @date 09/17/19
 * @author user
 */

#include <boost/test/unit_test.hpp>
#include "LE3_2D_MASS_WL_PEAK_COUNT/PeakParam.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "ElementsKernel/Auxiliary.h"
#include "ElementsServices/DataSync.h"
#include <boost/filesystem.hpp>
#include <ios>
#include <sstream>
#include <iostream>

using namespace LE3_2D_MASS_WL_PEAK_COUNT;

using namespace Euclid::WeakLensing::TwoDMass;
using namespace ElementsServices::DataSync;
using std::string;
using boost::filesystem::path;

//-----------------------------------------------------------------------------

struct PeakParamDataSyncFixture {
  DataSync sync; // This is the synchronizer
  // These are just shortcuts
  path xmlPConvFileName, xmlMAFileName;
  PeakParamDataSyncFixture () :
      sync(
          // Here is the connection configuration file
          "LE3_2D_MASS_WL_PEAK_COUNT/datasync_webdav.conf",
           // Here is the dependency configuration file
          "LE3_2D_MASS_WL_PEAK_COUNT/test_file_list.txt"),
         xmlPConvFileName(sync.absolutePath("PConvparam_test.xml")),
         xmlMAFileName(sync.absolutePath("PMassAp_test.xml"))
 {
    sync.download();
  }
};

//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_SUITE (PeakParam_test)
BOOST_FIXTURE_TEST_SUITE (PeakParam_test, PeakParamDataSyncFixture)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( Pparam_test ) {
 std::cout<<" --> PeakParam: Overall_test"<<std::endl;
 int nbScale(5);
 double MinPeakThresh(0.01);
 std::vector <double> ApPeakRadius(0.0);

 PeakParam Pparam( MinPeakThresh, ApPeakRadius, nbScale);
 BOOST_CHECK_EQUAL(Pparam.getnbScales(), 5);
 BOOST_CHECK_CLOSE(Pparam.getMinPeakThresh(), MinPeakThresh, 0.0001);
}
//-----------------------------------------------------------------------------
 BOOST_FIXTURE_TEST_CASE(xmlPConvFile_test, PeakParamDataSyncFixture) {
  std::cout << "-- PARAM: .xml file for Wavelet";
  PeakParam Param;
   if (true == fileHasField(xmlPConvFileName, "DpdTwoDMassParamsPeakCatalogConvergence")) {
    std::cout<<"Parameter file is for input Convergence Maps.."<<std::endl;
    Param.readPeakConvXML(xmlPConvFileName.native());
      int scale=Param.getnbScales();
      BOOST_CHECK_EQUAL(scale, 5);
   }
}

//-----------------------------------------------------------------------------
 BOOST_FIXTURE_TEST_CASE(xmlMAFile_test, PeakParamDataSyncFixture) {
  std::cout << "-- PARAM: .xml file for Mass Aperture";
  PeakParam Param;
   if (true == fileHasField(xmlMAFileName, "DpdTwoDMassParamsPeakCatalogMassAperture")) {
    std::cout<<"Parameter file is for Mass Aperture.."<<std::endl;
    Param.readPeakMassApertXML(xmlMAFileName.native());
    std::vector<double> nPScale=Param.getApPeakRadius();
    BOOST_CHECK_EQUAL(nPScale[0], 2);
   }
}

//-----------------------------------------------------------------------------
 BOOST_FIXTURE_TEST_CASE(readPeakParamFile_test, PeakParamDataSyncFixture) {
  std::cout << "-- PARAM: readPeakParamFile_test";
  PeakParam Param;
  readPeakParamFile (fs::path{xmlMAFileName}, Param);
  std::vector<double> nPScale=Param.getApPeakRadius();
  BOOST_CHECK_EQUAL(nPScale[0], 2);
}
//-----------------------------------------------------------------------------
 BOOST_FIXTURE_TEST_CASE(readPeakParamFile_test2, PeakParamDataSyncFixture) {
  std::cout << "-- PARAM: readPeakParamFile_test2";
  PeakParam Param;
  readPeakParamFile (fs::path{xmlPConvFileName}, Param);
  BOOST_CHECK_EQUAL(Param.getnbScales(), 5);
}
//-----------------------------------------------------------------------------
 BOOST_FIXTURE_TEST_CASE(readPeakParamFile_test3, PeakParamDataSyncFixture) {
  std::cout << "-- PARAM: readPeakParamFile_test3";
  auto fakeFilePath = fs::path("wrong_path/wrong_file.ini");
  PeakParam Param;
  //readPeakParamFile (fs::path{fakeFilePath}, Param);
  BOOST_CHECK_THROW(readPeakParamFile (fs::path{fakeFilePath}, Param), Elements::Exception);
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
