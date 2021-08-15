/**
 * @file tests/src/MassMapping_test.cpp
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
//#define BOOST_TEST_MODULE MassMapping_test
//#define private public
//#define protected public
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "ElementsKernel/Auxiliary.h"
#include "ElementsKernel/Temporary.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MassMapping.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MapMaker.h"
#include "ElementsKernel/Logging.h"
#include "ElementsServices/DataSync.h"
#include <ios>
#include <sstream>
#include <iostream>
using namespace Euclid::WeakLensing::TwoDMass;
using namespace LE3_2D_MASS_WL_CARTESIAN;
using LE3_2D_MASS_WL_CARTESIAN::CartesianParam;
using namespace ElementsServices::DataSync;
using std::string;

 // handle on created path names
 boost::filesystem::path test_path;

static Elements::Logging logger = Elements::Logging::getLogger("MassMapping_test");

//-----------------------------------------------------------------------------

struct MassMappingDataSyncFixture {
  DataSync sync; // This is the synchronizer
  // These are just shortcuts
  path testParamFile, shearMapFile;
  LE3_2D_MASS_WL_CARTESIAN::ShearMap *m_ShearMap = nullptr;
  LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap *m_ConvergenceMap = nullptr;
  MassMappingDataSyncFixture () :
      sync(
          // Here is the connection configuration file
          "LE3_2D_MASS_WL_CARTESIAN/datasync_webdav.conf",
           // Here is the dependency configuration file
          "LE3_2D_MASS_WL_CARTESIAN/test_file_list.txt"),
         testParamFile(sync.absolutePath("Cparam_test.xml")),
         shearMapFile(sync.absolutePath("data/shearMap_test.fits"))
 {
    sync.download();
  }
};

//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_SUITE (MassMapping_test)
BOOST_FIXTURE_TEST_SUITE (MassMapping_test, MassMappingDataSyncFixture)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( get_SheartoConv_viseversa_test ) {
  logger.info() <<"get_shearMap_test";
  Elements::TempDir one;
  test_path = one.path();
  CartesianParam params;
   // Case 1: check file is of XML format
   if (true == checkFileType(testParamFile, Euclid::WeakLensing::TwoDMass::signXML)) {
    // read XML, get parameter file type
    // check if it's to get a patch
    if (true == fileHasField(testParamFile, "DpdTwoDMassParamsConvergencePatch")) {
      params.ReadConvPatchXMLFile (testParamFile.native());
    }
   }
  MassMapping mass(params);
  boost::filesystem::path convMap= "Test_convMapFile_MassMapping.fits";
  m_ConvergenceMap = mass.getSheartoConv(shearMapFile.native());
  m_ConvergenceMap->writeMap((test_path/convMap).native(), params);
  BOOST_CHECK(m_ConvergenceMap!=nullptr);

  m_ShearMap = mass.getConvtoShear((test_path/convMap).native());
  BOOST_CHECK(m_ShearMap!=nullptr);
  delete m_ShearMap;
  delete m_ConvergenceMap;
  m_ShearMap = nullptr;
  m_ConvergenceMap = nullptr;
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
