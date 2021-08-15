/**
 * @file tests/src/DmInput_test.cpp
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
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"
#include "ElementsKernel/Auxiliary.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "ElementsServices/DataSync.h"

#include <fstream>
#include <string>

using namespace Euclid::WeakLensing::TwoDMass;
namespace fs = boost::filesystem;
using namespace ElementsServices::DataSync;

//-----------------------------------------------------------------------------

struct DmInputDataSyncFixture {
  DataSync sync; // This is the synchronizer
  // These are just shortcuts
  path filename, filename_ksb, filename_regauss, filename_momentml;
  path filename_le2, filename_cluster;
  DmInputDataSyncFixture () :
      sync(
          // Here is the connection configuration file
          "LE3_2D_MASS_WL_UTILITIES/datasync_webdav.conf",
           // Here is the dependency configuration file
          "LE3_2D_MASS_WL_UTILITIES/test_file_list.txt"),
         filename_ksb(sync.absolutePath("KSBCatalog.xml")),
         filename_regauss(sync.absolutePath("RegaussCatalog.xml")),
         filename_momentml(sync.absolutePath("MomentsMLCatalog.xml")),
         filename_le2(sync.absolutePath("InputLE2Catalog.xml")),
         filename_cluster(sync.absolutePath("ClusterCatalog.xml")),
         filename(sync.absolutePath("SampleInXML.xml"))
 {
    sync.download();
  }
};
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (DmInput_test, DmInputDataSyncFixture)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_SUITE (DmInput_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(readLensMCCatalogXMLFile_test ) {
  std::cout << "-- Read_XMLFile: readLensMCCatalogXMLFile"<<std::endl;
  try{
  //std::string file = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/SampleInXML.xml").generic_string();
   if (true == checkFileType(filename, Euclid::WeakLensing::TwoDMass::signXML)) {
     DmInput in_xml = DmInput::readLensMCCatalogXMLFile(filename);
     BOOST_CHECK( true);

     std::cout << "-- Read_XMLFile: getFitsCatalogFilename"<<std::endl;

     std::string inputCatalog = (in_xml.getFitsCatalogFilename()).string();
     if (false == inputCatalog.empty() && (true == Euclid::WeakLensing::TwoDMass::checkFileType(inputCatalog,
                                              Euclid::WeakLensing::TwoDMass::signFITS))) {
      BOOST_CHECK( true);
     } } } catch (...) {
     BOOST_FAIL("Exception in getFitsCatalogFilename");
  }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(readClusterCatalogXMLFile_test ) {
  std::cout << "-- Read_XMLFile: readClusterCatalogXMLFile"<<std::endl;
  try{

   if (true == checkFileType(filename_cluster, Euclid::WeakLensing::TwoDMass::signXML)) {
     DmInput in_xml = DmInput::readClusterCatalogXMLFile(filename_cluster);
     BOOST_CHECK( true);
   } } catch (...) {
     BOOST_FAIL("Exception in readClusterCatalogXMLFile");
  }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(readRegaussCatalogXMLFile_test ) {
  std::cout << "-- Read_XMLFile: readRegaussCatalogXMLFile"<<std::endl;
  try{
  //std::string file = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/RegaussCatalog.xml").generic_string();
  //const fs::path filename = file;
   if (true == checkFileType(filename_regauss, Euclid::WeakLensing::TwoDMass::signXML)) {
     DmInput in_xml = DmInput::readRegaussCatalogXMLFile(filename_regauss);
     BOOST_CHECK( in_xml.getFitsCatalogFilename().string() == std::string("Fake_LE2_Cat.fits"));
   } } catch (...) {
     BOOST_FAIL("Exception in readRegaussCatalogXMLFile");
  }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(readKSBCatalogXMLFile_test ) {
  std::cout << "-- Read_XMLFile: readKSBCatalogXMLFile"<<std::endl;
  try{
  //std::string file = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/KSBCatalog.xml").generic_string();
  //const fs::path filename = file;
   if (true == checkFileType(filename_ksb, Euclid::WeakLensing::TwoDMass::signXML)) {
     DmInput in_xml = DmInput::readKSBCatalogXMLFile(filename_ksb);
     BOOST_CHECK( in_xml.getFitsCatalogFilename().string() == std::string("Fake_LE2_Cat.fits"));
   } } catch (...) {
     BOOST_FAIL("Exception in readKSBCatalogXMLFile");
  }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(readMomentsMLCatalogXMLFile_test ) {
  std::cout << "-- Read_XMLFile: readMomentsMLCatalogXMLFile"<<std::endl;
  try{
  //std::string file = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/MomentsMLCatalog.xml").generic_string();
  //const fs::path filename = file;
   if (true == checkFileType(filename_momentml, Euclid::WeakLensing::TwoDMass::signXML)) {
     DmInput in_xml = DmInput::readMomentsMLCatalogXMLFile(filename_momentml);
     BOOST_CHECK( in_xml.getFitsCatalogFilename().string() == std::string("Fake_LE2_Cat.fits"));
   } } catch (...) {
     BOOST_FAIL("Exception in readMomentsMLCatalogXMLFile");
  }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(readLE2CatalogXMLFile_test ) {
  std::cout << "-- Read_XMLFile: readLE2CatalogXMLFile"<<std::endl;
  try{
   if (true == checkFileType(filename_le2, Euclid::WeakLensing::TwoDMass::signXML)) {
     DmInput in_xml = DmInput::readLE2CatalogXMLFile(filename_le2);
     BOOST_CHECK( in_xml.getFitsCatalogFilename().string() == std::string("lensmcCat.fits"));
   } } catch (...) {
     BOOST_FAIL("Exception in readLE2CatalogXMLFile");
  }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getMethodType_test ) {
  std::cout <<"-- DmInput: getMethodType_test"<<std::endl;
  try
  {
   if (true == checkFileType(filename_le2, Euclid::WeakLensing::TwoDMass::signXML)) {
     DmInput in_xml = DmInput::readLE2CatalogXMLFile(filename_le2);
     BOOST_CHECK( in_xml.getMethodType() == "LENSMC");
   }
  }
  catch (...)
  {
    BOOST_CHECK("Exception in getMethodType");
  }

}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
