/**
 * @file tests/src/SplitCatalog_test.cpp
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
#define BOOST_TEST_MODULE SplitCatalog_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "ElementsServices/DataSync.h"

#include "LE3_2D_MASS_WL_CARTESIAN/SplitCatalog.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "ElementsKernel/Auxiliary.h"
#include "ElementsKernel/Temporary.h"
#include <vector>
#include <memory>
#include <iostream>
#include <dirent.h>
#include <ios>
#include <sstream>
//#include "HeaderProvider/FilenameProvider.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
using namespace LE3_2D_MASS_WL_CARTESIAN;
using namespace Euclid::WeakLensing::TwoDMass;
using LE3_2D_MASS_WL_CARTESIAN::CartesianParam;
using namespace ElementsServices::DataSync;
using std::string;
namespace fs = boost::filesystem;

boost::filesystem::path test_path;

//-----------------------------------------------------------------------------

struct SplitCatalogDataSyncFixture {
  DataSync sync; // This is the synchronizer
  // These are just shortcuts
  path testParamFile, testParamFile1, Incatalog;
  SplitCatalogDataSyncFixture () :
      sync(
          // Here is the connection configuration file
          "LE3_2D_MASS_WL_CARTESIAN/datasync_webdav.conf",
           // Here is the dependency configuration file
          "LE3_2D_MASS_WL_CARTESIAN/test_file_list.txt"),
         testParamFile(sync.absolutePath("Cparam_test.xml")),
         testParamFile1(sync.absolutePath("Cparam_test2.xml")),
         Incatalog(sync.absolutePath("data/lensmcCat.fits"))
 {
    sync.download();
  }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (SplitCatalog_test, SplitCatalogDataSyncFixture)
//BOOST_AUTO_TEST_SUITE (SplitCatalog_test)

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( GetEqualBinSplittedCatalogs_test) {

   CartesianParam params;
   // Case 1: check file is of XML format
   if (true == checkFileType(testParamFile, Euclid::WeakLensing::TwoDMass::signXML)) {
    // read XML, get parameter file type
    // check if it's to get a patch
    if (true == fileHasField(testParamFile, "DpdTwoDMassParamsConvergencePatch")) {
      params.ReadConvPatchXMLFile (testParamFile.native());
    }
   }
  CatalogData read;
  std::vector<std::vector<double> > catData;
  std::string m_catType="LENSMC";
  //fs::path p = fs::path(Incatalog).parent_path();
  //fs::path f = fs::path(Incatalog).filename();
  read.readCatalog(Incatalog.native(), catData);
  std::cout << "-- SplitCatalog: GetEqualBinSplittedCatalogs"<<std::endl;
  //SplitCatalog split(catData, params);

  Elements::TempDir one{"Split-Catalogue-%%%%%%"};
  test_path = one.path();
  std::vector<std::string> Filenames_Cat;

  for(int i =0; i<params.getnbZBins(); i++) {
    boost::filesystem::path subCatalogName = fs::path ( "CatZ_Test_subcatalogue_Z0" +
                                 std::to_string(i) + ".fits");
     Filenames_Cat.push_back((subCatalogName).string());
  }
  Filenames_Cat.shrink_to_fit();
  try {
  SplitCatalog split(catData, params, m_catType);
  // Save these sub-catalogs
  split.writeSubCatalogs(test_path, Filenames_Cat);
   for(int i =0; i<params.getnbZBins(); i++) {
    BOOST_CHECK(boost::filesystem::exists(test_path.string() + "/" + Filenames_Cat[i]));
   }
  } catch (...) {
    BOOST_FAIL("Exception in GetEqualBinSplittedCatalogs_test");
  }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( GetSplittedCatalogs_test ) {

   CartesianParam params;
   // Case 1: check file is of XML format
   if (true == checkFileType(testParamFile1, Euclid::WeakLensing::TwoDMass::signXML)) {
    // read XML, get parameter file type
    // check if it's to get a patch
    if (true == fileHasField(testParamFile1, "DpdTwoDMassParamsConvergencePatch")) {
      params.ReadConvPatchXMLFile (testParamFile1.native());
    }
   }
  CatalogData read;
  std::vector<std::vector<double> > catData;
  //fs::path p = fs::path(Incatalog).parent_path();
  //fs::path f = fs::path(Incatalog).filename();
  //catData = read.getCatalogData(p, f);
  read.readCatalog(Incatalog.native(), catData);
  std::string m_catType=read.getCatalogType();
  std::cout << "-- SplitCatalog: GetSplittedCatalogs"<<std::endl;
  //SplitCatalog split(catData, params);

  Elements::TempDir one{"Split-Catalogue-%%%%%%"};
  test_path = one.path();
  std::vector<std::string> Filenames_Cat;

  for(int i =0; i<params.getnbZBins(); i++) {
    boost::filesystem::path subCatalogName = fs::path ( "CatZ_Test_subcatalogue_Z0" +
                                 std::to_string(i) + ".fits");
     Filenames_Cat.push_back((subCatalogName).string());
  }
  Filenames_Cat.shrink_to_fit();
  try {
  SplitCatalog split(catData, params, m_catType);
  // Save these sub-catalogs
  split.writeSubCatalogs(test_path, Filenames_Cat);
   for(int i =0; i<params.getnbZBins(); i++) {
    BOOST_CHECK(boost::filesystem::exists(test_path.string() + "/" + Filenames_Cat[i]));
   }
  } catch (...) {
    BOOST_FAIL("Exception in GetSplittedCatalogs_test");
  }
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
