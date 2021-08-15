/**
 * @file tests/src/ReadCatalog_test.cpp
 * @date 05/15/21
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
#include "ElementsKernel/Logging.h"
#include "LE3_2D_MASS_WL_UTILITIES/ReadCatalog.h"
#include "ElementsKernel/Auxiliary.h"
#include "ElementsServices/DataSync.h"

#include <ios>
#include <sstream>
#include <iostream>
#include <vector>

using std::string;
using boost::filesystem::path;
using namespace ElementsServices::DataSync;
using namespace Euclid::WeakLensing::TwoDMass;
using namespace Euclid;
using namespace FitsIO;

static Elements::Logging logger = Elements::Logging::getLogger("ReadCatalog_test");
//-----------------------------------------------------------------------------

struct ReadCatalogDataSyncFixture {
  DataSync sync; // This is the synchronizer
  // These are just shortcuts
  path fitsFilePath, xmlFilePath;
  ReadCatalogDataSyncFixture () :
      sync(
          // Here is the connection configuration file
          "LE3_2D_MASS_WL_UTILITIES/datasync_webdav.conf",
           // Here is the dependency configuration file
          "LE3_2D_MASS_WL_UTILITIES/test_file_list.txt"),
         fitsFilePath(sync.absolutePath("data/lensmcCat.fits")),
         xmlFilePath(sync.absolutePath("InputLE2Catalog.xml"))
 {
    sync.download();
  }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (ReadCatalog_test, ReadCatalogDataSyncFixture)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_SUITE (ReadCatalog_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( readShearCatalogFITS_test ) {
  logger.info() << "-- ReadCatalog: readShearCatalogFITS_test";
  try
  {
    // Retrieve fits file directory location and filemane
    //auto fitsFilePath = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/data/lensmcCat.fits");

    // fitsFilePath -> "/workdir/data/catalog.fits"
    auto inputCatalogueFileName = fitsFilePath.filename();

    // fitsFilePath.filename() -> "catalog.fits"
    // fitsFilePath.branch_path() -> "/workdir/data"
    // fitsFilePath.branch_path().branch_path() -> "/workdir"
    auto dataDirectoryPath = fitsFilePath.branch_path();
    auto workdirPath = dataDirectoryPath.branch_path();

    logger.info() << "Input catalogue filename: " << inputCatalogueFileName.string();
    logger.info() << "Input catalogue path: " << workdirPath.string();

    // Read the Shear catalog
    std::vector<std::vector<double>> testCatData;
    ReadCatalog read;
    read.readShearCatalog(workdirPath, inputCatalogueFileName, testCatData);

    MefFile b(fitsFilePath.native(), MefFile::Permission::Read);
    int expectedCatalogSize = b.access<BintableHdu>(1).readRowCount();
    logger.info() << "expected catalogue row count: " << expectedCatalogSize;

    int catalogSize = testCatData[0].size();
    logger.info() << "Input catalogue row count: " << catalogSize;
    BOOST_CHECK(catalogSize == expectedCatalogSize);
  }
  catch (...)
  {
    BOOST_CHECK(false);
  }

}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( readShearCatalogXML_test ) {
  logger.info() << "-- ReadCatalog: readShearCatalogXML_test";
  try
  {
    // Retrieve fits file directory location and filemane
    //auto xmlFilePath = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/InputLE2Catalog.xml");
    //auto fitsFilePath = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/data/lensmcCat.fits");

    auto inputCatalogueFileName = xmlFilePath.filename();

    // xmlFilePath.filename() -> "catalog.fits"
    // xmlFilePath.branch_path() -> "/workdir/data"
    // xmlFilePath.branch_path().branch_path() -> "/workdir"
    auto dataDirectoryPath = xmlFilePath.branch_path();
    auto workdirPath = dataDirectoryPath.branch_path();

    logger.info() << "Input catalogue filename: " << inputCatalogueFileName.string();
    logger.info() << "Input catalogue path: " << dataDirectoryPath.string();

    // Read the Shear catalog
    std::vector<std::vector<double>> testCatData;
    ReadCatalog read;
    read.readShearCatalog(dataDirectoryPath, inputCatalogueFileName, testCatData);

    MefFile b(fitsFilePath.native(), MefFile::Permission::Read);
    int expectedCatalogSize = b.access<BintableHdu>(1).readRowCount();
    logger.info() << "expected catalogue row count: " << expectedCatalogSize;

    int catalogSize = testCatData[0].size();
    logger.info() << "Input catalogue row count: " << catalogSize;
    BOOST_CHECK(catalogSize == expectedCatalogSize);
  }
  catch (...)
  {
    BOOST_CHECK(false);
  }

}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( getCatalogReadType_test ) {
  logger.info() <<"-- ReadCatalog: getCatalogReadType_test";
  try
  {
    // Retrieve fits file directory location and filemane
    //auto xmlFilePath = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/InputLE2Catalog.xml");

    auto inputCatalogueFileName = xmlFilePath.filename();

    // xmlFilePath.filename() -> "catalog.fits"
    // xmlFilePath.branch_path() -> "/workdir/data"
    // xmlFilePath.branch_path().branch_path() -> "/workdir"
    auto dataDirectoryPath = xmlFilePath.branch_path();
    auto workdirPath = dataDirectoryPath.branch_path();

    logger.info() << "Input catalogue filename: " << inputCatalogueFileName.string();
    logger.info() << "Input catalogue path: " << dataDirectoryPath.string();

    // Read the Shear catalog
    std::vector<std::vector<double>> testCatData;
    ReadCatalog read;
    read.readShearCatalog(dataDirectoryPath, inputCatalogueFileName, testCatData);

    std::string catType = read.getCatalogReadType();
    BOOST_CHECK(catType == "LENSMC");
  }
  catch (...)
  {
    BOOST_CHECK(false);
  }

}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( readShearCatalogXMLException_test ) {
  logger.info() << "-- ReadCatalog: readShearCatalogXMLException_test";
    // Retrieve fits file directory location and filemane
    auto fakeFilePath = fs::path("wrong_path/wrong_file.xml");
    auto inputCatalogueFileName = fakeFilePath.filename();
    auto dataDirectoryPath = fakeFilePath.branch_path();
    // Read the Shear catalog
    ReadCatalog read;
    std::vector<std::vector<double>> testCatData;
    BOOST_CHECK_THROW( read.readShearCatalog(dataDirectoryPath, inputCatalogueFileName, testCatData),
                                                                                Elements::Exception);
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( readShearCatalogFitsException_test ) {
  logger.info() << "-- ReadCatalog: readShearCatalogFitsException_test";
    // Retrieve fits file directory location and filemane
    auto fakeFilePath = fs::path("wrong_path/wrong_file.fits");
    auto inputCatalogueFileName = fakeFilePath.filename();
    auto dataDirectoryPath = fakeFilePath.branch_path();
    // Read the Shear catalog
    ReadCatalog read;
    std::vector<std::vector<double>> testCatData;
    BOOST_CHECK_THROW( read.readShearCatalog(dataDirectoryPath, inputCatalogueFileName, testCatData),
                                                                                Elements::Exception);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( readShearCatalogEmptyCatalogFileName_test ) {
  logger.info() << "-- ReadCatalog: readShearCatalogEmptyCatalogFileName_test";
    // Retrieve fits file directory location and filemane
    auto fakeFilePath = fs::path("wrong_path/");
    auto inputCatalogueFileName = fakeFilePath.filename();
    auto dataDirectoryPath = fakeFilePath.branch_path();
    // Read the Shear catalog
    ReadCatalog read;
    std::vector<std::vector<double>> testCatData;
    BOOST_CHECK_THROW( read.readShearCatalog(dataDirectoryPath, inputCatalogueFileName, testCatData),
                                                                                Elements::Exception);
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
