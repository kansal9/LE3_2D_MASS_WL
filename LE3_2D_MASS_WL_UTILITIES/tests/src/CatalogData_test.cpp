/**
 * @file tests/src/CatalogData_test.cpp
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
#include "ElementsKernel/Logging.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "ElementsKernel/Auxiliary.h"

#include "ElementsServices/DataSync.h"

#include <ios>
#include <fstream>
#include <string>
#include <vector>

using std::string;
using boost::filesystem::path;
using namespace ElementsServices::DataSync;
using namespace Euclid::WeakLensing::TwoDMass;

static Elements::Logging logger = Elements::Logging::getLogger("CatalogData_test");

//-----------------------------------------------------------------------------

struct CatalogDataSyncFixture {
  DataSync sync; // This is the synchronizer
  // These are just shortcuts
  path fitsFilePath_lensmc, filename_lensmc, filename_ksb, filename_regauss, filename_momentml, file_le2Cat;
  path filename_le2, filename_cluster, fitsFilename_cluster;
  CatalogDataSyncFixture () :
      sync(
          // Here is the connection configuration file
          "LE3_2D_MASS_WL_UTILITIES/datasync_webdav.conf",
           // Here is the dependency configuration file
          "LE3_2D_MASS_WL_UTILITIES/test_file_list.txt"),
         fitsFilePath_lensmc(sync.absolutePath("data/lensmcCat.fits")),
         filename_ksb(sync.absolutePath("KSBCatalog.xml")),
         filename_regauss(sync.absolutePath("RegaussCatalog.xml")),
         filename_momentml(sync.absolutePath("MomentsMLCatalog.xml")),
         filename_le2(sync.absolutePath("InputLE2Catalog.xml")),
         filename_cluster(sync.absolutePath("ClusterCatalog.xml")),
         fitsFilename_cluster(sync.absolutePath("data/Cluster_Catalog.fits")),
         file_le2Cat(sync.absolutePath("data/Fake_LE2_Cat.fits")),
         filename_lensmc(sync.absolutePath("SampleInXML.xml"))
 {
    sync.download();
  }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (CatalogData_test, CatalogDataSyncFixture)
//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_SUITE (CatalogData_test)

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getCatalogData_LENSMC_test ) {
  logger.info() << "-- CatalogData: getCatalogData_LENSMC_test";
  try
  {
    //auto filename = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/SampleInXML.xml");
    //auto fitsFilePath = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/data/lensmcCat.fits");
    std::vector<std::vector<double> > testCatData;
    auto incat = filename_lensmc.filename();
    auto test_path = filename_lensmc.branch_path();
    logger.info() << "Input catalogue filename: "<< incat;
    logger.info() << "Input catalogue path: "<< test_path;
    CatalogData cd;
    cd.getCatalogData(test_path, incat, testCatData);

    MefFile b(fitsFilePath_lensmc.native(), MefFile::Permission::Read);
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
BOOST_AUTO_TEST_CASE( getCatalogData_KSB_test ) {
  logger.info() << "-- CatalogData: getCatalogData_KSB_test";
  try
  {
    //auto filename = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/KSBCatalog.xml");
    //auto fitsFilePath = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/data/Fake_LE2_Cat.fits");
    std::vector<std::vector<double> > testCatData;
    auto incat = filename_ksb.filename();
    auto test_path = filename_ksb.branch_path();
    logger.info() << "Input catalogue filename: "<< incat;
    logger.info() << "Input catalogue path: "<< test_path;
    CatalogData cd;
    cd.getCatalogData(test_path, incat, testCatData);

    MefFile b(file_le2Cat.native(), MefFile::Permission::Read);
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
BOOST_AUTO_TEST_CASE( getCatalogData_REGAUSS_test ) {
  logger.info() << "-- CatalogData: getCatalogData_REGAUSS_test";
  try
  {
    //auto filename = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/RegaussCatalog.xml");
    //auto fitsFilePath = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/data/Fake_LE2_Cat.fits");
    std::vector<std::vector<double> > testCatData;
    auto incat = filename_regauss.filename();
    auto test_path = filename_regauss.branch_path();
    logger.info() << "Input catalogue filename: "<< incat;
    logger.info() << "Input catalogue path: "<< test_path;
    CatalogData cd;
    cd.getCatalogData(test_path, incat, testCatData);

    MefFile b(file_le2Cat.native(), MefFile::Permission::Read);
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
BOOST_AUTO_TEST_CASE( getCatalogData_MOMENTML_test ) {
  logger.info() << "-- CatalogData: getCatalogData_MOMENTML_test";
  try
  {
    //auto filename = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/MomentsMLCatalog.xml");
    //auto fitsFilePath = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/data/Fake_LE2_Cat.fits");
    std::vector<std::vector<double> > testCatData;
    auto incat = filename_momentml.filename();
    auto test_path = filename_momentml.branch_path();
    logger.info() << "Input catalogue filename: "<< incat;
    logger.info() << "Input catalogue path: "<< test_path;
    CatalogData cd;
    cd.getCatalogData(test_path, incat, testCatData);

    MefFile b(file_le2Cat.native(), MefFile::Permission::Read);
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
  logger.info() <<"-- CatalogData: getCatalogReadType_test";
  try
  {
    // Retrieve fits file directory location and filemane
    //auto xmlFilePath = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/InputLE2Catalog.xml");

    auto inputCatalogueFileName = filename_le2.filename();

    auto dataDirectoryPath = filename_le2.branch_path();

    logger.info() << "Input catalogue filename: " << inputCatalogueFileName.string();
    logger.info() << "Input catalogue path: " << dataDirectoryPath.string();

    // Read the Shear catalog
    std::vector<std::vector<double>> testCatData;
    CatalogData cd;
    cd.getCatalogData(dataDirectoryPath, inputCatalogueFileName, testCatData);

    std::string catType = cd.getCatalogType();
    BOOST_CHECK(catType == "LENSMC");
  }
  catch (...)
  {
    BOOST_CHECK(false);
  }

}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( readShearCatalog_test ) {
  logger.info() <<"-- CatalogData: readShearCatalog_test";
  try
  {
    // Read the Shear catalog
    std::vector<std::vector<double>> testCatData;
    CatalogData cd;
    cd.readShearCatalog(file_le2Cat.native(), testCatData);

    MefFile b(file_le2Cat.native(), MefFile::Permission::Read);
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
BOOST_AUTO_TEST_CASE( readClusterCatalog_test ) {
  logger.info() <<"-- CatalogData: readClusterCatalog_test";
  try
  {
    // Retrieve fits file directory location and filemane
    //auto fitsfilepath = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/data/Cluster_Catalog.fits");
    // Read the Shear catalog
    std::vector<std::vector<double>> testCatData;
    CatalogData cd("CLUSTER");
    cd.readClusterCatalog(fitsFilename_cluster.native(), testCatData);

    MefFile b(fitsFilename_cluster.native(), MefFile::Permission::Read);
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
BOOST_AUTO_TEST_CASE( getCatalogData_Cluster_test ) {
  logger.info() << "-- CatalogData: getCatalogData_Cluster_test";
  try
  {
    //auto filename = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/ClusterCatalog.xml");
    //auto fitsFilePath = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/data/Cluster_Catalog.fits");
    std::vector<std::vector<double> > testCatData;
    auto incat = filename_cluster.filename();
    auto test_path = filename_cluster.branch_path();
    logger.info() << "Input catalogue filename: "<< incat;
    logger.info() << "Input catalogue path: "<< test_path;
    CatalogData cd("CLUSTER");
    cd.getCatalogData(test_path, incat, testCatData);

    MefFile b(fitsFilename_cluster.native(), MefFile::Permission::Read);
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
BOOST_AUTO_TEST_CASE( getCatalogFitsName_test ) {
  logger.info() <<"-- CatalogData: getCatalogFitsName_test";
  try
  {
    // Retrieve fits file directory location and filemane
    //auto xmlFilePath = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/ClusterCatalog.xml");
    //auto fitsFilePath = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/data/Cluster_Catalog.fits");

    auto inputCatalogueFileName = filename_cluster.filename();

    auto dataDirectoryPath = filename_cluster.branch_path();

    logger.info() << "Input catalogue filename: " << inputCatalogueFileName.string();
    logger.info() << "Input catalogue path: " << dataDirectoryPath.string();
    //auto fitsfilename = fitsFilePath.filename();
    // Read the Shear catalog
    std::vector<std::vector<double>> testCatData;
    CatalogData cd("CLUSTER");
    cd.getCatalogData(dataDirectoryPath, inputCatalogueFileName, testCatData);

    std::string catFitsFileName = cd.getCatalogFitsName();
    BOOST_CHECK(catFitsFileName == fitsFilename_cluster.string());
  }
  catch (...)
  {
    BOOST_CHECK(false);
  }

}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
