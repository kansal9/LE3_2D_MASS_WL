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
#include "ElementsServices/DataSync.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"

#include <ios>
#include <fstream>
#include <string>
#include <vector>

using boost::filesystem::path;
using namespace ElementsServices::DataSync;
using namespace LE3_2D_MASS_WL_UTILITIES;

static Elements::Logging logger = Elements::Logging::getLogger(
        "CatalogData_test");

//-----------------------------------------------------------------------------

struct CatalogDataSyncFixture
{
    DataSync sync;
    CatalogData cat;
    path fitsFilePath_lensmc, filename_lensmc, filename_ksb, filename_regauss,
         filename_momentml, file_le2Cat, filename_le2, filename_cluster,
         fitsFilename_cluster;
    CatalogDataSyncFixture() :
        sync("LE3_2D_MASS_WL_UTILITIES/datasync_webdav.conf",
             "LE3_2D_MASS_WL_UTILITIES/test_file_list.txt"),
        fitsFilePath_lensmc(sync.absolutePath("data/lensmcCat.fits")),
        filename_lensmc(sync.absolutePath("SampleInXML.xml")),
        filename_ksb(sync.absolutePath("KSBCatalog.xml")),
        filename_regauss(sync.absolutePath("RegaussCatalog.xml")),
        filename_momentml(sync.absolutePath("MomentsMLCatalog.xml")),
        file_le2Cat(sync.absolutePath("data/Fake_LE2_Cat.fits")),
        filename_le2(sync.absolutePath("InputLE2Catalog.xml")),
        filename_cluster(sync.absolutePath("ClusterCatalog.xml")),
        fitsFilename_cluster(sync.absolutePath("data/Cluster_Catalog.fits"))
    {
        sync.download();
    }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (CatalogData_test, CatalogDataSyncFixture)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( getCatalog_lensmc_test )
{
    logger.info() << "-- CatalogData: getCatalogData_LENSMC_test";
    fs::path workdir = filename_lensmc.parent_path();
    fs::path catalogFileName = filename_lensmc.filename();
    logger.info() << "Workdir " << workdir;
    logger.info() << "Filename " << catalogFileName;
    cat.getCatalogData(workdir, catalogFileName);
    logger.info() << "Input catalogue size: " << cat.getNentries();
    BOOST_CHECK_CLOSE(cat["dec"](0), 26.553995132446289, 0.0001);
    BOOST_CHECK(cat.getNentries() == 344749);
    BOOST_CHECK(cat.getNcols() == 7);
    BOOST_CHECK(cat.getCatalogType() == "LENSMC");
    BOOST_CHECK(cat.getCatalogFitsName() == fitsFilePath_lensmc.filename());

}

BOOST_AUTO_TEST_CASE( getCatalog_le2_test )
{
    logger.info() << "-- CatalogData: getCatalogData_LE2_test";
    fs::path workdir = file_le2Cat.parent_path().parent_path();
    fs::path catalogFileName = file_le2Cat.filename();
    logger.info() << "Workdir " << workdir;
    logger.info() << "Filename " << catalogFileName;
    cat.getCatalogData(workdir, catalogFileName);
    logger.info() << "Input catalogue size: " << cat.getNentries();
    BOOST_CHECK_CLOSE(cat["dec"](0), 25.1138706207275, 0.0001);
    BOOST_CHECK(cat.getNentries() == 2743);
    BOOST_CHECK(cat.getNcols() == 7);
    BOOST_CHECK(cat.getCatalogType() == "LENSMC");
    BOOST_CHECK(cat.getCatalogFitsName() == file_le2Cat.filename());
}

BOOST_AUTO_TEST_CASE( copyCtor_test )
{
    fs::path workdir = file_le2Cat.parent_path().parent_path();
    fs::path catalogFileName = file_le2Cat.filename();
    cat.getCatalogData(workdir, catalogFileName);
    CatalogData copycat(cat);
    BOOST_CHECK_CLOSE(copycat["dec"](0), 25.1138706207275, 0.0001);
    BOOST_CHECK(copycat.getNentries() == 2743);
}

BOOST_AUTO_TEST_CASE( fillTest_test )
{
    CatalogData catFillTest;
    int N = 10;
    catFillTest.fillTest(N);
    BOOST_CHECK(catFillTest.getNentries() == N);
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
