/**
 * @file tests/src/VisibilityMask_test.cpp
 * @date 06/16/21
 * @author Vanshika Kansal
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

#include "LE3_2D_MASS_WL_UTILITIES/VisibilityMask.h"

#include <ios>
#include <fstream>
#include <string>
#include <vector>

using boost::filesystem::path;
using namespace ElementsServices::DataSync;
using namespace LE3_2D_MASS_WL_UTILITIES;

static Elements::Logging logger = Elements::Logging::getLogger(
        "VisibilityMask_test");
//-----------------------------------------------------------------------------

struct VisibilityMaskDataSyncFixture
{
    DataSync sync; // This is the synchronizer
    // These are just shortcuts
    path fitsFilePath, filename;
    VisibilityMaskDataSyncFixture() :
            sync(// Here is the connection configuration file
                 "LE3_2D_MASS_WL_UTILITIES/datasync_webdav.conf",
                 // Here is the dependency configuration file
                 "LE3_2D_MASS_WL_UTILITIES/test_file_list.txt"),
            fitsFilePath(sync.absolutePath("data/mask_combined_256.fits")),
            filename(sync.absolutePath("vMask.xml"))
    {
        sync.download();
    }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (VisibilityMask_test, VisibilityMaskDataSyncFixture)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( readVisibilityMask_test )
{
    logger.info() << "-- VisibilityMask: readVisibilityMask_test";
    try
    {
        std::vector<std::vector<double> > testCatData;
        auto incat = filename.filename();
        auto test_path = filename.branch_path();
        logger.info() << "Input catalogue filename: " << incat;
        logger.info() << "Input catalogue path: " << test_path;
        VisibilityMask vk(256);
        vk.readVisibilityMask(test_path, incat, testCatData);

        MefFile b(fitsFilePath.native(), FileMode::Read);
        int expectedCatalogSize = b.access<BintableHdu>(1).readRowCount();
        logger.info() << "expected catalogue row count: "
                << expectedCatalogSize;

        int catalogSize = testCatData[0].size();
        logger.info() << "Input catalogue row count: " << catalogSize;
        BOOST_CHECK(catalogSize == expectedCatalogSize);
    } catch (...)
    {
        BOOST_CHECK(false);
    }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( getMaskData_test )
{
    logger.info() << "-- VisibilityMask: getMaskData_test";
    try
    {
        std::vector<std::vector<double> > testCatData;
        auto incat = filename.filename();
        auto test_path = filename.branch_path();
        logger.info() << "Input catalogue filename: " << incat;
        logger.info() << "Input catalogue path: " << test_path;
        VisibilityMask vk(256);
        vk.getMaskData(test_path, incat, testCatData);

        MefFile b(fitsFilePath.native(), FileMode::Read);
        int expectedCatalogSize = b.access<BintableHdu>(1).readRowCount();
        logger.info() << "expected catalogue row count: "
                << expectedCatalogSize;

        int catalogSize = testCatData[0].size();
        logger.info() << "Input catalogue row count: " << catalogSize;
        BOOST_CHECK(catalogSize == expectedCatalogSize);
    } catch (...)
    {
        BOOST_CHECK(false);
    }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( getMaskData2_test )
{
    logger.info() << "-- VisibilityMask: getMaskData2_test";
    try
    {
        //auto fitsFilePath = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/data/mask_combined_256.fits");
        std::vector<std::vector<double> > testCatData;

        VisibilityMask vk(256);
        vk.getMaskData(fitsFilePath.native(), testCatData);

        MefFile b(fitsFilePath.native(), FileMode::Read);
        int expectedCatalogSize = b.access<BintableHdu>(1).readRowCount();
        logger.info() << "expected catalogue row count: "
                << expectedCatalogSize;

        int catalogSize = testCatData[0].size();
        logger.info() << "Input catalogue row count: " << catalogSize;
        BOOST_CHECK(catalogSize == expectedCatalogSize);
    } catch (...)
    {
        BOOST_CHECK(false);
    }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( changeResolutionUpgrad_test )
{
    logger.info() << "-- VisibilityMask: changeResolutionUpgrad_test ";
    try
    {
        std::vector<std::vector<double> > testCatData;

        VisibilityMask vk(512);
        vk.getMaskData(fitsFilePath.native(), testCatData);

        int catalogSize = testCatData[0].size();
        logger.info() << "Input catalogue row count: " << catalogSize;
        BOOST_CHECK(catalogSize == 12 * 512 * 512);
    } catch (...)
    {
        BOOST_CHECK(false);
    }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( changeResolutionDowngrad_test )
{
    logger.info() << "-- VisibilityMask: changeResolutionDowngrad_test";
    try
    {
        std::vector<std::vector<double> > testCatData;

        VisibilityMask vk(128);
        vk.getMaskData(fitsFilePath.native(), testCatData);

        int catalogSize = testCatData[0].size();
        logger.info() << "Input catalogue row count: " << catalogSize;
        BOOST_CHECK(catalogSize == 12 * 128 * 128);
    } catch (...)
    {
        BOOST_CHECK(false);
    }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( getMaskDataCatalogXMLException_test )
{
    logger.info() << "-- VisibilityMask: getMaskDataCatalogXMLException_test";
    // Retrieve fits file directory location and filemane
    auto fakeFilePath = fs::path("wrong_path/wrong_file.xml");
    auto inputFileName = fakeFilePath.filename();
    auto dataDirectoryPath = fakeFilePath.branch_path();
    // Read Mask
    VisibilityMask vk(256);
    std::vector<std::vector<double>> testData;
    BOOST_CHECK_THROW(
            vk.readVisibilityMask(dataDirectoryPath, inputFileName, testData),
            Elements::Exception);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( getMaskDataCatalogFitsException_test )
{
    logger.info() << "-- VisibilityMask: getMaskDataCatalogFitsException_test";
    // Retrieve fits file directory location and filemane
    auto fakeFilePath = fs::path("wrong_path/wrong_file.fits");
    auto inputFileName = fakeFilePath.filename();
    auto dataDirectoryPath = fakeFilePath.branch_path();
    // Read Mask
    VisibilityMask vk(256);
    std::vector<std::vector<double>> testData;
    BOOST_CHECK_THROW(
            vk.readVisibilityMask(dataDirectoryPath, inputFileName, testData),
            Elements::Exception);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
