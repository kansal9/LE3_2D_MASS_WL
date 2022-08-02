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

#include "ElementsServices/DataSync.h"

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"

#include <fstream>
#include <string>

using namespace ElementsServices::DataSync;
using namespace LE3_2D_MASS_WL_UTILITIES;
namespace fs = boost::filesystem;

//-----------------------------------------------------------------------------
struct DmInputDataSyncFixture
{
    DmInput inp;
    DataSync sync;
    path filename_ksb, filename_regauss, filename_momentml, filename_le2,
         filename_cluster, filename_lensmc;
    DmInputDataSyncFixture() :
            sync("LE3_2D_MASS_WL_UTILITIES/datasync_webdav.conf",
                 "LE3_2D_MASS_WL_UTILITIES/test_file_list.txt"),
            filename_ksb(sync.absolutePath("KSBCatalog.xml")),
            filename_regauss(sync.absolutePath("RegaussCatalog.xml")),
            filename_momentml(sync.absolutePath("MomentsMLCatalog.xml")),
            filename_le2(sync.absolutePath("InputLE2Catalog.xml")),
            filename_cluster(sync.absolutePath("ClusterCatalog.xml")),
            filename_lensmc(sync.absolutePath("SampleInXML.xml"))
    {
        sync.download();
    }
};
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (DmInput_test, DmInputDataSyncFixture)
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( dmInput_le2catalog_test )
{
    inp.readXmlFile(filename_le2);
    std::cout << inp.getCatalogClassName() << " "
              << inp.getCatalogFilename() << std::endl;
    BOOST_CHECK(inp.getCatalogClassName() == "DpdWLLE2Catalog");
    BOOST_CHECK(inp.getCatalogFilename() == "lensmcCat.fits");
}

BOOST_AUTO_TEST_CASE( dmInput_ksbcatalog_test )
{
    inp.readXmlFile(filename_ksb);
    std::cout << inp.getCatalogClassName() << " "
              << inp.getCatalogFilename() << std::endl;
    BOOST_CHECK(inp.getCatalogClassName() == "DpdTwoDMassKSBCatalog");
    BOOST_CHECK(inp.getCatalogFilename() == "Fake_LE2_Cat.fits");
}

BOOST_AUTO_TEST_CASE( dmInput_regausscatalog_test )
{
    inp.readXmlFile(filename_regauss);
    std::cout << inp.getCatalogClassName() << " "
              << inp.getCatalogFilename() << std::endl;
    BOOST_CHECK(inp.getCatalogClassName() == "DpdTwoDMassRegaussCatalog");
    BOOST_CHECK(inp.getCatalogFilename() == "Fake_LE2_Cat.fits");
}

BOOST_AUTO_TEST_CASE( dmInput_momentmlcatalog_test )
{
    inp.readXmlFile(filename_momentml);
    std::cout << inp.getCatalogClassName() << " "
              << inp.getCatalogFilename() << std::endl;
    BOOST_CHECK(inp.getCatalogClassName() == "DpdTwoDMassMomentsMLCatalog");
    BOOST_CHECK(inp.getCatalogFilename() == "Fake_LE2_Cat.fits");
}

BOOST_AUTO_TEST_CASE( dmInput_clustercatalog_test )
{
    inp.readXmlFile(filename_cluster);
    std::cout << inp.getCatalogClassName() << " "
              << inp.getCatalogFilename() << std::endl;
    BOOST_CHECK(inp.getCatalogClassName() == "DpdTwoDMassClusterCatalog");
    BOOST_CHECK(inp.getCatalogFilename() == "Cluster_Catalog.fits");
}

BOOST_AUTO_TEST_CASE( dmInput_lensmccatalog_test )
{
    inp.readXmlFile(filename_lensmc);
    std::cout << inp.getCatalogClassName() << " "
              << inp.getCatalogFilename() << std::endl;
    BOOST_CHECK(inp.getCatalogClassName() == "DpdTwoDMassLensMCCatalog");
    BOOST_CHECK(inp.getCatalogFilename() == "lensmcCat.fits");
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
