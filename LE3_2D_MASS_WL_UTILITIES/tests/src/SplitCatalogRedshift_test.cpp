
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
#define BOOST_TEST_MODULE SplitCatalogRedshift_test

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_UTILITIES/SplitCatalogRedshift.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"

#include "ElementsServices/DataSync.h"
#include "ElementsKernel/Temporary.h"

#include <vector>
#include <memory>
#include <iostream>
#include <dirent.h>
#include <ios>
#include <sstream>

using namespace LE3_2D_MASS_WL_UTILITIES;
using namespace ElementsServices::DataSync;

using LE3_2D_MASS_WL_CARTESIAN::CartesianParam;
using LE3_2D_MASS_WL_UTILITIES::SplitCatalogRedshift;

namespace fs = boost::filesystem;

boost::filesystem::path test_path;

//-----------------------------------------------------------------------------

struct SplitCatalogDataSyncFixture
{
    DataSync sync; // This is the synchronizer
    // These are just shortcuts
    path testParamFile, testParamFile1, Incatalog;
    SplitCatalogDataSyncFixture() :
            sync(// Here is the connection configuration file
                 "LE3_2D_MASS_WL_UTILITIES/datasync_webdav.conf",
                 // Here is the dependency configuration file
                 "LE3_2D_MASS_WL_UTILITIES/test_file_list.txt"),
             testParamFile(sync.absolutePath("Cparam_test.xml")),
             testParamFile1(sync.absolutePath("Cparam_test2.xml")),
             Incatalog(sync.absolutePath("data/lensmcCat.fits"))
    {
        sync.download();
    }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (SplitCatalogRedshift_test, SplitCatalogDataSyncFixture)
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getEqualBinSplittedCatalogs_test )
{
    CartesianParam params;
    CatalogData cat;
    fs::path workdir = Incatalog.parent_path().parent_path();
    fs::path cataloFileName = Incatalog.filename();
    cat.getCatalogData(workdir, cataloFileName);
    std::string m_catType = cat.getCatalogType();
    std::cout << "-- SplitCatalog: GetEqualBinSplittedCatalogs" << std::endl;

    if (checkFileType(testParamFile, signXML))
    {
        // read XML, get parameter file type
        // check if it's to get a patch
        if (fileHasField(testParamFile, "DpdTwoDMassParamsConvergencePatch"))
        {
            params.readConvPatchXMLFile(testParamFile.native(), cat);
        }
    }

    Elements::TempDir one { "Split-Catalogue-%%%%%%" };
    test_path = one.path();
    std::vector<fs::path> Filenames_Cat;

    for (int i = 0; i < params.getNbins(); i++)
    {
        boost::filesystem::path subCatalogName = fs::path(
                "CatZ_Test_subcatalogue_Z0" + std::to_string(i) + ".fits");
        Filenames_Cat.push_back(subCatalogName);
    }
    Filenames_Cat.shrink_to_fit();

    SplitCatalogRedshift split(cat, params, m_catType);
    // Save these sub-catalogs
    split.writeSubCatalogs(test_path, Filenames_Cat);
    for (int i = 0; i < params.getNbins(); i++)
    {
        BOOST_CHECK(
                boost::filesystem::exists(test_path / Filenames_Cat[i]));
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getSplittedCatalogs_test )
{
    CartesianParam params;

    CatalogData cat;
    fs::path workdir = Incatalog.parent_path().parent_path();
    fs::path cataloFileName = Incatalog.filename();
    cat.getCatalogData(workdir, cataloFileName);
    std::string m_catType = cat.getCatalogType();
    std::cout << "-- SplitCatalog: GetSplittedCatalogs" << std::endl;

    if (checkFileType(testParamFile1, signXML))
    {
        // read XML, get parameter file type
        // check if it's to get a patch
        if (fileHasField(testParamFile1, "DpdTwoDMassParamsConvergencePatch"))
        {
            params.readConvPatchXMLFile(testParamFile1.native(), cat);
        }
    }

    Elements::TempDir one { "Split-Catalogue-%%%%%%" };
    test_path = one.path();
    std::vector<fs::path> Filenames_Cat;

    for (int i = 0; i < params.getNbins(); i++)
    {
        boost::filesystem::path subCatalogName = fs::path(
                "CatZ_Test_subcatalogue_Z0" + std::to_string(i) + ".fits");
        Filenames_Cat.push_back(subCatalogName);
    }
    Filenames_Cat.shrink_to_fit();
    SplitCatalogRedshift split(cat, params, m_catType);
    // Save these sub-catalogs
    split.writeSubCatalogs(test_path, Filenames_Cat);
    for (int i = 0; i < params.getNbins(); i++)
    {
        BOOST_CHECK(
                boost::filesystem::exists(test_path / Filenames_Cat[i]));
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( split_test )
{
    std::vector<double> v = {1, 2, 3, 4, 5};
    std::vector<double> vsplit;
    unsigned int n = 2;
    vsplit = LE3_2D_MASS_WL_UTILITIES::split(v, n, 0);
    BOOST_CHECK(vsplit.size() == n + 1);
    BOOST_CHECK(vsplit[0] == v[0]);
    BOOST_CHECK(vsplit[1] == v[1]);
    BOOST_CHECK(vsplit[2] == v[2]);
    vsplit = LE3_2D_MASS_WL_UTILITIES::split(v, n, 1);
    BOOST_CHECK(vsplit.size() == n);
    BOOST_CHECK(vsplit[0] == v[3]);
    BOOST_CHECK(vsplit[1] == v[4]);
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
