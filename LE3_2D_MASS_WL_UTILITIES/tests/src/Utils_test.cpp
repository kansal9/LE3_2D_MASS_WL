/**
 * @file tests/src/Utils_test.cpp
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
#include <boost/filesystem/fstream.hpp>

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

#include "ElementsServices/DataSync.h"
#include "ElementsKernel/Temporary.h"

#include <ios>
#include <sstream>
#include <iostream>
#include <vector>

using namespace LE3_2D_MASS_WL_UTILITIES;
using namespace ElementsServices::DataSync;
using std::string;

// handle on created path names
boost::filesystem::path test_path;
//-----------------------------------------------------------------------------

struct UtilsDataSyncFixture
{
    DataSync sync;
    path filename_json, fitsFilePath_lensmc, filename_Sph, filename;
    UtilsDataSyncFixture() :
            sync(// Here is the connection configuration file
                 "LE3_2D_MASS_WL_UTILITIES/datasync_webdav.conf",
                 // Here is the dependency configuration file
                 "LE3_2D_MASS_WL_UTILITIES/test_file_list.txt"),
            filename_json(sync.absolutePath("filenames.json")),
            fitsFilePath_lensmc(sync.absolutePath("data/lensmcCat.fits")),
            filename_Sph(sync.absolutePath("data/Gamma_NSide16_Sphere.fits")),
            filename(sync.absolutePath("SampleInXML.xml"))
    {
        sync.download();
    }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (Utils_test, UtilsDataSyncFixture)
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getShearColNamesAndProxy_test )
{
    std::cout << "getShearColNamesAndProxy_test" << std::endl;
    auto v = getShearColNamesAndProxy("LENSMC");
    BOOST_CHECK(v.size() == 7);
    BOOST_CHECK(v[0].first == "SHE_LENSMC_UPDATED_RA");
    BOOST_CHECK(v[0].second == "ra");
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getClusterColNamesAndProxy_test )
{
    std::cout << "getShearColNamesAndProxy_test" << std::endl;
    auto v = getClusterColNamesAndProxy();
    BOOST_CHECK(v.size() == 6);
    BOOST_CHECK(v[0].first == "ID_DET_CLUSTER");
    BOOST_CHECK(v[0].second == "id");
}
//-----------------------------------------------------------------------------
// test generate time string
BOOST_AUTO_TEST_CASE(getDateTimeString_test)
{
    std::cout << "-- namespace LE3_2D_MASS_WL_UTILITIES: getDateTimeString_test"
            << std::endl;
    try
    {
        std::string time = getDateTimeString();
        BOOST_CHECK(true);
    } catch (...)
    {
        BOOST_CHECK(false);
    }
}
//-----------------------------------------------------------------------------
// test check file type
BOOST_AUTO_TEST_CASE(checkFileType_Fits_test)
{
    std::cout
            << "-- namespace LE3_2D_MASS_WL_UTILITIES: checkFileType_fits_test"
            << std::endl;
    try
    {
        if (checkFileType(fitsFilePath_lensmc, signFITS))
        {
            BOOST_CHECK(true);
        }
    } catch (...)
    {
        BOOST_CHECK(false);
    }
}
//-----------------------------------------------------------------------------
// test check file type
BOOST_AUTO_TEST_CASE(checkFileType_XML_test)
{
    std::cout << "-- namespace LE3_2D_MASS_WL_UTILITIES: checkFileType_XML_test"
            << std::endl;
    try
    {
        if (checkFileType(filename, signXML))
        {
            BOOST_CHECK(true);
        }
    } catch (...)
    {
        BOOST_CHECK(false);
    }
}
//-----------------------------------------------------------------------------
// test check file has field
BOOST_AUTO_TEST_CASE(checkFilehasfield_test)
{
    std::cout << "-- namespace LE3_2D_MASS_WL_UTILITIES: checkFilehasfield_test"
            << std::endl;
    try
    {
        if (fileHasField(filename.native(), "DpdTwoDMassLensMCCatalog"))
        {
            BOOST_CHECK(true);
        }
    } catch (...)
    {
        BOOST_CHECK(false);
    }
}
//-----------------------------------------------------------------------------
// test splits string
BOOST_AUTO_TEST_CASE(splitString_test)
{
    std::cout << "-- namespace LE3_2D_MASS_WL_UTILITIES: splits_string_test"
            << std::endl;
    try
    {
        const std::string str = "le3.wl.2dmass.output.patchconvergence";
        auto results = splitString(str, ".");
        if (results[2] == "2dmass")
        {
            BOOST_CHECK(true);
        }
    } catch (...)
    {
        BOOST_CHECK(false);
    }
}
//-----------------------------------------------------------------------------
// test read filenames from json files
BOOST_AUTO_TEST_CASE(check_read_filenames_test)
{
    std::cout
            << "-- namespace LE3_2D_MASS_WL_UTILITIES: check_read_filenames_test"
            << std::endl;
    try
    {
        fs::path workdir = filename_json.parent_path();
        fs::path _filename = filename_json.filename();
        std::vector<fs::path> filenames = readFilenamesInJson(
                workdir / _filename);
        if (false == filenames[0].string().empty())
        {
            BOOST_CHECK(true);
        }
    } catch (...)
    {
        BOOST_CHECK(false);
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( vecMinMax_test )
{
    std::cout << "-- Utils: VecMinMax" << std::endl;
    std::vector<double> vector = { 1, 5, 6, 7, 8, 30 };
    double min, max;
    try
    {
        vecMinMax(vector, min, max);
        BOOST_CHECK_EQUAL(min, 1);
        BOOST_CHECK_EQUAL(max, 30);



    } catch (...)
    {
        BOOST_FAIL("Exception in vecMinMax_test");
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getIndexCol_test )
{
    std::cout << "-- Utils: getIndexCol " << std::endl;
    std::vector<std::string> vector =
    { "RA", "DEC", "G1", "G2", "Z", "Weight" };
    int index;
    try
    {
        index = getIndexCol(vector, "G1");
        BOOST_CHECK_EQUAL(index, 2);
    } catch (...)
    {
        BOOST_FAIL("Exception in getIndexCol_test");
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getColName_test )
{
    std::cout << "-- Utils: getColName " << std::endl;
    std::vector<std::string> vector =
    { "KSB_RA", "KSB_DEC", "KSB_Gamma1", "KSB_Gamma2", "Z", "KSB_Weight" };
    std::string name;
    try
    {
        name = getColName(vector, "DEC");
        BOOST_CHECK(name == std::string("KSB_DEC"));
    } catch (...)
    {
        BOOST_FAIL("Exception in getColName_test");
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getColName2_test )
{
    std::cout << "-- Utils: getColName2 " << std::endl;
    std::vector<std::string> vector =
    { "RA", "DEC", "G1", "G2", "Z", "Mask" };
    std::string name;
    try
    {
        name = getColName(vector, "Mask", "Weight");
        BOOST_CHECK(name == std::string("Mask"));
    } catch (...)
    {
        BOOST_FAIL("Exception in getColName2_test");
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( readHealpixMap_test )
{
    std::cout << "-- Utils: readHealpixMap_Test" << std::endl;
    try
    {
        std::pair<Healpix_Map<double>, Healpix_Map<double> > mapPair;
        mapPair = readHealpixMap(filename_Sph.native());
        int m_npix = mapPair.first.Npix();
        int nside = mapPair.first.Nside();
        BOOST_CHECK(m_npix == (12 * nside * nside));
    } catch (...)
    {
        BOOST_CHECK(false);
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getIndex2RaDec_test )
{
    std::cout << "-- Utils: getIndex2RaDec_Test" << std::endl;
    try
    {
        std::pair<Healpix_Map<double>, Healpix_Map<double> > mapPair;
        mapPair = readHealpixMap(filename_Sph.native());
        std::pair<double, double> radec = getIndex2RaDec(mapPair.first, 5);
        BOOST_CHECK_CLOSE(radec.first, 67.5, 0.1);
    } catch (...)
    {
        BOOST_CHECK(false);
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getJsonOutput_test )
{
    std::cout << "-- Utils: getJsonOutput_Test" << std::endl;
    try
    {
        using Elements::TempDir;
        TempDir one;
        test_path = one.path();
        fs::path outFilename = boost::filesystem::path("outfile.json");
        std::vector<fs::path> filenames =
        { "abc.fits", "bcd.fits", "cde.fits", "def.fits" };
        fillJsonOutput(test_path, outFilename, filenames);
        BOOST_CHECK(fs::exists(test_path / outFilename));
    } catch (...)
    {
        BOOST_CHECK(false);
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( fillVecColumnLinspace_test )
{
    std::cout << "-- Utils: fillVecColumnLinspace_test" << std::endl;
    size_t N = 10;
    double start = 0;
    double stop = 10;
    double step;
    VecColumn<double> col({"", ""}, N);

    fillVecColumnLinspace(col, N, start, stop);
    step = (stop - start) / (N - 1);
    for(size_t i = 0; i < N; i++)
    {
        BOOST_CHECK_EQUAL(col[i], start + i * step);
    }

    fillVecColumnLinspace(col, N, start, stop, false);
    step = (stop - start) / N;
    for(size_t i = 0; i < N; i++)
    {
        BOOST_CHECK_EQUAL(col[i], start + i * step);
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getXmlProductType_test )
{
    std::cout << "-- Utils: getXmlProductType_test" << std::endl;
    std::string product_type = getXmlProductType(filename);
    BOOST_CHECK(product_type == "DpdTwoDMassLensMCCatalog");
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
