/**
 * @file tests/src/GetMap_test.cpp
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

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "ElementsKernel/Logging.h"
#include "ElementsServices/DataSync.h"
#include "ElementsKernel/Temporary.h"

#include "LE3_2D_MASS_WL_CARTESIAN/GenericMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"

#include <fstream>
#include <string>

using boost::filesystem::path;
using namespace LE3_2D_MASS_WL_CARTESIAN;
using namespace ElementsServices::DataSync;

static Elements::Logging logger = Elements::Logging::getLogger("GetMap_test");
//----------------------------------------------------------------------------
struct GetMapFixture
{
    DataSync sync;
    path testMapFilename;
    GetMapFixture() :
        sync("LE3_2D_MASS_WL_CARTESIAN/datasync_webdav.conf",
             "LE3_2D_MASS_WL_CARTESIAN/test_file_list.txt"),
        testMapFilename(sync.absolutePath("data/convMap_test.fits")),
        myArrayTestMap(16, 16, 2, 100000000),
        myArrayTestMapNoGal(16, 16, 2),
        ramin(0.), ramax(0.), decmin(0.), decmax(0.),
        mapUniformValueZ0(1.f), mapUniformValueZ1(2.f)
    {
        sync.download();

        // Set values to the array
        for (size_t i = 0; i < myArrayTestMap.getXdim(); i++)
        {
            for (size_t j = 0; j < myArrayTestMap.getYdim(); j++)
            {
                for (size_t k = 0; k < myArrayTestMap.getZdim(); k++)
                {
                    myArrayTestMap.setBinValue(i, j, k,
                            mapUniformValueZ0 + k * mapUniformValueZ1);
                    myArrayTestMapNoGal.setBinValue(i, j, k,
                            mapUniformValueZ0 + k * mapUniformValueZ1);
                }
            }
        }
    }

    GenericMap myArrayTestMap;
    GenericMap myArrayTestMapNoGal;
    float ramin, ramax, decmin, decmax;
    float mapUniformValueZ0, mapUniformValueZ1;
};

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE (GetMap_test)

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( arrayConstructor_test, GetMapFixture)
{
    // Check the values in the map are corresponding to the values in the input array
    for (size_t i = 0; i < myArrayTestMap.getXdim(); i++)
    {
        for (size_t j = 0; j < myArrayTestMap.getYdim(); j++)
        {
            for (size_t k = 0; k < myArrayTestMap.getZdim(); k++)
            {
                // Check the difference between the map value and the input array is lower than 0.0001%
                BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, k),
                        mapUniformValueZ0 + k * mapUniformValueZ1, 0.000001);
            }
        }
    }
}

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( copyConstructor_test, GetMapFixture)
{
    // Create a Map using the copy constructor
    GenericMap myCopyMap = GenericMap(myArrayTestMap);

    // Check this copied map has the same properties as the original map
    BOOST_CHECK(myCopyMap.getXdim() == myArrayTestMap.getXdim());
    BOOST_CHECK(myCopyMap.getYdim() == myArrayTestMap.getYdim());
    BOOST_CHECK(myCopyMap.getZdim() == myArrayTestMap.getZdim());

    // Check the values in the map are corresponding to the values of the original map
    for (size_t i = 0; i < myArrayTestMap.getXdim(); i++)
    {
        for (size_t j = 0; j < myArrayTestMap.getYdim(); j++)
        {
            for (size_t k = 0; k < myArrayTestMap.getZdim(); k++)
            {
                // Check the difference between the map value and the input array is lower than 0.0001%
                BOOST_CHECK_CLOSE(myCopyMap.getBinValue(i, j, k),
                        myArrayTestMap.getBinValue(i, j, k), 0.000001);
            }
        }
    }

    // Check the values in the original map are still available and correct
    for (size_t i = 0; i < myArrayTestMap.getXdim(); i++)
    {
        for (size_t j = 0; j < myArrayTestMap.getYdim(); j++)
        {
            for (size_t k = 0; k < myArrayTestMap.getZdim(); k++)
            {
                BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, k),
                        mapUniformValueZ0 + k * mapUniformValueZ1, 0.000001);
            }
        }
    }
}

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( arrayConstructorFromFile_test, GetMapFixture)
{
    GenericMap map(testMapFilename);
    // check that matrix is load from file with its dimensions
    BOOST_CHECK(map.getXdim() == 128);
    BOOST_CHECK(map.getYdim() == 128);
    BOOST_CHECK(map.getZdim() == 3);

    // try with wrong filename
    map = GenericMap((std::string)"wrong_filename.fits");
    BOOST_CHECK(map.getXdim() == 0);
    BOOST_CHECK(map.getYdim() == 0);
    BOOST_CHECK(map.getZdim() == 0);
}

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( getBinValue_test, GetMapFixture)
{
    // Check the values outside of map return well the edge of the map
    BOOST_CHECK_CLOSE(
            myArrayTestMap.getBinValue(myArrayTestMap.getXdim(), 0, 0),
            myArrayTestMap.getBinValue(myArrayTestMap.getXdim() - 1, 0, 0),
            0.001);
    BOOST_CHECK_CLOSE(
            myArrayTestMap.getBinValue(0, myArrayTestMap.getYdim(), 0),
            myArrayTestMap.getBinValue(0, myArrayTestMap.getYdim() - 1, 0),
            0.001);
    BOOST_CHECK_CLOSE(
            myArrayTestMap.getBinValue(0, 0, myArrayTestMap.getZdim()),
            myArrayTestMap.getBinValue(0, 0, myArrayTestMap.getZdim() - 1),
            0.001);
    BOOST_CHECK_CLOSE(
            myArrayTestMap.getBinValue(myArrayTestMap.getXdim(),
                    myArrayTestMap.getYdim(), myArrayTestMap.getZdim()),
            myArrayTestMap.getBinValue(myArrayTestMap.getXdim() - 1,
                    myArrayTestMap.getYdim() - 1, myArrayTestMap.getZdim() - 1),
            0.001);

    // Check the values in the map are corresponding to the values in the input array
    for (size_t i = 0; i < myArrayTestMap.getXdim(); i++)
    {
        for (size_t j = 0; j < myArrayTestMap.getYdim(); j++)
        {
            for (size_t k = 0; k < myArrayTestMap.getZdim(); k++)
            {
                // Check the difference between the map value and the input array is lower than 0.0001%
                BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, k),
                        mapUniformValueZ0 + k * mapUniformValueZ1, 0.000001);
            }
        }
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( getAxisDim_test, GetMapFixture)
{
    // Check the values of the dimensions are the same
    BOOST_CHECK(myArrayTestMap.getXdim() == myArrayTestMapNoGal.getXdim());
    BOOST_CHECK(myArrayTestMap.getYdim() == myArrayTestMapNoGal.getYdim());
    BOOST_CHECK(myArrayTestMap.getZdim() == myArrayTestMapNoGal.getZdim());

    // Try to rebin the map
    bool isPixelated = myArrayTestMap.pixelate(1, 1);

    // Check pixelation worked well
    BOOST_CHECK(isPixelated == true);

    // Check the values of the dimensions are well updated
    BOOST_CHECK(myArrayTestMap.getXdim() == myArrayTestMapNoGal.getXdim() / 2);
    BOOST_CHECK(myArrayTestMap.getYdim() == myArrayTestMapNoGal.getYdim() / 2);
    BOOST_CHECK(myArrayTestMap.getZdim() == myArrayTestMapNoGal.getZdim());
}

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( getMeanValues_test, GetMapFixture)
{
    // Perform the check for the map with galaxy info
    std::vector<double> myMeans = myArrayTestMap.getMeanValues();
    for (size_t k = 0; k < myArrayTestMap.getZdim(); k++)
    {
        BOOST_CHECK_CLOSE(myMeans[k], mapUniformValueZ0 + k*mapUniformValueZ1, 0.001);
    }
}

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( removeOffset_test, GetMapFixture)
{
    // define offsets to remove then to the map
    std::vector<double> offsets;
    offsets.push_back(0.3 * mapUniformValueZ0);
    offsets.push_back(0.45 * mapUniformValueZ1);

    // Remove those offsets to the map
    myArrayTestMap.removeOffset(offsets);
    // Check that the values of the map are in agreement
    for (size_t i = 0; i < myArrayTestMap.getXdim(); i++)
    {
        for (size_t j = 0; j < myArrayTestMap.getYdim(); j++)
        {
            BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, 0),
                    0.7 * mapUniformValueZ0, 0.001);
            BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, 1),
                    mapUniformValueZ0 + 0.55 * mapUniformValueZ1, 0.001);
        }
    }
}

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( pixelatePerBin_test, GetMapFixture)
{
    // Try to rebin with same binning
    bool isPixelated = myArrayTestMap.pixelate(0, 0);
    // Should return false
    BOOST_CHECK(isPixelated == false);

    // Check the values of the dimensions are unchanged
    BOOST_CHECK(myArrayTestMap.getXdim() == myArrayTestMapNoGal.getXdim());
    BOOST_CHECK(myArrayTestMap.getYdim() == myArrayTestMapNoGal.getYdim());
    BOOST_CHECK(myArrayTestMap.getZdim() == myArrayTestMapNoGal.getZdim());

    // Try to rebin so that the number of bin would be less or equal to 1
    isPixelated = myArrayTestMap.pixelate(8, 8);
    BOOST_CHECK(isPixelated == false);

    isPixelated = myArrayTestMap.pixelate(0, 8);
    BOOST_CHECK(isPixelated == false);

    isPixelated = myArrayTestMap.pixelate(8, 0);
    BOOST_CHECK(isPixelated == false);

    // Try to rebin only X axis
    isPixelated = myArrayTestMap.pixelate(1, 0);
    // Should return true
    BOOST_CHECK(isPixelated == true);

    // Check the values of the dimensions are well updated
    BOOST_CHECK(myArrayTestMap.getXdim() == myArrayTestMapNoGal.getXdim() / 2);
    BOOST_CHECK(myArrayTestMap.getYdim() == myArrayTestMapNoGal.getYdim());
    BOOST_CHECK(myArrayTestMap.getZdim() == myArrayTestMapNoGal.getZdim());

    // Check the values of the map are well updated
    // Set values to the array
    for (size_t i = 0; i < myArrayTestMap.getXdim(); i++)
    {
        for (size_t j = 0; j < myArrayTestMap.getYdim(); j++)
        {
            for (size_t k = 0; k < myArrayTestMap.getZdim(); k++)
            {
                BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, k),
                        myArrayTestMapNoGal.getBinValue(i, j, k), 0.001);
            }
        }
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( pixelatePerBin_test2, GetMapFixture)
{
    // Try to rebin only Y axis
    bool isPixelated = myArrayTestMap.pixelate(0, 1);
    // Should return true
    BOOST_CHECK(isPixelated == true);

    // Check the values of the dimensions are well updated
    BOOST_CHECK(myArrayTestMap.getXdim() == myArrayTestMapNoGal.getXdim());
    BOOST_CHECK(myArrayTestMap.getYdim() == myArrayTestMapNoGal.getYdim() / 2);
    BOOST_CHECK(myArrayTestMap.getZdim() == myArrayTestMapNoGal.getZdim());

    // Check the values of the map are well updated
    // Set values to the array
    for (size_t i = 0; i < myArrayTestMap.getXdim(); i++)
    {
        for (size_t j = 0; j < myArrayTestMap.getYdim(); j++)
        {
            for (size_t k = 0; k < myArrayTestMap.getZdim(); k++)
            {
                BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, k),
                        mapUniformValueZ0 + k * mapUniformValueZ1, 0.001);
            }
        }
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( pixelatePerBin_test3, GetMapFixture)
{
    // Try to rebin both X and Y axis
    bool isPixelated = myArrayTestMap.pixelate(2, 2);
    // Should return true
    BOOST_CHECK(isPixelated == true);

    // Check the values of the dimensions are well updated
    BOOST_CHECK(myArrayTestMap.getXdim() == myArrayTestMapNoGal.getXdim() / 4);
    BOOST_CHECK(myArrayTestMap.getYdim() == myArrayTestMapNoGal.getYdim() / 4);
    BOOST_CHECK(myArrayTestMap.getZdim() == myArrayTestMapNoGal.getZdim());

    // Check the values of the map are well updated
    // Set values to the array
    for (size_t i = 0; i < myArrayTestMap.getXdim(); i++)
    {
        for (size_t j = 0; j < myArrayTestMap.getYdim(); j++)
        {
            for (size_t k = 0; k < myArrayTestMap.getZdim(); k++)
            {
                BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, k),
                        mapUniformValueZ0 + k * mapUniformValueZ1, 0.001);
            }
        }
    }
}

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( addBorders_test, GetMapFixture)
{
    std::cout << "-- GetMap: addBorders_test" << std::endl;
    // Add borders
    myArrayTestMap.add_borders();

    // Check the size of the map is multiplied by 2 in X and Y
    BOOST_CHECK(myArrayTestMap.getXdim() == 2 * myArrayTestMapNoGal.getXdim());
    BOOST_CHECK(myArrayTestMap.getYdim() == 2 * myArrayTestMapNoGal.getYdim());
    BOOST_CHECK(myArrayTestMap.getZdim() == myArrayTestMapNoGal.getZdim());

    // Check the borders added are zeros
    for (size_t i = 0; i < myArrayTestMap.getXdim(); i++)
    {
        for (size_t j = 0; j < myArrayTestMap.getYdim(); j++)
        {
            if (i < myArrayTestMap.getXdim() / 4
                    || i >= 3 * myArrayTestMap.getXdim() / 4
                    || j < myArrayTestMap.getYdim() / 4
                    || j >= 3 * myArrayTestMap.getYdim() / 4)
            {
                BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, 0), 0, 0.1);
                BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, 1), 0, 0.1);
            }
            else
            {
                BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, 0),
                        myArrayTestMapNoGal.getBinValue(
                                i - myArrayTestMap.getXdim() / 4,
                                j - myArrayTestMap.getYdim() / 4, 0), 0.1);
                BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, 1),
                        myArrayTestMapNoGal.getBinValue(
                                i - myArrayTestMap.getXdim() / 4,
                                j - myArrayTestMap.getYdim() / 4, 1), 0.1);
            }
        }
    }

    // Remove borders
    myArrayTestMap.remove_borders();

    // Check the size of the map is set back to original values
    BOOST_CHECK(myArrayTestMap.getXdim() == myArrayTestMapNoGal.getXdim());
    BOOST_CHECK(myArrayTestMap.getYdim() == myArrayTestMapNoGal.getYdim());
    BOOST_CHECK(myArrayTestMap.getZdim() == myArrayTestMapNoGal.getZdim());

    // Check the values are the original ones in the map
    for (size_t i = 0; i < myArrayTestMap.getXdim(); i++)
    {
        for (size_t j = 0; j < myArrayTestMap.getYdim(); j++)
        {
            for (size_t k = 0; k < myArrayTestMap.getZdim(); k++)
            {
                // Check the difference between the map value and the input array is lower than 0.0001%
                BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, k),
                        myArrayTestMapNoGal.getBinValue(i, j, k), 0.000001);
            }
        }
    }
}

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( thresholding_test, GetMapFixture)
{
    std::cout << "-- GetMap: thresholding_test" << std::endl;

    double th = 0.5;
    int k = 0;

    // Apply a threshold
    myArrayTestMap.applyThreshold(th);

    // Check the values below the threshold are well set to zero
    double val = 0;
    for (size_t i = 0; i < myArrayTestMap.getXdim(); i++)
    {
        for (size_t j = 0; j < myArrayTestMap.getYdim(); j++)
        {
            val = myArrayTestMapNoGal.getBinValue(i, j, k);
            if (val < th)
            {
                BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, k), 0, 1);
            }
            else
            {
                BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, k), val, 1);
            }

        }
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( getSigma_test, GetMapFixture)
{
    // Get the standard deviation
    double stdev = myArrayTestMap.getSigma(0);

    // Check the max is the right value
    BOOST_CHECK_CLOSE(stdev, 0, 1);
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( getMax_test, GetMapFixture)
{
    for (size_t k = 0; k < myArrayTestMap.getZdim(); k++)
    {
        myArrayTestMap.setBinValue(0, 0, k, 100);
        double max = myArrayTestMap.getMax(k);
        BOOST_CHECK_CLOSE(max, myArrayTestMap.getBinValue(0, 0, k), 0.000001);
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( getMin_test, GetMapFixture)
{
    for (size_t k = 0; k < myArrayTestMap.getZdim(); k++)
    {
        myArrayTestMap.setBinValue(0, 0, k, -100);
        double min = myArrayTestMap.getMin(k);
        BOOST_CHECK_CLOSE(min, myArrayTestMap.getBinValue(0, 0, k), 0.000001);
    }
}
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( getFlux_test, GetMapFixture)
{
    for (size_t k = 0; k < myArrayTestMap.getZdim(); k++)
    {
        double flux = myArrayTestMap.getFlux(k);
        double expected_flux = (mapUniformValueZ0 + k * mapUniformValueZ1)
                * myArrayTestMap.getXdim() * myArrayTestMap.getYdim();
        BOOST_CHECK_CLOSE(flux, expected_flux, 0.000001);
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
