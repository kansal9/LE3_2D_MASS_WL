/**
 * @file tests/src/ShearMap_test.cpp
 * @date 06/09/20
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

#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GenericMap.h"

#include "ElementsServices/DataSync.h"
#include "ElementsKernel/Temporary.h"
#include <fstream>
#include <string>

using boost::filesystem::path;

using namespace LE3_2D_MASS_WL_UTILITIES;
using namespace LE3_2D_MASS_WL_CARTESIAN;
using namespace ElementsServices::DataSync;

//----------------------------------------------------------------------------
struct ShearMapFixture
{
    ShearMapFixture() :
        myArrayTestMap(64, 64, 2), myArrayTestMapRef(64, 64, 2),
        ramin(0.), ramax(0.), decmin(0.), decmax(0.),
        mapUniformValueZ0(1.f), mapUniformValueZ1(2.f)
    {
        // Set values to the array
        for (size_t i = 0; i < myArrayTestMap.getXdim(); i++)
        {
            for (size_t j = 0; j < myArrayTestMap.getYdim(); j++)
            {
                for (size_t k = 0; k < myArrayTestMap.getZdim(); k++)
                {
                    myArrayTestMap.setBinValue(i, j, k,
                            mapUniformValueZ0 + k * mapUniformValueZ1);
                    myArrayTestMapRef.setBinValue(i, j, k,
                            mapUniformValueZ0 + k * mapUniformValueZ1);
                }
            }
        }
    }

    ShearMap myArrayTestMap;
    ShearMap myArrayTestMapRef;
    float ramin, ramax, decmin, decmax;
    float mapUniformValueZ0, mapUniformValueZ1;
};
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE (ShearMap_test)

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( arrayConstructor_test, ShearMapFixture )
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
                        myArrayTestMapRef.getBinValue(i, j, k), 0.000001);
            }
        }
    }
}

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( copyConstructor_test, ShearMapFixture)
{
    // Create a Map using the copy constructor
    ShearMap myCopyMap = ShearMap(myArrayTestMap);

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

    // Check that if one map changes, the other one is no modified (i.e. not sharing same smart pointer)
    myCopyMap.setBinValue(0, 0, 0, 42);

    // Check the values in the original map are still the same and correct
    for (size_t i = 0; i < myArrayTestMap.getXdim(); i++)
    {
        for (size_t j = 0; j < myArrayTestMap.getYdim(); j++)
        {
            for (size_t k = 0; k < myArrayTestMap.getZdim(); k++)
            {
                BOOST_CHECK_CLOSE(myArrayTestMap.getBinValue(i, j, k),
                        myArrayTestMapRef.getBinValue(i, j, k), 0.000001);
            }
        }
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
