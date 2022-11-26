/**
 * @file tests/src/Matrix_test.cpp
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
#include "ElementsKernel/Logging.h"
#include "LE3_2D_MASS_WL_UTILITIES/Matrix.h"
#include <iostream>

using LE3_2D_MASS_WL_UTILITIES::Matrix;
static Elements::Logging logger = Elements::Logging::getLogger("Matrix_test");
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE (Matrix_test)
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( constructor_tests )
{
    logger.info() << "-- Matrix: Constructor_test";
    // Create an image initialized with zeros
    unsigned int imSize(32);
    Matrix myImage(imSize, imSize);
    myImage.clear();
    // Check all values are well zeros
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            BOOST_CHECK_CLOSE(myImage(i, j), 0, 1);
        }
    }
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( values_tests )
{
    logger.info() << "-- Matrix: values_tests";
    unsigned int imSize(32);

    // Create an image initialized with values
    Matrix myImage(imSize, imSize);
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            myImage(i, j) = double(i - j);
        }
    }

    // Check values are well initialized
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            BOOST_CHECK_CLOSE(myImage(i, j), i - j, 0.01);
        }
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( copyConstructor_tests )
{
    logger.info() << "-- Matrix: copyConstructor_test";
    unsigned int imSize(32);

    // Create an image initialized with values
    Matrix myImage(imSize, imSize);
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            myImage(i, j) = double(i - j);
        }
    }

    // Create another image from this one
    Matrix myImage2(myImage);

    // Check values are well the same
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            BOOST_CHECK_CLOSE(myImage(i, j), myImage2(i, j), 0.01);
        }
    }
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( operatorEqual_tests )
{
    logger.info() << "-- Matrix: operatorEqual_test";
    unsigned int imSize(32);

    // Create an image initialized with values
    Matrix myImage(imSize, imSize);
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            myImage(i, j) = double(i - j);
        }
    }

    // Create another image from this one
    Matrix myImage2(myImage);

    // Check both images do not share same memory
    myImage(0, 0) = 10;
    BOOST_CHECK_CLOSE(myImage(0, 0) - 10, myImage2(0, 0), 0.01);

    // Check the operator= now
    myImage2 = myImage;
    // Check values are well the same
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            BOOST_CHECK_CLOSE(myImage(i, j), myImage2(i, j), 0.01);
        }
    }
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getXY_tests )
{
    logger.info() << "-- Matrix: getXY_test";
    // Create an image initialized with zeros
    unsigned int imSize(32);
    Matrix myImage(imSize, imSize * 2);
    // Check the dimensions X and Y are correct
    BOOST_CHECK_EQUAL(myImage.getXdim(), imSize);
    BOOST_CHECK_EQUAL(myImage.getYdim(), imSize * 2);
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( add_substract_multiply_tests )
{
    logger.info() << "-- Matrix: add_substract_multiply_test";
    unsigned int imSize(32);

    // Create two images
    Matrix myImage(imSize, imSize);
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            myImage(i, j) = double(i - j);
        }
    }

    Matrix myImage2(myImage);

    // Get the added image
    Matrix myAddedImage = myImage + myImage2;
    // Get the substracted image
    Matrix mySubImage = myImage - myImage2;
    // Get the multiplied image
    double factor = 2.58;
    Matrix myMulImage = myImage * factor;

    // Check the values
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            BOOST_CHECK_CLOSE(myAddedImage(i, j), 2. * double(i - j), 0.01);
            BOOST_CHECK_CLOSE(mySubImage(i, j), 0., 0.01);
            BOOST_CHECK_CLOSE(myMulImage(i, j), factor * double(i - j), 0.01);
        }
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getMax_tests )
{
    logger.info() << "-- Matrix: getMax_test";
    unsigned int imSize(32);

    // Create two images
    Matrix myImage(imSize, imSize);
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            myImage(i, j) = double(i + j);
        }
    }

    // Get the max
    double max = myImage.getMax();

    // Check the max is the right value
    BOOST_CHECK_CLOSE(max, 2 * (imSize - 1), 1);
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getMin_tests )
{
    logger.info() << "-- Matrix: getMin_test";
    unsigned int imSize(32);

    Matrix myImage(imSize, imSize);
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            myImage(i, j) = double(i + j);
        }
    }

    // Get the max
    double Min = myImage.getMin();

    // Check the max is the right value
    BOOST_CHECK_CLOSE(Min, 0, 1);
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getFlux_tests )
{
    logger.info() << "-- Matrix: getFlux_test";
    unsigned int imSize(32);

    Matrix myImage(imSize, imSize);
    double tot = 0;
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            double val = i * 1.5;
            myImage(i, j) = val;
            tot += val;
        }
    }

    double flux = myImage.getFlux();
    logger.info() << "-- Matrix: getFlux_test" << flux;
    BOOST_CHECK_CLOSE(flux, tot, 1);
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( reset_tests )
{
    logger.info() << "-- Matrix: reset_test";
    unsigned int imSize(32);

    Matrix myImage(imSize, imSize);
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            myImage(i, j) = double(i + j);
        }
    }

    // reset image
    myImage.clear();

    // Check the values
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            BOOST_CHECK_CLOSE(myImage(i, j), 0., 0.01);
        }
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getStandardDeviation_tests )
{
    logger.info() << "-- Matrix: getStandardDeviation_test";
    unsigned int imSize(32);

    Matrix myImage(imSize, imSize);
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            myImage(i, j) = double(i + j);
        }
    }

    // Get the standard deviation
    double stdev = myImage.getSigma();

    // Check the max is the right value
    BOOST_CHECK_CLOSE(stdev, 13.06394529, 1);
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( applyThreshold_tests )
{
    logger.info() << "-- Matrix: applyThreshold_test";
    unsigned int imSize(32);

    Matrix myImage(imSize, imSize);
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            myImage(i, j) = double(i + j);
        }
    }

    // Apply a threshold
    myImage.applyThreshold(10.);

    // Check the values below the threshold are well set to zero
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            if (i + j < 10.)
            {
                BOOST_CHECK_CLOSE(myImage(i, j), 0, 1);
            }
            else
            {
                BOOST_CHECK_CLOSE(myImage(i, j), i + j, 1);
            }
        }
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( localMax_tests )
{
    logger.info() << "-- Matrix: localMax_test";
    unsigned int imSize(32);

    Matrix myImage(imSize, imSize);
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            myImage(i, j) = double(i + j);
        }
    }

    // Check there is no local max in this context
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            BOOST_CHECK(myImage.isLocalMax(i, j) == false);
        }
    }
    // Now add arbitrary local maximums
    myImage(18, 12) = 70;
    myImage(3, 17) = 85;

    // Check those pixels are now local maximum
    BOOST_CHECK(myImage.isLocalMax(18, 12) == true);
    BOOST_CHECK(myImage.isLocalMax(3, 17) == true);
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getandsetValues_tests )
{
    logger.info() << "-- Matrix: getandsetValues_test";
    // Create an image initialized with zeros
    unsigned int imSize(32);
    Matrix myImage(imSize, imSize);

    // Set the values
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            myImage(i, j) = 2 * i + j;
        }
    }
    // Check the values are well set
    for (unsigned int i = 0; i < imSize; i++)
    {
        for (unsigned int j = 0; j < imSize; j++)
        {
            BOOST_CHECK_CLOSE(myImage(i, j), 2 * i + j, 0.01);
        }
    }
    // Check the values are well set
    for (unsigned int i = 0; i < imSize; i++)
    {
        BOOST_CHECK_CLOSE(myImage(i, imSize - 1), 2 * i + (imSize - 1), 0.01);
    }
    for (unsigned int j = 0; j < imSize; j++)
    {
        BOOST_CHECK_CLOSE(myImage(imSize - 1, j), 2 * (imSize - 1) + j, 0.01);
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
