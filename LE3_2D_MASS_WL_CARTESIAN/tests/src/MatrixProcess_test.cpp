/**
 * @file tests/src/MatrixProcess_test.cpp
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

#include "LE3_2D_MASS_WL_CARTESIAN/MatrixProcess.h"
#include "LE3_2D_MASS_WL_UTILITIES/Matrix.h"

#include <iostream>

using namespace LE3_2D_MASS_WL_CARTESIAN;
using LE3_2D_MASS_WL_UTILITIES::Matrix;

static Elements::Logging logger = Elements::Logging::getLogger(
        "MatrixProcess_test");
struct MatrixProcessTestEnv
{
    static const unsigned int m_imSize = 32;
    Matrix m_values;
    Matrix m_valuesRef;

    MatrixProcessTestEnv() :
            m_values(m_imSize, m_imSize)
    {
        // Define arbitraty values for the image
        for (unsigned int i = 0; i < m_imSize; i++)
        {
            for (unsigned int j = 0; j < m_imSize; j++)
            {
                m_values(i, j) = double(3 + i + 2 * j);
            }
        }

        m_valuesRef = m_values;
    }
};
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (MatrixProcess_test, MatrixProcessTestEnv)
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( DCT_test )
{
    logger.info() << "-- MatrixProcess: DCT_test";

    // Create the image
    Matrix myDCTImage(m_imSize, m_imSize);

    // Create the IP object
    MatrixProcess myMP(m_imSize, m_imSize);

    // Perform DCT on the image
    myMP.performDCT(m_values.data().begin(), myDCTImage.data().begin());

    // And then perform IDCT on the output DCT image
    myMP.performIDCT(myDCTImage.data().begin(), m_values.data().begin());

    // Check the values of the original and back image are close
    for (unsigned int i = 0; i < m_imSize; i++)
    {
        for (unsigned int j = 0; j < m_imSize; j++)
        {
            BOOST_CHECK_CLOSE(m_values(i, j), m_valuesRef(i, j), 0.1);
        }
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( spline_test )
{
    logger.info() << "-- MatrixProcess: spline_test";

    // Create the MP object
    MatrixProcess myMP(m_imSize, m_imSize);

    // Perform an arbitrary wavelet decomposition
    int nScales = 5;
    std::vector<Matrix> myBand(nScales);
    myMP.transformBspline(m_values, myBand, nScales);

    // Reconstruct the image
    myMP.reconsBspline(myBand, m_values);

    // Check the values of the original and back image are close
    for (unsigned int i = 0; i < m_imSize; i++)
    {
        for (unsigned int j = 0; j < m_imSize; j++)
        {
            printf("(%d,%d): %f %f\n", i, j, m_values(i, j), m_valuesRef(i, j));
            BOOST_CHECK_CLOSE(m_values(i, j), m_valuesRef(i, j), 0.1);
        }
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
