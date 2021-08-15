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
#include "LE3_2D_MASS_WL_CARTESIAN/Matrix.h"

#include <iostream>

using namespace LE3_2D_MASS_WL_CARTESIAN;
static Elements::Logging logger = Elements::Logging::getLogger("MatrixProcess_test");
struct MatrixProcessTestEnv {
  unsigned int imSize = 32;
  double *values = new double[imSize*imSize];
  MatrixProcessTestEnv  () {
  // Define arbitraty values for the image
  for (unsigned int i=0; i<imSize; i++)
  {
    for (unsigned int j=0; j<imSize; j++)
    {
      values[j*imSize + i] = double(3+i+2*j);
    }
  }
  }
  // delete all test products
  ~MatrixProcessTestEnv  () {};
};
//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_SUITE (MatrixProcess_test)
BOOST_FIXTURE_TEST_SUITE (MatrixProcess_test, MatrixProcessTestEnv)
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( DCT_test ) {
  logger.info() << "-- MatrixProcess: DCT_test";

  // Create the image
  Matrix myImage(imSize, imSize, values);

  // Create the IP object
  MatrixProcess myMP(imSize, imSize);

  // Perform DCT on the image
  Matrix myDCTImage = myMP.performDCT(myImage);

  // And then perform IDCT on the output DCT image
  Matrix myImageBack = myMP.performIDCT(myDCTImage);

  // Check the values of the original and back image are close
  for (unsigned int i=0; i<imSize; i++)
  {
    for (unsigned int j=0; j<imSize; j++)
    {
      BOOST_CHECK_CLOSE(myImage.getValue(i, j), myImageBack.getValue(i, j), 0.1);
    }
  }

  delete [] values;
  values = nullptr;
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( DCTblocks_test ) {
  logger.info() << "-- MatrixProcess: DCTblocks_test";

  // Create the image
  Matrix myImage(imSize, imSize, values);

  // Create the IP object
  MatrixProcess myMP(imSize, imSize);

  // Perform DCT on the image
  unsigned int blockSize(8);
  Matrix myDCTImage = myMP.performDCT(myImage, blockSize, blockSize, true);

  // And then perform IDCT on the output DCT image
  Matrix myImageBack = myMP.performDCT(myDCTImage, blockSize, blockSize, false);

  // Check the values of the original and back image are close
  for (unsigned int i=0; i<imSize; i++)
  {
    for (unsigned int j=0; j<imSize; j++)
    {
      BOOST_CHECK_CLOSE(myImage.getValue(i, j), myImageBack.getValue(i, j), 0.1);
    }
  }

  delete [] values;
  values = nullptr;
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( spline_test ) {
  logger.info() << "-- MatrixProcess: spline_test";
  // Create the image
  Matrix myImage(imSize, imSize, values);

  // Create the IP object
  MatrixProcess myMP(imSize, imSize);

  // Perform an arbitrary wavelet decomposition
  std::vector<Matrix> myBand = myMP.transformBspline(myImage, 5);

  // Reconstruct the image
  Matrix myImageBack = myMP.reconsBspline(myBand);

  // Check the values of the original and back image are close
  for (unsigned int i=0; i<imSize; i++)
  {
    for (unsigned int j=0; j<imSize; j++)
    {
      BOOST_CHECK_CLOSE(myImage.getValue(i, j), myImageBack.getValue(i, j), 0.1);
    }
  }

  delete [] values;
  values = nullptr;
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
