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
#include "LE3_2D_MASS_WL_CARTESIAN/Matrix.h"
#include <iostream>

using LE3_2D_MASS_WL_CARTESIAN::Matrix;
static Elements::Logging logger = Elements::Logging::getLogger("Matrix_test");
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE (Matrix_test)

//-----------------------------------------------------------------------------
// Constructor Test with values
BOOST_AUTO_TEST_CASE( Constructor_tests ) {
  logger.info() << "-- Matrix: Constructor_test";
 // Create an image initialized with zeros
 unsigned int imSize(32);
 Matrix myImage(imSize, imSize);
 // Check all values are well zeros
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   BOOST_CHECK_CLOSE(myImage.getValue(i, j), 0, 1);
  }
  }
 }

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( values_tests ) {
  logger.info() << "-- Matrix: values_tests";
 unsigned int imSize(32);
 double *values = new double[imSize*imSize];
 // Define arbitraty values for an image
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   values[j*imSize + i] = double(i-j);
  }
 }
 // Create an image initialized with values
 Matrix myImage(imSize, imSize, values);
 // Check values are well initialized
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   BOOST_CHECK_CLOSE(myImage.getValue(i, j), i-j, 0.01);
  }
 }
 delete [] values;
 values = nullptr;
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( copyConstructor_tests ) {
  logger.info() << "-- Matrix: copyConstructor_test";
 unsigned int imSize(32);
 double *values = new double[imSize*imSize];
 // Define arbitraty values for an image
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   values[j*imSize + i] = double(i-j);
  }
 }
 // Create an image initialized with values
 Matrix myImage(imSize, imSize, values);
 // Create another image from this one
 Matrix myImage2(myImage);
 // Check values are well the same
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   BOOST_CHECK_CLOSE(myImage.getValue(i, j), myImage2.getValue(i, j), 0.01);
  }
 }
 delete [] values;
 values = nullptr;
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( operatorEqual_tests ) {
  logger.info() << "-- Matrix: operatorEqual_test";
 unsigned int imSize(32);
 double *values = new double[imSize*imSize];
 // Define arbitraty values for an image
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++){
   values[j*imSize + i] = double(i-j);
  }
 }
 // Create an image initialized with values
 Matrix myImage(imSize, imSize, values);
 // Create another image from this one
 Matrix myImage2(myImage);
 // Check both images do not share same memory
 myImage.setValue(0, 0, 10);
 BOOST_CHECK_CLOSE(myImage.getValue(0, 0)-10, myImage2.getValue(0, 0), 0.01);

 // Check the operator= now
 myImage2 = myImage;
 // Check values are well the same
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   BOOST_CHECK_CLOSE(myImage.getValue(i, j), myImage2.getValue(i, j), 0.01);
  }
 }
 delete [] values;
 values = nullptr;
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getXY_tests ) {
  logger.info() << "-- Matrix: getXY_test";
 // Create an image initialized with zeros
 unsigned int imSize(32);
 Matrix myImage(imSize, imSize*2);
 // Check the dimensions X and Y are correct
 BOOST_CHECK_EQUAL(myImage.getXdim(), imSize);
 BOOST_CHECK_EQUAL(myImage.getYdim(), imSize*2);
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( add_substract_multiply_tests ) {
  logger.info() << "-- Matrix: add_substract_multiply_test";
 unsigned int imSize(32);
 double *values = new double[imSize*imSize];
 // Define arbitraty values for an image
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   values[j*imSize + i] = double(i-j);
  }
 }
 // Create two images
 Matrix myImage(imSize, imSize, values);
 Matrix myImage2(imSize, imSize, values);

 // Get the added image
 Matrix myAddedImage = myImage.add(myImage2);
 // Get the substracted image
 Matrix mySubImage = myImage.substract(myImage2);
 // Get the multiplied image
 double factor = 2.58;
 Matrix myMulImage = myImage.multiply(factor);

 // Check the values
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   BOOST_CHECK_CLOSE(myAddedImage.getValue(i, j), 2.*double(i-j), 0.01);
   BOOST_CHECK_CLOSE(mySubImage.getValue(i, j), 0., 0.01);
   BOOST_CHECK_CLOSE(myMulImage.getValue(i, j), factor*double(i-j), 0.01);
  }
 }
 delete [] values;
 values = nullptr;
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( AddSubstract_differentSizeImage_tests ) {
  logger.info() << "-- Matrix: Add_Substract_different_Size_Image_test";
 unsigned int imSize(32);
 double *values = new double[imSize*imSize];
 // Define arbitraty values for an image
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   values[j*imSize + i] = double(i-j);
  }
 }
 // Create two images
 Matrix myImage(imSize, imSize, values);
 Matrix myImage2(imSize/2, imSize/2, values);

 // Get the added image
 Matrix myAddedImage = myImage.add(myImage2);
 // Get the substracted image
 Matrix mySubImage = myImage.substract(myImage2);

 BOOST_CHECK_EQUAL(myAddedImage.getXdim(), 0);
 BOOST_CHECK_EQUAL(myAddedImage.getYdim(), 0);
 BOOST_CHECK_EQUAL(mySubImage.getXdim(), 0);
 BOOST_CHECK_EQUAL(mySubImage.getYdim(), 0);

 delete [] values;
 values = nullptr;
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getMax_tests ) {
 logger.info() << "-- Matrix: getMax_test";
 unsigned int imSize(32);
 double *values = new double[imSize*imSize];
 // Define arbitraty values for an image
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   values[j*imSize + i] = double(i+j);
  }
 }
 // Create two images
 Matrix myImage(imSize, imSize, values);

 // Get the max
 double max = myImage.getMax();

 // Check the max is the right value
 BOOST_CHECK_CLOSE(max, 2*(imSize-1), 1);

 delete [] values;
 values = nullptr;
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getMin_tests ) {
 logger.info() << "-- Matrix: getMin_test";
 unsigned int imSize(32);
 double *values = new double[imSize*imSize];
 // Define arbitraty values for an image
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   values[j*imSize + i] = double(i+j);
  }
 }
 // Create two images
 Matrix myImage(imSize, imSize, values);

 // Get the max
 double Min = myImage.getMin();

 // Check the max is the right value
 BOOST_CHECK_CLOSE(Min, 1, 1);

 delete [] values;
 values = nullptr;
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getFlux_tests ) {
 logger.info() << "-- Matrix: getFlux_test";
 unsigned int imSize(32);
 double *values = new double[imSize*imSize];
 // Define arbitraty values for an image
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   values[j*imSize + i] = double(i);
  }
 }
 // Create two images
 Matrix myImage(imSize, imSize, values);

 // Get the max
 double flux = myImage.getFlux();
 double fac = (double (imSize-1)/2)*(double (imSize))*32;
 logger.info() << "-- Matrix: getFlux_test" << fac;
 logger.info() << "-- Matrix: getFlux_test" <<flux;
 // Check the max is the right value
 BOOST_CHECK_CLOSE(flux, fac, 1);

 delete [] values;
 values = nullptr;
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( reset_tests ) {
 logger.info() << "-- Matrix: reset_test";
 unsigned int imSize(32);
 double *values = new double[imSize*imSize];
 // Define arbitraty values for an image
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   values[j*imSize + i] = double(i+j);
  }
 }
 // Create two images
 Matrix myImage(imSize, imSize, values);

 // reset image
 myImage.reset();

 // Check the values
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   BOOST_CHECK_CLOSE(myImage.getValue(i, j), 0., 0.01);
  }
 }

 delete [] values;
 values = nullptr;
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getStandardDeviation_tests ) {
 logger.info() << "-- Matrix: getStandardDeviation_test";
 unsigned int imSize(32);
 double *values = new double[imSize*imSize];
 // Define arbitraty values for an image
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   values[j*imSize + i] = double(i+j);
  }
 }
 // Create two images
 Matrix myImage(imSize, imSize, values);
 // Get the standard deviation
 double stdev = myImage.getStandardDeviation();

 // Check the max is the right value
 BOOST_CHECK_CLOSE(stdev, 13.06394529, 1);

 delete [] values;
 values = nullptr;
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( applyThreshold_tests ) {
 logger.info() << "-- Matrix: applyThreshold_test";
 unsigned int imSize(32);
 double *values = new double[imSize*imSize];
 // Define arbitraty values for an image
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   values[j*imSize + i] = double(i+j);
  }
 }

 // Create an image
 Matrix myImage(imSize, imSize, values);
 // Apply a threshold
 myImage.applyThreshold(10.);

 // Check the values below the threshold are well set to zero
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++){
   if (i+j<10.) {
    BOOST_CHECK_CLOSE(myImage.getValue(i, j), 0, 1);
   } else {
    BOOST_CHECK_CLOSE(myImage.getValue(i, j), i+j, 1);
   }
  }
 }
 delete [] values;
 values = nullptr;
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( localMax_tests ) {
 logger.info() << "-- Matrix: localMax_test";
 unsigned int imSize(32);
 double *values = new double[imSize*imSize];

 // Define arbitraty values for an image
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   values[j*imSize + i] = double(i+j);
  }
 }
 // Create an image
 Matrix myImage(imSize, imSize, values);

 // Check there is no local max in this context
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   BOOST_CHECK(myImage.isLocalMax(i, j)==false);
  }
 }
 // Now add arbitrary local maximums
 myImage.setValue(18, 12, 70);
 myImage.setValue(3, 17, 85);

 // Check those pixels are now local maximum
 BOOST_CHECK(myImage.isLocalMax(18, 12)==true);
 BOOST_CHECK(myImage.isLocalMax(3, 17)==true);
 delete [] values;
 values = nullptr;
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getandsetValues_tests ) {
  logger.info() << "-- Matrix: getandsetValues_test";
 // Create an image initialized with zeros
 unsigned int imSize(32);
 Matrix myImage(imSize, imSize);

 // Set the values
 for (unsigned int i=0; i<imSize; i++){
  for (unsigned int j=0; j<imSize; j++){
   myImage.setValue(i, j, 2*i+j);
  }
 }
 // Check the values are well set
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   BOOST_CHECK_CLOSE(myImage.getValue(i, j), 2*i+j, 0.01);
  }
 }
 // Check the values are well set
 for (unsigned int i=0; i<imSize; i++) {
  BOOST_CHECK_CLOSE(myImage.getValue(i, imSize), 2*i+imSize-1, 0.01);
 }
 for (unsigned int j=0; j<imSize; j++) {
  BOOST_CHECK_CLOSE(myImage.getValue(imSize, j), 2*(imSize-1)+j, 0.01);
 }
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( getArray_tests ) {
 logger.info() << "-- Matrix: getArray_test";
 unsigned int imSize(32);
 double *values = new double[imSize*imSize];
 // Define arbitraty values for an image
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   values[j*imSize + i] = double(i-j);
  }
 }
 // Create two images
 Matrix myImage(imSize, imSize, values);
 // Check the getArray method provided the good values
 double *myArray = myImage.getArray();
 // Define arbitraty values for an image
 for (unsigned int i=0; i<imSize; i++) {
  for (unsigned int j=0; j<imSize; j++) {
   BOOST_CHECK_CLOSE(myArray[j*imSize + i], values[j*imSize + i], 0.01);
  }
 }

 delete [] values;
 values = nullptr;
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
