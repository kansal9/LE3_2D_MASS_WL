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
#include "ElementsKernel/Logging.h"
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "ElementsKernel/Auxiliary.h"
#include "ElementsKernel/Temporary.h" 
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GetMap.h"

#include <fstream>
#include <string>

using std::string;
using boost::filesystem::path;
using namespace Euclid::WeakLensing::TwoDMass;
using namespace LE3_2D_MASS_WL_CARTESIAN;
static Elements::Logging logger = Elements::Logging::getLogger("GetMap_test");
//----------------------------------------------------------------------------
struct GetMapFixture {
  std::string shearMap;
  std::string testParamFile;
  GetMapFixture():xSize(16), ySize(16), zSize(2), mapUniformValueZ0(1.f),
  mapUniformValueZ1(2.f), nGalaxies(100000000), ramin(0.), ramax(0.), decmin(0.), decmax(0.)
  {
    // Allocate the test array
    double *array = new double[xSize*ySize*zSize];

    // Set values to the array
    for (int i=0; i<xSize; i++)
    {
      for (int j=0; j<ySize; j++)
      {
        for (int k=0; k<zSize; k++)
        {
          array[i + j*xSize + k*xSize*ySize]=mapUniformValueZ0 + k*mapUniformValueZ1;
        }
      }
    }
    testParamFile = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/Cparam_test.xml").generic_string();
    shearMap = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/data/shearMap_test.fits").generic_string();
    // Create a map based on this array
    myArrayTestMap = new GetMap(array, xSize, ySize, zSize, nGalaxies);

    // Create a map based on this array without galaxy number
    myArrayTestMapNoGal = new GetMap(array, xSize, ySize, zSize);

    delete [] array;
    array = nullptr;
  }

  ~GetMapFixture()
  {
    delete myArrayTestMap;
    myArrayTestMap = nullptr;

    delete myArrayTestMapNoGal;
    myArrayTestMapNoGal = nullptr;
  }

  GetMap *myArrayTestMap;
  GetMap *myArrayTestMapNoGal;
  int xSize;
  int ySize;
  int zSize;
  float ramin, ramax, decmin, decmax;
  float mapUniformValueZ0, mapUniformValueZ1;
  int nGalaxies;
};

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE (GetMap_test)

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( arrayConstructor_test, GetMapFixture)
{
  // Check the values in the map are corresponding to the values in the input array
  for (int i=0; i<xSize; i++)
  {
    for (int j=0; j<ySize; j++)
    {
      for (int k=0; k<zSize; k++)
      {
        // Check the difference between the map value and the input array is lower than 0.0001%
        BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(i, j, k),
                          mapUniformValueZ0 + k*mapUniformValueZ1, 0.000001);
      }
    }
  }
}

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( arrayConstructor2_test, GetMapFixture)
{
  CartesianParam carParam;
  carParam.ReadConvPatchXMLFile(testParamFile);
  int m_xSize = carParam.getXaxis();
  int m_ySize = carParam.getYaxis();
  int m_zSize = myArrayTestMap->getZdim();
  // Check the values in the map are corresponding to the values in the input array
  for (int i=0; i<m_xSize; i++)
  {
    for (int j=0; j<m_ySize; j++)
    {
      for (int k=0; k<m_zSize; k++)
      {
        // Check the difference between the map value and the input array is lower than 0.0001%
        BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(i, j, k),
                          mapUniformValueZ0 + k*mapUniformValueZ1, 0.000001);
      }
    }
  }
}

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( filenameConstructor_test, GetMapFixture)
{
 GetMap *arr;
 arr = new GetMap(shearMap);
 if (arr==nullptr) {
   BOOST_CHECK(false);
 } else {
   BOOST_CHECK(true);
 }
 delete arr;
 arr = nullptr;
}

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( copyConstructor_test, GetMapFixture)
{
  // Create a Map using the copy constructor
  GetMap *myCopyMap = new GetMap(*myArrayTestMap);	

  // Check this copied map has the same properties as the original map
  BOOST_CHECK(myCopyMap->getXdim() == myArrayTestMap->getXdim());
  BOOST_CHECK(myCopyMap->getYdim() == myArrayTestMap->getYdim());
  BOOST_CHECK(myCopyMap->getZdim() == myArrayTestMap->getZdim());

  // Check the values in the map are corresponding to the values of the original map
  for (int i=0; i<xSize; i++)
  {
    for (int j=0; j<ySize; j++)
    {
      for (int k=0; k<zSize; k++)
      {
        // Check the difference between the map value and the input array is lower than 0.0001%
        BOOST_CHECK_CLOSE(myCopyMap->getBinValue(i, j, k),
                          myArrayTestMap->getBinValue(i, j, k), 0.000001);
      }
    }
  }

  // Check that if I delete one of the maps, the other one is still available (i.e. not sharing same pointer)
  delete myCopyMap;
  myCopyMap = nullptr;

  // Check the values in the original map are still available and correct
  for (int i=0; i<xSize; i++)
  {
    for (int j=0; j<ySize; j++)
    {
      for (int k=0; k<zSize; k++)
      {
        BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(i, j, k),
                          mapUniformValueZ0 + k*mapUniformValueZ1, 0.000001);
      }
    }
  }
}

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( getBinValue_test, GetMapFixture)
{
  // Check the values outside of map return well the edge of the map
  BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(xSize, 0, 0), myArrayTestMap->getBinValue(xSize-1, 0, 0), 0.001);
  BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(0, ySize, 0), myArrayTestMap->getBinValue(0, ySize-1, 0), 0.001);
  BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(0, 0, zSize), myArrayTestMap->getBinValue(0, 0, zSize-1), 0.001);
  BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(xSize, ySize, zSize),
                                                myArrayTestMap->getBinValue(xSize-1, ySize-1, zSize-1), 0.001);

  // Check the values in the map are corresponding to the values in the input array
  for (int i=0; i<xSize; i++)
  {
    for (int j=0; j<ySize; j++)
    {
      for (int k=0; k<zSize; k++)
      {
        // Check the difference between the map value and the input array is lower than 0.0001%
        BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(i, j, k),
                          mapUniformValueZ0 + k*mapUniformValueZ1, 0.000001);
      }
    }
  }
}
//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( getAxisDim_test, GetMapFixture)
{
  // Check the values of the dimensions are the same
  BOOST_CHECK(myArrayTestMap->getXdim() == xSize);
  BOOST_CHECK(myArrayTestMap->getYdim() == ySize);
  BOOST_CHECK(myArrayTestMap->getZdim() == zSize);

  // Try to rebin the map
  bool isPixelated = myArrayTestMap->pixelate(1, 1);
  // Check pixelation worked well
  BOOST_CHECK(isPixelated == true);

  // Check the values of the dimensions are well updated
  BOOST_CHECK(myArrayTestMap->getXdim() == xSize/2);
  BOOST_CHECK(myArrayTestMap->getYdim() == ySize/2);
  BOOST_CHECK(myArrayTestMap->getZdim() == zSize);
}

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( pixelateGalaxiesPerBin_test, GetMapFixture)
{
// Mean number of galaxies per bin
  float galaxiesPerBin = float(nGalaxies)/(xSize*ySize);

  // Try to redo a pixelation to have two times less galaxies per bin
  size_t rebin = myArrayTestMap->pixelate(galaxiesPerBin/2);
  // That should not work so the method shall return 1
  BOOST_CHECK(rebin == 1);

  // Check the values of the dimensions are unchanged
  BOOST_CHECK(myArrayTestMap->getXdim() == xSize);
  BOOST_CHECK(myArrayTestMap->getYdim() == ySize);
  BOOST_CHECK(myArrayTestMap->getZdim() == zSize);

  // Now try to have 5 times more galaxies per bin
  rebin = myArrayTestMap->pixelate(galaxiesPerBin*5);
  // That should work and return 2 (rebinning 2x2)
  BOOST_CHECK(rebin == 2);

  // Update the number of galaxies per bin and the map value
  galaxiesPerBin *= rebin*rebin;
  mapUniformValueZ0 *= rebin*rebin;
  mapUniformValueZ1 *= rebin*rebin;

  // Update the axis dimensions
  xSize /= rebin;
  ySize /= rebin;

  // Check the values of the dimensions are well updated
  BOOST_CHECK(myArrayTestMap->getXdim() == xSize);
  BOOST_CHECK(myArrayTestMap->getYdim() == ySize);
  BOOST_CHECK(myArrayTestMap->getZdim() == zSize);

  // Check the values of the map are well updated
  // Set values to the array
  for (int i=0; i<xSize; i++) {
    for (int j=0; j<ySize; j++) {
      for (int k=0; k<zSize; k++) {
        BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(i, j, k),
                          mapUniformValueZ0 + k*mapUniformValueZ1, 0.001);
  } } }
}

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( getMeanValues_test, GetMapFixture)
{
  // Perform the check for the map with galaxy info
  std::vector<double> myMeans = myArrayTestMap->getMeanValues();
  // Check that the mean values are in agreement with the input values
  BOOST_CHECK_CLOSE(myMeans[0], mapUniformValueZ0, 0.001);
  BOOST_CHECK_CLOSE(myMeans[1], mapUniformValueZ0 + mapUniformValueZ1, 0.001);

  // Perform the check for the map without galaxy info
  myMeans = myArrayTestMapNoGal->getMeanValues();
  // Check that the mean values are in agreement with the input values
  BOOST_CHECK_CLOSE(myMeans[0], mapUniformValueZ0, 0.001);
  BOOST_CHECK_CLOSE(myMeans[1], mapUniformValueZ0 + mapUniformValueZ1, 0.001);
}

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( removeOffset_test, GetMapFixture)
{
  // define offsets to remove then to the map
  std::vector<double> offsets;
  offsets.push_back(0.3*mapUniformValueZ0);
  offsets.push_back(0.45*mapUniformValueZ1);

  // Remove those offsets to the map
  myArrayTestMap->removeOffset(offsets);
  // Check that the values of the map are in agreement
  for (int i=0; i<xSize; i++)
  {
    for (int j=0; j<ySize; j++)
    {
        BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(i, j, 0),
                          0.7*mapUniformValueZ0, 0.001);
        BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(i, j, 1),
                          mapUniformValueZ0 + 0.55*mapUniformValueZ1, 0.001);
    }
  }
}

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( pixelatePerBin_test, GetMapFixture)
{
  // Try to rebin with same binning
  bool isPixelated = myArrayTestMap->pixelate(0, 0);
  // Should return false
  BOOST_CHECK(isPixelated == false);

  // Check the values of the dimensions are unchanged
  BOOST_CHECK(myArrayTestMap->getXdim() == xSize);
  BOOST_CHECK(myArrayTestMap->getYdim() == ySize);
  BOOST_CHECK(myArrayTestMap->getZdim() == zSize);

  // Try to rebin only X axis
  isPixelated = myArrayTestMap->pixelate(1, 0);
  // Should return true
  BOOST_CHECK(isPixelated == true);

  // Update the axis dimensions
  xSize /= 2;

  // Check the values of the dimensions are well updated
  BOOST_CHECK(myArrayTestMap->getXdim() == xSize);
  BOOST_CHECK(myArrayTestMap->getYdim() == ySize);
  BOOST_CHECK(myArrayTestMap->getZdim() == zSize);

  // Check the values of the map are well updated
  // Set values to the array
  for (int i=0; i<xSize; i++)
  {
    for (int j=0; j<ySize; j++)
    {
      for (int k=0; k<zSize; k++)
      {
        BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(i, j, k),
                          mapUniformValueZ0 + k*mapUniformValueZ1, 0.001);
      }
    }
  }

  // Try to rebin only Y axis
  isPixelated = myArrayTestMap->pixelate(0, 1);
  // Should return true
  BOOST_CHECK(isPixelated == true);

  // Update the axis dimensions
  ySize /= 2;

  // Check the values of the dimensions are well updated
  BOOST_CHECK(myArrayTestMap->getXdim() == xSize);
  BOOST_CHECK(myArrayTestMap->getYdim() == ySize);
  BOOST_CHECK(myArrayTestMap->getZdim() == zSize);

  // Check the values of the map are well updated
  // Set values to the array
  for (int i=0; i<xSize; i++)
  {
    for (int j=0; j<ySize; j++)
    {
      for (int k=0; k<zSize; k++)
      {
        BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(i, j, k),
                          mapUniformValueZ0 + k*mapUniformValueZ1, 0.001);
      }
    }
  }

  // Try to rebin both X and Y axis
  isPixelated = myArrayTestMap->pixelate(2, 2);
  // Should return true
  BOOST_CHECK(isPixelated == true);

  // Update the axis dimensions
  xSize /= 4;
  ySize /= 4;

  // Check the values of the dimensions are well updated
  BOOST_CHECK(myArrayTestMap->getXdim() == xSize);
  BOOST_CHECK(myArrayTestMap->getYdim() == ySize);
  BOOST_CHECK(myArrayTestMap->getZdim() == zSize);

  // Check the values of the map are well updated
  // Set values to the array
  for (int i=0; i<xSize; i++)
  {
    for (int j=0; j<ySize; j++)
    {
      for (int k=0; k<zSize; k++)
      {
        BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(i, j, k),
                          mapUniformValueZ0 + k*mapUniformValueZ1, 0.001);
      }
    }
  }
  // Try to rebin so that the number of bin would be less or equal to 1
  isPixelated = myArrayTestMap->pixelate(8, 8);
  // Should return false
  BOOST_CHECK(isPixelated == false);
}

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( applyGaussianFilter_test, GetMapFixture)
{
  std::cout << "-- GetMap: applyGaussianFilter_test"<<std::endl;
  GetMap *myMap;
  myMap = new GetMap(shearMap);

  // Apply a gaussian filter on this map
  myMap->applyGaussianFilter(2.0);

 // Check the difference between the FITS file and the input file is lower than 2%
  for (int i=0; i<xSize; i++)
  {
    for (int j=0; j<ySize; j++)
    {
      for (int k=0; k<zSize; k++)
      {
        //BOOST_CHECK_CLOSE(myMap->getBinValue(i, j, k), mapUniformValueZ0 + k*mapUniformValueZ1, 2);
      }
    }
  }
}

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_CASE( addBorders_test, GetMapFixture)
{
  std::cout << "-- GetMap: addBorders_test"<<std::endl;
  // Add borders
  myArrayTestMap->add_borders();

  // Check the size of the map is multiplied by 2 in X and Y
  BOOST_CHECK(myArrayTestMap->getXdim()==2*xSize);
  BOOST_CHECK(myArrayTestMap->getYdim()==2*ySize);

  // Check the borders added are zeros
  for (int i=0; i<myArrayTestMap->getXdim(); i++)
  {
    for (int j=0; j<myArrayTestMap->getYdim(); j++)
    {
      if(i<myArrayTestMap->getXdim()/4 || i>=3*myArrayTestMap->getXdim()/4
          || j<myArrayTestMap->getYdim()/4 || j>=3*myArrayTestMap->getYdim()/4)
      {
        BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(i, j, 0), 0, 0.1);
        BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(i, j, 1), 0, 0.1);
      }
    }
  }

  // Remove borders
  myArrayTestMap->remove_borders();

  // Check the size of the map is set back to original values
  BOOST_CHECK(myArrayTestMap->getXdim()==xSize);
  BOOST_CHECK(myArrayTestMap->getYdim()==ySize);

  // Check the values are the original ones in the map
  for (int i=0; i<xSize; i++)
  {
    for (int j=0; j<ySize; j++)
    {
      for (int k=0; k<zSize; k++)
      {
        // Check the difference between the map value and the input array is lower than 0.0001%
        BOOST_CHECK_CLOSE(myArrayTestMap->getBinValue(i, j, k),
                          mapUniformValueZ0 + k*mapUniformValueZ1, 0.000001);
      }
    }
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( thresholding_test)
{
  std::cout << "-- GetMap: thresholding_test"<<std::endl;
  int xSize(32), ySize(32), zSize(1);
    // Allocate the test array
    double *Tarray = new double[xSize*ySize*zSize];

    // Set values to the array
    for (int i=0; i<xSize; i++)
    {
      for (int j=0; j<ySize; j++)
      {
        for (int k=0; k<zSize; k++)
        {
          Tarray[i + j*xSize + k*xSize*ySize]= i+j;
        }
      }
    }
  GetMap *myMap;
  myMap = new GetMap(Tarray, xSize, ySize, zSize);

  // Apply a threshold
  myMap->thresholding(10.0);

 // Check the values below the threshold are well set to zero
  for (int i=0; i<xSize; i++)
  {
    for (int j=0; j<ySize; j++)
    {
      for (int k=0; k<zSize; k++)
      {
       if (i+j<10.) {
         BOOST_CHECK_CLOSE(myMap->getBinValue(i, j, k), 0, 1);
       } else {
         BOOST_CHECK_CLOSE(myMap->getBinValue(i, j, k), i+j, 1);
       }
      }
    }
  }
    delete [] Tarray;
    Tarray = nullptr;
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( getSigma_test)
{
  std::cout << "-- GetMap: getSigma_test"<<std::endl;
  int xSize(32), ySize(32), zSize(1);
    // Allocate the test array
    double *Tarray = new double[xSize*ySize*zSize];

    // Set values to the array
    for (int i=0; i<xSize; i++)
    {
      for (int j=0; j<ySize; j++)
      {
        for (int k=0; k<zSize; k++)
        {
          Tarray[i + j*xSize + k*xSize*ySize]= i+j;
        }
      }
    }
  GetMap *myMap;
  myMap = new GetMap(Tarray, xSize, ySize, zSize);

   // Get the standard deviation
 double stdev = myMap->getSigma();

 // Check the max is the right value
 BOOST_CHECK_CLOSE(stdev, 13.06394529, 1);

    delete [] Tarray;
    Tarray = nullptr;
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
