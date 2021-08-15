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
#include "LE3_2D_MASS_WL_CARTESIAN/GetMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "ElementsKernel/Auxiliary.h"
#include "ElementsKernel/Temporary.h"
#include <fstream>
#include <string>

using std::string;
using boost::filesystem::path;
using namespace Euclid::WeakLensing::TwoDMass;
using namespace LE3_2D_MASS_WL_CARTESIAN;
//----------------------------------------------------------------------------
struct ShearMapFixture {
  std::string shearMap;
  std::string testParamFile;
  CartesianParam params;
  ShearMapFixture():xSize(16), ySize(16), zSize(2), mapUniformValueZ0(1.f),
  mapUniformValueZ1(2.f), nGalaxies(0), ramin(0.), ramax(0.), decmin(0.), decmax(0.)
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
    if (true == fileHasField(testParamFile, "DpdTwoDMassParamsConvergencePatch")) {
        std::cout<<"Parameter file is for Convergence Patch.."<<std::endl;
        params.ReadConvPatchXMLFile (testParamFile);
    }
    shearMap = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/data/shearMap_test.fits").generic_string();
    // Create a map based on this array
    myArrayTestMap = new ShearMap(array, xSize, ySize, zSize, nGalaxies);

    // Create a map based on this array without galaxy number
    myArrayTestMapNoGal = new ShearMap(array, xSize, ySize, zSize);

    delete [] array;
    array = nullptr;
  }

  ~ShearMapFixture()
  {
    delete myArrayTestMap;
    myArrayTestMap = nullptr;

    delete myArrayTestMapNoGal;
    myArrayTestMapNoGal = nullptr;
  }

  ShearMap *myArrayTestMap;
  ShearMap *myArrayTestMapNoGal;
  int xSize;
  int ySize;
  int zSize;
  float ramin, ramax, decmin, decmax;
  float mapUniformValueZ0, mapUniformValueZ1;
  int nGalaxies;
};
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE (ShearMap_test)

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( arrayConstructor_test, ShearMapFixture)
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

BOOST_FIXTURE_TEST_CASE( arrayConstructor2_test, ShearMapFixture)
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

BOOST_FIXTURE_TEST_CASE( filenameConstructor_test, ShearMapFixture)
{
 ShearMap *arr;
 arr = new ShearMap(shearMap);
 if (arr==nullptr) {
   BOOST_CHECK(false);
 } else {
   BOOST_CHECK(true);
 }
 delete arr;
 arr = nullptr;
}
//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( copyConstructor_test, ShearMapFixture)
{
  // Create a Map using the copy constructor
  ShearMap *myCopyMap = new ShearMap(*myArrayTestMap);	

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

BOOST_FIXTURE_TEST_CASE( getConvergence_test, ShearMapFixture)
{
 ShearMap *m_shearMap(nullptr);
 m_shearMap = new ShearMap(shearMap);
 ConvergenceMap *m_convMap(nullptr);
 m_convMap = new ConvergenceMap(m_shearMap->getConvMap());
 if (m_convMap==nullptr) {
   BOOST_CHECK(false);
 } else {
   BOOST_CHECK(true);
 }
 delete m_shearMap;
 delete m_convMap;
 m_shearMap = nullptr;
 m_convMap = nullptr;
}
//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( computeReducedShear_test, ShearMapFixture)
{
 ShearMap *m_reducedShear(nullptr);
 m_reducedShear = new ShearMap(shearMap);
 ConvergenceMap *m_convMap(nullptr);
 m_convMap = new ConvergenceMap(m_reducedShear->getConvMap());
 int reducedIter = params.getNItReducedShear();
 for (int it = 0; it < reducedIter; it++) {
  m_reducedShear->computeReducedShear(*m_convMap);
 }

 if (m_reducedShear==nullptr) {
   BOOST_CHECK(false);
 } else {
   BOOST_CHECK(true);
 }

 delete m_reducedShear;
 delete m_convMap;
 m_reducedShear = nullptr;
 m_convMap = nullptr;
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
