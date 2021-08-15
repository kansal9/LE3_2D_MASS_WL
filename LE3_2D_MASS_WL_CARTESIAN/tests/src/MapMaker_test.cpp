/**
 * @file tests/src/MapMaker_test.cpp
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
//#define BOOST_TEST_MODULE MapMaker_test
//#define private public
//#define protected public

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "ElementsKernel/Logging.h"
#include "ElementsKernel/Auxiliary.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MapMaker.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GetMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CoordinateBound.h"
#include <ios>
#include <sstream>
#include <iostream>

using namespace Euclid::WeakLensing::TwoDMass;
using namespace LE3_2D_MASS_WL_CARTESIAN;
using LE3_2D_MASS_WL_CARTESIAN::CartesianParam;
using LE3_2D_MASS_WL_CARTESIAN::CoordinateBound;
using std::string;
using boost::filesystem::path;
static Elements::Logging logger = Elements::Logging::getLogger("MapMaker_test");

//-----------------------------------------------------------------------------


struct MapMakerTestEnv {
  std::string testParamFile;
  CartesianParam params;
  LE3_2D_MASS_WL_CARTESIAN::ShearMap *m_ShearMap = nullptr;
  LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap *m_ConvergenceMap = nullptr;
  MapMakerTestEnv () {
    testParamFile = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/Cparam_test.xml").generic_string();
    if (true == fileHasField(testParamFile, "DpdTwoDMassParamsConvergencePatch")) {
        logger.info() <<"Parameter file is for Convergence Patch..";
        params.ReadConvPatchXMLFile (testParamFile);
    }
  // Assign arbiratary values to required fields
   for (int i = 0; i<= 36; i++) {
    for (int j = 0; j<= 18; j++) {
     ra.push_back(10. * i);
     dec.push_back(10. * j);
     z.push_back(0.5 * j);
     weight.push_back(1.);
     kappa.push_back(0.);
     gamma1.push_back(i/90.);
     gamma2.push_back(j/90.);
    }
   }
   testCatData.push_back(ra);
   testCatData.push_back(dec);
   testCatData.push_back(kappa);
   testCatData.push_back(gamma1);
   testCatData.push_back(gamma2);
   testCatData.push_back(z);
   testCatData.push_back(weight);
   centerX = params.getMapCenterX();
   centerY = params.getMapCenterY();
   centerX.shrink_to_fit();
   centerY.shrink_to_fit();
   rmin = params.getRaMin(centerX[0]);
   dmin = params.getDecMin(centerY[0]);
   zmax = *max_element(testCatData[5].begin(), testCatData[5].end());
   zmin = *min_element(testCatData[5].begin(), testCatData[5].end());

  }
  // delete all test products
  ~MapMakerTestEnv () {};
  std::vector <double> ra;
  std::vector <double> dec;
  std::vector <double> z;
  std::vector <double> kappa;
  std::vector <double> gamma1;
  std::vector <double> gamma2;
  std::vector <double> weight;
  std::vector<double> centerX;
  std::vector<double> centerY;
  std::vector <std::vector <double> > testCatData;
  double rmin, zmax;
  double dmin, zmin;
  //LE3_2D_MASS_WL_CARTESIAN::CoordinateBound m_CB;
};
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE (MapMaker_test)
//BOOST_FIXTURE_TEST_SUITE (MapMaker_test, MapMakerDataSyncFixture)
//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_CASE( get_shearMap_test) {
BOOST_FIXTURE_TEST_CASE( get_shearMap_test, MapMakerTestEnv) {
  logger.info() <<"get_shearMap_test";
  MapMaker map(testCatData,  params);
  std::cout <<"rmin" << rmin<<std::endl;
  std::cout <<"dmin" << rmin<<std::endl;
  CoordinateBound m_CB(rmin, params.getRaMax(rmin), dmin, params.getDecMax(dmin), zmin, zmax);
  m_ShearMap = map.getShearMap(m_CB);
  if (m_ShearMap==nullptr) {
   BOOST_CHECK(false);
  } else {
   BOOST_CHECK(true);
  }
  delete m_ShearMap;
  m_ShearMap = nullptr;
}
//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_CASE( get_convMap_test) {
BOOST_FIXTURE_TEST_CASE( get_convMap_test, MapMakerTestEnv) {
  logger.info() <<"get_convMap_test";
  MapMaker map(testCatData,  params);
  CoordinateBound m_CB(rmin, params.getRaMax(rmin), dmin, params.getDecMax(dmin), zmin, zmax);
  m_ConvergenceMap = map.getConvMap(m_CB);
  if (m_ConvergenceMap==nullptr) {
   BOOST_CHECK(false);
  } else {
   BOOST_CHECK(true);
  }
  delete m_ConvergenceMap;
  m_ConvergenceMap = nullptr;
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
