/**
 * @file tests/src/GetCartesianMCMaps_test.cpp
 * @date 10/13/20
 * @author vkansal
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

#include "LE3_2D_MASS_WL_CARTESIAN/GetCartesianMCMaps.h"
#include "ElementsKernel/Logging.h"
#include <boost/filesystem.hpp>
#include "ElementsKernel/Auxiliary.h"
#include "ElementsKernel/Temporary.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MassMapping.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CoordinateBound.h"
#include <ios>
#include <sstream>
#include <iostream>


using namespace Euclid::WeakLensing::TwoDMass;
using namespace LE3_2D_MASS_WL_CARTESIAN;
using LE3_2D_MASS_WL_CARTESIAN::CoordinateBound;
using std::string;

 // handle on created path names
 boost::filesystem::path test_path;
static Elements::Logging logger = Elements::Logging::getLogger("CartesianMCMaps_Test");
struct GetCartesianMCMapsTestEnv {
  std::vector<fs::path> MCShearMaps;
  std::string testParamFile;
  CartesianParam params;
  GetCartesianMCMapsTestEnv () {
    testParamFile = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/Cparam_test.xml").generic_string();
    if (true == fileHasField(testParamFile, "DpdTwoDMassParamsConvergencePatch")) {
        std::cout<<"Parameter file is for Convergence Patch.."<<std::endl;
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
  ~GetCartesianMCMapsTestEnv () {};
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
};
//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_SUITE (GetCartesianMCMaps_test)
 BOOST_FIXTURE_TEST_SUITE (GetCartesianMCMaps_test, GetCartesianMCMapsTestEnv)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( getDeNoisedShearMap_test ) {
  std::cout <<"-- GetCartesianMCMaps: getDeNoisedShearMap_test" << std::endl;
  CoordinateBound m_CB(rmin, params.getRaMax(rmin), dmin, params.getDecMax(dmin), zmin, zmax);
  GetCartesianMCMaps mc(testCatData, params, m_CB);
  ShearMap *DeNoisedShearMap;
  DeNoisedShearMap = mc.getDeNoisedShearMap();
  if (DeNoisedShearMap==nullptr) {
   BOOST_CHECK(false);
  } else {
   BOOST_CHECK(true);
  }
  delete DeNoisedShearMap;
  DeNoisedShearMap = nullptr;
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( getNoisedShearMaps_test ) {
  std::cout <<"-- GetCartesianMCMaps: getNoisedShearMaps_test" << std::endl;
  CoordinateBound m_CB(rmin, params.getRaMax(rmin), dmin, params.getDecMax(dmin), zmin, zmax);
  GetCartesianMCMaps mc(testCatData, params, m_CB);
  std::vector<LE3_2D_MASS_WL_CARTESIAN::ShearMap*> NoisedShearMapList;
  NoisedShearMapList = mc.getNoisedShearMaps();
  for (const auto &iter: NoisedShearMapList){
    if (iter==nullptr) {
     BOOST_CHECK(false);
    } else {
     BOOST_CHECK(true);
    }
  }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( performAddition_test ) {
  std::cout <<"-- GetCartesianMCMaps: performAddition_test" << std::endl;
  CoordinateBound m_CB(rmin, params.getRaMax(rmin), dmin, params.getDecMax(dmin), zmin, zmax);
  GetCartesianMCMaps mc(testCatData, params, m_CB);
  ShearMap *DeNoisedShearMap;
  DeNoisedShearMap = mc.getDeNoisedShearMap();
  std::vector<LE3_2D_MASS_WL_CARTESIAN::ShearMap*> NoisedShearMapList;
  NoisedShearMapList = mc.getNoisedShearMaps();
  Elements::TempDir one;
  test_path = one.path();
  for (size_t i = 0; i<NoisedShearMapList.size(); i++) {
   boost::filesystem::path shearMapName =  (fs::path ("EUC_LE3_WL_ShearMap_0" +
                                 std::to_string(i)+ ".fits"));
   MCShearMaps.push_back(test_path/shearMapName);
  }

  for (const auto &iter: NoisedShearMapList){
    auto index = (&iter - &NoisedShearMapList[0]);
    mc.performAddition(*DeNoisedShearMap, *iter, MCShearMaps[index].native() );
  }
   for(size_t i =0; i<MCShearMaps.size(); i++) {
    BOOST_CHECK(boost::filesystem::exists(MCShearMaps[i]));
   }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
