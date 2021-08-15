/**
 * @file tests/src/GetSphericalMCMaps_test.cpp
 * @date 10/13/20
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

#include "LE3_2D_MASS_WL_SPHERICAL/GetSphericalMCMaps.h"
#include "ElementsKernel/Logging.h"
#include <boost/filesystem.hpp>
#include "ElementsKernel/Auxiliary.h"
#include "ElementsKernel/Temporary.h"

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_SPHERICAL/Sph_map_maker.h"
#include "LE3_2D_MASS_WL_SPHERICAL/Sph_mass_mapping.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalUtils.h"

#include <ios>
#include <sstream>
#include <iostream>

using std::string;
using namespace Euclid::WeakLensing::TwoDMass;
using namespace LE3_2D_MASS_WL_SPHERICAL;
using namespace Euclid::WeakLensing::TwoDMass::Spherical;
using Euclid::WeakLensing::TwoDMass::Spherical::SphericalIO;

 // handle on created path names
 boost::filesystem::path test_path;
using boost::filesystem::path;
static Elements::Logging logger = Elements::Logging::getLogger("SphericalMCMaps_Test");

struct GetSphericalMCMapsTestEnv {
  std::vector<fs::path> MCShearMaps;
  std::string testParamFile;
  SphericalParam params;
  GetSphericalMCMapsTestEnv () {
    testParamFile = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/Sparam_test.xml").generic_string();
    if (true == fileHasField(testParamFile, "DpdTwoDMassParamsConvergenceSphere")) {
        std::cout<<"Parameter file is for Convergence Sphere.."<<std::endl;
        params.getConvergenceSphereParam (testParamFile);
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
  }
  // delete all test products
  ~GetSphericalMCMapsTestEnv () {};
  std::vector <double> ra;
  std::vector <double> dec;
  std::vector <double> z;
  std::vector <double> kappa;
  std::vector <double> gamma1;
  std::vector <double> gamma2;
  std::vector <double> weight;
  std::vector <std::vector <double> > testCatData;
};
//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_SUITE (GetSphericalMCMaps_test)
BOOST_FIXTURE_TEST_SUITE (GetSphericalMCMaps_test, GetSphericalMCMapsTestEnv)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( getDeNoisedShearMap_test  ) {
  logger.info() << "-- GetSphericalMCMaps: getDeNoisedShearMap_test";
  GetSphericalMCMaps smc(testCatData, params);
  std::pair<Healpix_Map<double>, Healpix_Map<double> > DenoisedShearMapPair = smc.getDeNoisedShearMap();
  BOOST_CHECK(true);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( getNoisedShearMaps_test  ) {
  logger.info() << "-- GetSphericalMCMaps: getNoisedShearMaps_test";
  GetSphericalMCMaps smc(testCatData, params);
  std::vector<std::pair<Healpix_Map<double>, Healpix_Map<double> > > ShearMapList;
  ShearMapList = smc.getNoisedShearMaps();
  BOOST_CHECK(true);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( performAddition_test  ) {
  logger.info() << "-- GetSphericalMCMaps: performAddition_test";
  GetSphericalMCMaps smc(testCatData, params);

  std::pair<Healpix_Map<double>, Healpix_Map<double> > DenoisedShearMapPair = smc.getDeNoisedShearMap();

  std::vector<std::pair<Healpix_Map<double>, Healpix_Map<double> > > NoisedShearMapList;
  NoisedShearMapList = smc.getNoisedShearMaps();

  Elements::TempDir one;
  test_path = one.path();

  for (size_t i = 0; i<NoisedShearMapList.size(); i++) {
   boost::filesystem::path shearMapName =  (fs::path ("EUC_LE3_WL_ShearMap_0" +
                                 std::to_string(i)+ ".fits"));
   MCShearMaps.push_back(test_path/shearMapName);
  }

  for (size_t i = 0; i < NoisedShearMapList.size(); i++) {
    smc.performAddition(DenoisedShearMapPair, NoisedShearMapList[i], MCShearMaps[i].native()); 
  }

   for(size_t i =0; i<MCShearMaps.size(); i++) {
    BOOST_CHECK(boost::filesystem::exists(MCShearMaps[i]));
   }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
