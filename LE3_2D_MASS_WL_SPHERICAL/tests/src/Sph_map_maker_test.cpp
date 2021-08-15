/**
 * @file tests/src/Sph_map_maker_test.cpp
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
#define BOOST_TEST_MODULE Sph_map_maker_test
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "ElementsKernel/Auxiliary.h"
#include "ElementsKernel/Temporary.h"
#include <iostream>
#include "LE3_2D_MASS_WL_SPHERICAL/Sph_map_maker.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"

using LE3_2D_MASS_WL_SPHERICAL::Sph_map_maker;
using LE3_2D_MASS_WL_SPHERICAL::SphericalParam;
using namespace Euclid::WeakLensing::TwoDMass;
 // handle on created path names
 boost::filesystem::path test_path;
//-----------------------------------------------------------------------------

struct SphMapMakerFixture {
  std::string testParamFile;
  SphericalParam params;
  SphMapMakerFixture()
  {
   testParamFile = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/Sparam_test.xml").generic_string();
   if (true == fileHasField(testParamFile, "DpdTwoDMassParamsConvergenceSphere")) {
    std::cout<<"Parameter file is for Convergence Sphere.."<<std::endl;
    params.getConvergenceSphereParam(testParamFile);
   }
   // Assign arbiratary values to required fields
   for (int i = 0; i<= 36; i++) {
    for (int j = 0; j<= 18; j++) {
     ra.push_back(10. * i);
     dec.push_back(10. * j);
     z.push_back(0.5 * j);
     weight.push_back(0.);
     kappa.push_back(0.);
     gamma1.push_back(i/90.);
     gamma2.push_back(j/90.);
    }
   }
   catData.push_back(ra);
   catData.push_back(dec);
   catData.push_back(kappa);
   catData.push_back(gamma1);
   catData.push_back(gamma2);
   catData.push_back(z);
   catData.push_back(weight);
  }
  ~SphMapMakerFixture ()
  { }
  int nside;
  std::vector <double> ra;
  std::vector <double> dec;
  std::vector <double> z;
  std::vector <double> kappa;
  std::vector <double> gamma1;
  std::vector <double> gamma2;
  std::vector <double> weight;
  std::vector <std::vector <double> > catData;
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_SUITE (Sph_map_maker_test)
BOOST_FIXTURE_TEST_SUITE (Sph_map_maker_test, SphMapMakerFixture)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( SphMapMaker_test ) {
  std::cout << "-- SphMapMaker_Test"<<std::endl;
  Sph_map_maker mapMaker(params);
  //std::pair<Healpix_Map<double>, Healpix_Map<double> > shearPair = mapMaker.create_ShearMap(catData);
  auto [ Shear1, Shear2, GalCount ] = mapMaker.create_ShearMap(catData);
  BOOST_CHECK(true);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
