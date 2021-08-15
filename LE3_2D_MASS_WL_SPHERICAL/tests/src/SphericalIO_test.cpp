/**
 * @file tests/src/SphericalIO_test.cpp
 * @date 12/24/20
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

#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"
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
using Euclid::WeakLensing::TwoDMass::Spherical::SphericalIO;
using namespace Euclid::WeakLensing::TwoDMass;
using namespace Euclid::WeakLensing::TwoDMass::Spherical;
 // handle on created path names
 boost::filesystem::path test_path;
//-----------------------------------------------------------------------------

struct SphericalIOFixture {
  std::string testParamFile;
  SphericalParam params;
  SphericalIOFixture()
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
  ~SphericalIOFixture ()
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

//BOOST_AUTO_TEST_SUITE (SphericalIO_test)
BOOST_FIXTURE_TEST_SUITE (SphericalIO_test, SphericalIOFixture)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( wrtieMap_test ) {
 std::cout << "-- SphericalIO: writeMap_Test"<<std::endl;

 using Elements::TempDir;
 TempDir one;
 test_path = one.path();

 Sph_map_maker mapMaker(params);
 //std::pair<Healpix_Map<double>, Healpix_Map<double> > shearPair; 
 auto [ Shear1, Shear2, GalCount ] = mapMaker.create_ShearMap(catData);
 SphericalIO SphericalIO(params);
// Save this map back in a new file$
 boost::filesystem::path filename = test_path / "GMap.fits";
 std::cout << "-- SphericalIO: writing Map1"<<std::endl;
 //bool status = SphericalIO.write_Map(filename.native(), shearPair.first, "GAMMA1");
 bool status = SphericalIO.write_Map(filename.native(), Shear1, "GAMMA1");
 BOOST_REQUIRE(status = true);
 std::cout << "-- SphericalIO: writing Map2"<<std::endl;
 //status = SphericalIO.write_Map(filename.native(), shearPair.second, "GAMMA2");
 status = SphericalIO.write_Map(filename.native(), Shear2, "GAMMA2");
 BOOST_REQUIRE(status = true);
 std::cout << "-- SphericalIO: Complete test"<<std::endl;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
