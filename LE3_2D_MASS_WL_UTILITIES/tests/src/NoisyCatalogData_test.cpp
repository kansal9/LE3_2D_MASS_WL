/**
 * @file tests/src/NoisyCatalogData_test.cpp
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

#include "LE3_2D_MASS_WL_UTILITIES/NoisyCatalogData.h"
#include <functional>
#include "ElementsKernel/Logging.h"
#include <fstream>
#include <string>
#include <vector>

using namespace Euclid::WeakLensing::TwoDMass;

static Elements::Logging logger = Elements::Logging::getLogger("Radomization_test");

#define BOOST_CHECK_NOT_EQUAL( a, b ) BOOST_CHECK_PREDICATE( _1 != _2, (a)(b))

//-----------------------------------------------------------------------------

struct NoisyCatalogDataFixture {
  NoisyCatalogDataFixture()
  {
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
   testCatData.push_back(ra);
   testCatData.push_back(dec);
   testCatData.push_back(kappa);
   testCatData.push_back(gamma1);
   testCatData.push_back(gamma2);
   testCatData.push_back(z);
   testCatData.push_back(weight);
  }
  ~NoisyCatalogDataFixture ()
  { }
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
BOOST_FIXTURE_TEST_SUITE (NoisyCatalogData_test, NoisyCatalogDataFixture)
//BOOST_AUTO_TEST_SUITE (NoisyCatalogData_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( createNoisyData_test ) {
  logger.info() << "-- NoisyCatalogData: createNoisyData";
  NoisyCatalogData nData;
  logger.info()<<"size: "<<testCatData.size();
  std::vector<std::vector<double> > testNoisyData;
  testNoisyData = nData.create_noisy_data(testCatData);
  try {
   for (size_t i = 0; i<testNoisyData[0].size(); i++){
      //BOOST_CHECK_CLOSE( testNoisyData[3][i], testCatData[3][i], 0.01 );
     if (testNoisyData[3][i] == 0 && testCatData[3][i] == 0) {
        BOOST_CHECK_EQUAL( testNoisyData[3][i], testCatData[3][i] ); //Check better way
     } else {
     BOOST_CHECK_PREDICATE( std::not_equal_to<double>(), (testNoisyData[3][i])(testCatData[3][i]) );
     //BOOST_CHECK_NOT_EQUAL( testNoisyData[3][i], testCatData[3][i] );
     }
   }
  } catch (...) {
    BOOST_FAIL("Exception in createNoisyData_test");
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
