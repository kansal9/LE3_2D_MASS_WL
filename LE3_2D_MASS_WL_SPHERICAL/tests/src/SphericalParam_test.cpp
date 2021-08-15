/*
 * Copyright (C) 2012-2020 Euclid Science Ground Segment
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
 */

/**
 * @file tests/src/SphericalParam_test.cpp
 * @date 07/14/19
 * @author user
 */

#include <boost/test/unit_test.hpp>
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"
#include "ElementsKernel/Auxiliary.h"
#include <boost/filesystem.hpp>
#include <ios>
#include <sstream>
#include <iostream>

using namespace LE3_2D_MASS_WL_SPHERICAL;
using namespace Euclid::WeakLensing::TwoDMass;
using std::string;
using boost::filesystem::path;

struct SphericalParamTestEnv {
  std::string xmlFileName;
  SphericalParamTestEnv() {
    xmlFileName = Elements::getAuxiliaryPath("LE3_2D_MASS_WL_UTILITIES/Sparam_test.xml").generic_string();
  }
  // delete all test products
  ~SphericalParamTestEnv() {};
};
//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_SUITE (SphericalParam_test)
BOOST_FIXTURE_TEST_SUITE (SphericalParam_test, SphericalParamTestEnv)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( Overall_test ) {
  std::cout<<" --> SphericalParam: Overall_test"<<std::endl;
  int NInpaint(5);
  int NItReducedShear(3);
  int NInpScales(2);
  int nside(256);
  int Nbins(5);
  long BmodesZeros(0);
  long EqualVarPerScale(1);
  float sigmaGauss(4.2);
  float RSsigmaGauss(1.2);
  float threshold(1.3);
  double Zmin(4.32);
  double Zmax(8.54);
  long balancedBins(1);
  int NResamples(30);
  std::string extname = "KAPPA_SPHERE";

  SphericalParam param(nside, NItReducedShear, NInpaint, BmodesZeros, EqualVarPerScale,
                       NInpScales, Nbins, Zmin, Zmax, balancedBins, extname, RSsigmaGauss, sigmaGauss,
                        threshold, NResamples);
  BOOST_CHECK_EQUAL(param.getNside(), 256);
  BOOST_CHECK_EQUAL(param.getNItReducedShear(), 3);
  BOOST_CHECK_EQUAL(param.getNInpaint(), 5);
  if (param.getBmodesZeros() == 0) {
  BOOST_CHECK(true);
  } else {
  BOOST_CHECK(false);
  }
  if (param.getEqualVarPerScale() == 1) {
  BOOST_CHECK(true);
  } else {
  BOOST_CHECK(false);
  }
  if (param.getBalancedBins() == 1) {
  BOOST_CHECK(true);
  } else {
  BOOST_CHECK(false);
  }
  BOOST_CHECK_EQUAL(param.getNInpScales(), 2);
  BOOST_CHECK_CLOSE(param.getSigmaGauss(), sigmaGauss, 1.0e-07);
  BOOST_CHECK_CLOSE(param.getThresholdFDR(), threshold, 0.1);
  BOOST_CHECK_CLOSE(param.getRSSigmaGauss(), RSsigmaGauss, 1.0e-07);
  BOOST_CHECK_EQUAL(param.getNbins(), 5);
  BOOST_CHECK_EQUAL(param.getNResamples(), 30);
  BOOST_CHECK_CLOSE(param.getZMin(), Zmin, 0.1);
  BOOST_CHECK_CLOSE(param.getZMax(), Zmax, 0.1);

}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
 BOOST_FIXTURE_TEST_CASE(xmlFile_test, SphericalParamTestEnv) {
  std::cout << "-- Spherical_PARAM: .xml file";
  SphericalParam sphParam;
   if (true == fileHasField(xmlFileName, "DpdTwoDMassParamsConvergenceSphere")) {
    std::cout<<"Parameter file is for Convergence Sphere.."<<std::endl;
    sphParam.getConvergenceSphereParam(xmlFileName);
     int nside=sphParam.getNside();
     BOOST_CHECK_EQUAL(nside, 16);
   }
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
