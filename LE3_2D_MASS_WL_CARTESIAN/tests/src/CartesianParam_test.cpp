/**
 * @file tests/src/CartesianParam_test.cpp
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
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "ElementsKernel/Auxiliary.h"
#include "ElementsServices/DataSync.h"
#include <boost/filesystem.hpp>
#include <ios>
#include <sstream>
#include <iostream>

using namespace LE3_2D_MASS_WL_CARTESIAN;
namespace fs = boost::filesystem;
using namespace Euclid::WeakLensing::TwoDMass;
using namespace ElementsServices::DataSync;
using std::string;
using boost::filesystem::path;

//-----------------------------------------------------------------------------

struct CartesianParamDataSyncFixture {
  DataSync sync; // This is the synchronizer
  // These are just shortcuts
  path xmlFileName, clusterxmlFileName, p2sphxmlFileName;
  CartesianParamDataSyncFixture () :
      sync(
          // Here is the connection configuration file
          "LE3_2D_MASS_WL_CARTESIAN/datasync_webdav.conf",
           // Here is the dependency configuration file
          "LE3_2D_MASS_WL_CARTESIAN/test_file_list.txt"),
         xmlFileName(sync.absolutePath("Cparam_test.xml")),
         clusterxmlFileName(sync.absolutePath("clusterParam.xml")),
         p2sphxmlFileName(sync.absolutePath("P2Sph.xml"))
 {
    sync.download();
  }
};

//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_SUITE (CartesianParam_test)
BOOST_FIXTURE_TEST_SUITE (CartesianParam_test, CartesianParamDataSyncFixture)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( Overall_test ) {
  std::cout<<" --> CartesianParam: Overall_test"<<std::endl;
  std::vector<double> zMin;
  zMin.push_back(0.5);
  double zMax(1.2);
  int nbZBins(2);
  int NPatches(2);
  float Pixelsize(0.586);
  float sigmaGauss(0.000001);
  float RSsigmaGauss(0.000001);
  long removeOffset(0);
  long add_borders(0);
  int NInpaint(2);
  int NItReducedShear(2);
  long ForceBMode(0);
  long EqualVarPerScale(0);
  float thresholdFDR(-12.);
  int nbScales(2);
  bool squareMap(false);
  float PatchWidth(10.);
  std::vector<double> mapCenterX;
  mapCenterX.push_back(20.);
  std::vector<double> mapCenterY;
  mapCenterY.push_back(30.);
  long BalancedBins(1);
  int nbSamples(1);
  double zMargin(0.4);
  double massThreshold(0.2);
  std::string extName = "TEST_PATCH";
  std::string ParaFileType = "Conv_Patch";

  CartesianParam carParam(NItReducedShear, NPatches, Pixelsize, PatchWidth, mapCenterX,
                          mapCenterY, nbZBins, zMin, zMax, zMargin, BalancedBins, NInpaint,
                          EqualVarPerScale, ForceBMode, nbScales, add_borders, RSsigmaGauss, sigmaGauss,
                          extName, ParaFileType, massThreshold, thresholdFDR, nbSamples,
                          removeOffset, squareMap);
  std::vector <double> zmin = carParam.getZMin();
  //BOOST_CHECK_EQUAL(zmin[0], 0.1);
  BOOST_CHECK_CLOSE(zmin[0], zMin[0], 0.1);
  //double zmax = carParam.getZMax();
  //BOOST_CHECK_EQUAL(zmax, 0.1);
  //BOOST_CHECK_CLOSE(carParam.getZMin(), zMin, 0.1);
  BOOST_CHECK_CLOSE(carParam.getZMax(), zMax, 0.1);
  BOOST_CHECK_CLOSE(carParam.getSigmaGauss(), sigmaGauss, 0.0000001);
  BOOST_CHECK_CLOSE(carParam.getRSSigmaGauss(), RSsigmaGauss, 1.0e-07);
  BOOST_CHECK_EQUAL(carParam.getPatchWidth(), 10.);
  BOOST_CHECK_EQUAL(carParam.getnbZBins(), 2);
  BOOST_CHECK_EQUAL(carParam.getnbPatches(), 2);
  BOOST_CHECK_EQUAL(carParam.getNInpaint(), 2);
  BOOST_CHECK_EQUAL(carParam.getNSamples(), 1);
  BOOST_CHECK_EQUAL(carParam.getNItReducedShear(), 2);
  BOOST_CHECK_EQUAL(carParam.getnbScales(), 2);
  std::vector <double> MCX = carParam.getMapCenterX();
  BOOST_CHECK_EQUAL(MCX[0], 20.);
  std::vector <double> MCY = carParam.getMapCenterY();
  BOOST_CHECK_EQUAL(MCY[0], 30.);
  BOOST_CHECK_EQUAL(carParam.getSquareMap(), false);
  BOOST_CHECK_EQUAL(carParam.getThreshold(), -12.);
  BOOST_CHECK_EQUAL(carParam.getMassThreshold(), 0.2);
  BOOST_CHECK_EQUAL(carParam.getZMargin(), 0.4);
  double raminn = carParam.getRaMin(MCX[0]);
  BOOST_CHECK_EQUAL(raminn, 15.);
  double decminn = carParam.getDecMin(MCY[0]);
  BOOST_CHECK_EQUAL(decminn, 25.);
  BOOST_CHECK_EQUAL(carParam.getRaMax(raminn), 25.);
  BOOST_CHECK_EQUAL(carParam.getDecMax(decminn), 35.);
  BOOST_CHECK_CLOSE(carParam.getXaxis(), ceil(carParam.getPatchWidth()/carParam.getPixelsize()), 0.000001);
  BOOST_CHECK_CLOSE(carParam.getYaxis(), ceil(carParam.getPatchWidth()/carParam.getPixelsize()), 0.000001);
  if (carParam.getEqualVarPerScale() == 0) {
  BOOST_CHECK(true);
  } else {
  BOOST_CHECK(false);
  }
  if (carParam.get_BalancedBins() == 1) {
  BOOST_CHECK(true);
  } else {
  BOOST_CHECK(false);
  }
  if (carParam.getForceBMode() == 0) {
  BOOST_CHECK(true);
  } else {
  BOOST_CHECK(false);
  }
  if (carParam.get_removeOffset() == 0) {
  BOOST_CHECK(true);
  } else {
  BOOST_CHECK(false);
  }
  if (carParam.get_addBorders() == 0) {
  BOOST_CHECK(true);
  } else {
  BOOST_CHECK(false);
  }

  BOOST_CHECK_CLOSE(carParam.getPixelsize(), (Pixelsize/60.), 0.0001);

  if (carParam.getExtName() == "TEST_PATCH") {
   BOOST_CHECK(true);
  } else {
   BOOST_CHECK(false);
  }
  std::string str = "NEW_TEST_PATCH";
  carParam.setExtName(str);
  if (carParam.getExtName() == "NEW_TEST_PATCH") {
   BOOST_CHECK(true);
  } else {
   BOOST_CHECK(false);
  }

  if (carParam.getParaFileType() == "Conv_Patch") {
   BOOST_CHECK(true);
  } else {
   BOOST_CHECK(false);
  }
}

//-----------------------------------------------------------------------------
 BOOST_FIXTURE_TEST_CASE (ConvergencePatchxmlFile_test, CartesianParamDataSyncFixture) {
  std::cout << "-- CARTESIANPARAM: ConvergencePatch.xml file";
  CartesianParam carParam;
   if (true == fileHasField(xmlFileName, "DpdTwoDMassParamsConvergencePatch")) {
    std::cout<<"Parameter file is for Convergence Patch.."<<std::endl;
    carParam.ReadConvPatchXMLFile(xmlFileName.native());
    BOOST_CHECK_EQUAL(carParam.getPatchWidth(), 10.);
   }
}
//-----------------------------------------------------------------------------
 BOOST_FIXTURE_TEST_CASE (ClusterxmlFile_test, CartesianParamDataSyncFixture) {
  std::cout << "-- CARTESIANPARAM: Cluster.xml file";
  CartesianParam carParam;
   if (true == fileHasField(clusterxmlFileName, "DpdTwoDMassParamsConvergenceClusters")) {
    std::cout<<"Parameter file is for cluster Convergence Patch.."<<std::endl;
    carParam.readConvClustersXMLFile(clusterxmlFileName.native());
    BOOST_CHECK_EQUAL(carParam.getPatchWidth(), 10.);
   }
}
//-----------------------------------------------------------------------------
 BOOST_FIXTURE_TEST_CASE (Patches2SpherexmlFile_test, CartesianParamDataSyncFixture) {
  std::cout << "-- CARTESIANPARAM: Patches2Sphere.xml file";
  CartesianParam carParam;
   if (true == fileHasField(p2sphxmlFileName, "DpdTwoDMassParamsConvergencePatchesToSphere")) {
    std::cout<<"Parameter file is for Patches to sphere.."<<std::endl;
    double placeHolder = 0.;
    carParam.readConvPatchesToSphereXMLFile(p2sphxmlFileName.native(), placeHolder, placeHolder, placeHolder,
                                                                                                 placeHolder);
    BOOST_CHECK_EQUAL(carParam.getNside(), 512);
   }
}

//-----------------------------------------------------------------------------
 BOOST_FIXTURE_TEST_CASE(readParameterFile_test, CartesianParamDataSyncFixture) {
  std::cout << "-- CARTESIANPARAM: readParameterFile_test";
  CartesianParam carParam;
  readParameterFile (fs::path{xmlFileName}, carParam);
  BOOST_CHECK_EQUAL(carParam.getPatchWidth(), 10.);
}

//-----------------------------------------------------------------------------
 BOOST_FIXTURE_TEST_CASE(readParameterFileCluster_test, CartesianParamDataSyncFixture) {
  std::cout << "-- CARTESIANPARAM: readParameterFileCluster_test";
  CartesianParam carParam;
  readParameterFile (fs::path{clusterxmlFileName}, carParam);
  BOOST_CHECK_EQUAL(carParam.getPatchWidth(), 10.);
}

//-----------------------------------------------------------------------------

 BOOST_FIXTURE_TEST_CASE(readParameterFilePatches2Sphere_test, CartesianParamDataSyncFixture) {
  std::cout << "-- CARTESIANPARAM: readParameterFilePatches2Sphere_test";
  CartesianParam carParam;
  readParameterFile (fs::path{p2sphxmlFileName}, carParam, 70., 90., 350., 360.);
  BOOST_CHECK_EQUAL(carParam.getNside(), 512);
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
