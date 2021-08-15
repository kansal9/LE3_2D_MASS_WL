/**
 * @file tests/src/DmOutput_test.cpp
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

#include "LE3_2D_MASS_WL_UTILITIES/DmOutput.h"

#include <boost/test/unit_test.hpp>

#include <boost/filesystem.hpp>
#include "ElementsKernel/Auxiliary.h"
#include "ElementsKernel/Temporary.h"
#include "ElementsServices/DataSync.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include <fstream>
#include <string>

using namespace Euclid::WeakLensing::TwoDMass;
using std::string;
using namespace ElementsServices::DataSync;
namespace fs = boost::filesystem;
using Elements::TempDir;
 // handle on created path names
 boost::filesystem::path test_path;

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE (DmOutput_test)

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(createNoisedPatchXml_test ) {

std::cout << "-- DmOutput:createNoisedPatchXml_test"<<std::endl;

 try {
Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
  fs::path fitsFile {"testFile.fits"};
  const fs::path file = test_path / "NoisedConvergenceMapPatch.xml";
    std::string product_type = "DpdTwoDMassConvergencePatch";
    int NResamples = 0;

    auto product = initProduct<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch,
    							   pro::le3::wl::twodmass::twoDMassCollectConvergencePatch,
			                       int>(product_type, file, NResamples);

  dm.createNoisedPatchXml(product, fitsFile);
  writeProduct<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch>(product, file);
  BOOST_CHECK(boost::filesystem::is_regular_file(file));
 } catch (...) {
    BOOST_FAIL("Exception in createNoisedPatchXml");
  }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(createDenoisedPatchXml_test ) {

std::cout << "-- DmOutput:createDenoisedPatchXml_test"<<std::endl;
 try {
Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
  fs::path fitsFile {"testFile.fits"};
  const fs::path filename = test_path / "DenoisedConvergenceMapPatch.xml";
    std::string product_type = "DpdTwoDMassConvergencePatch";
    int NResamples = 0;

    auto product = initProduct<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch,
    							   pro::le3::wl::twodmass::twoDMassCollectConvergencePatch,
			                       int>(product_type, filename, NResamples);
  dm.createDenoisedPatchXml(product, fitsFile);
  writeProduct<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch>(product, filename);
  BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createDenoisedPatchXml");
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(createSingleClusterXml_test ) {
std::cout << "-- DmOutput:createSingleClusterXml_test"<<std::endl;
 try {
Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
  fs::path fitsFile {"testFile.fits"};
 const fs::path filename = test_path / "ClusterConvergenceMapPatch.xml";
 std::string product_type = "DpdTwoDMassConvergenceClusters";
    std::string ClusterID = "198";
       auto product = initProduct<dpd::le3::wl::twodmass::out::convergencesinglecluster::dpdTwoDMassConvergenceSingleCluster,
    							    pro::le3::wl::twodmass::twoDMassConvergenceSingleCluster,
			                        std::string>(product_type, filename, ClusterID);
     dm.createSingleClusterXml(product, fitsFile);
  writeProduct<dpd::le3::wl::twodmass::out::convergencesinglecluster::dpdTwoDMassConvergenceSingleCluster>(product, filename);
  BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createSingleClusterXml");
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(createSNRPatchOutputXml_test ) {
std::cout << "-- DmOutput:createSNRPatchOutputXml_test"<<std::endl;
 try {
Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
  fs::path fitsFile {"testFile.fits"};
 const fs::path filename = test_path / "SNRPatch.xml";
    std::string product_type = "DpdTwoDMassConvergencePatch";
    int NResamples = 2;

    auto product = initProduct<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch,
    							   pro::le3::wl::twodmass::twoDMassCollectConvergencePatch,
			                       int>(product_type, filename, NResamples);
  dm.createSNRPatchOutputXml(product, fitsFile);
  writeProduct<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch>(product, filename);
  BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createSNRPatchOutputXml");
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(createPatchOutputXml_test ) {
std::cout << "-- DmOutput:createPatchOutputXml_test"<<std::endl;
 try {
Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
 const fs::path filename = test_path / "Patch.xml";
    std::string product_type = "DpdTwoDMassConvergencePatch";
    int NResamples = 2;

    auto product = initProduct<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch,
    							   pro::le3::wl::twodmass::twoDMassCollectConvergencePatch,
			                       int>(product_type,filename, NResamples);
    fs::path out_fits_file_1 {"some_file_1.fits"};
    dm.createNoisedPatchXml(product, out_fits_file_1);
    // Add something else
    fs::path out_fits_file_2 {"some_file_2.fits"};
    dm.createDenoisedPatchXml(product, out_fits_file_2);
    // Add something else
    fs::path out_fits_file_3 {"some_file_3.fits"};
    dm.createSNRPatchOutputXml(product, out_fits_file_3);
  writeProduct<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch>(product, filename);
   BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createSNRPatchOutputXml");
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(createSphereOutputXml_test ) {
std::cout << "-- DmOutput:createSphereOutputXml_test"<<std::endl;
 try {
Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
 const fs::path filename = test_path / "Sphere.xml";
    std::string product_type = "DpdTwoDMassConvergenceSphere";
    int NResamples = 2;

    auto product = initProduct<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere,
    							    pro::le3::wl::twodmass::twoDMassCollectConvergenceSphere,
			                        int>(product_type, filename, NResamples);
    fs::path out_fits_file_1 {"some_file_1.fits"};
    dm.createNoisedSphereXml(product, out_fits_file_1);

    fs::path out_fits_file_2 {"some_file_2.fits"};
    dm.createDenoisedSphereXml(product, out_fits_file_2);

    fs::path out_fits_file_3 {"some_file_3.fits"};
    dm.createSNRSphereOutputXml(product, out_fits_file_3);

    fs::path out_fits_file_4 {"some_file_4.fits"};
    dm.createSphereGalCountXml(product, out_fits_file_4);
  writeProduct<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere>(product, filename);
   BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createSphereOutputXml");
  }
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(createNoisedSphereXml_test ) {
std::cout << "-- DmOutput:createNoisedSphereXml_test"<<std::endl;
 try {
 Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
  fs::path fitsFileSphere {"testFile.fits"};
 const fs::path filename = test_path / "NoisedConvergenceSphere.xml";
    std::string product_type = "DpdTwoDMassConvergenceSphere";
    int NResamples = 0;
    auto product = initProduct<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere,
    							    pro::le3::wl::twodmass::twoDMassCollectConvergenceSphere,
			                        int>(product_type, filename, NResamples);

  dm.createNoisedSphereXml(product, fitsFileSphere);
  writeProduct<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere>(product, filename);
  BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createNoisedSphereXml");
  }
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(createDenoisedSphereXml_test ) {
std::cout << "-- DmOutput:createDenoisedSphereXml_test"<<std::endl;
 try {
 Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
 fs::path fitsFileSphere {"testFile.fits"};
 const fs::path filename = test_path / "DeNoisedConvergenceSphere.xml";
    std::string product_type = "DpdTwoDMassConvergenceSphere";
    int NResamples = 0;
    auto product = initProduct<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere,
    							    pro::le3::wl::twodmass::twoDMassCollectConvergenceSphere,
			                        int>(product_type, filename, NResamples);

  dm.createDenoisedSphereXml(product, fitsFileSphere);
  writeProduct<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere>(product, filename);
  BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createDenoisedSphereXml");
  }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(createSphereGalCountXml_test ) {
std::cout << "-- DmOutput:createSphereGalCountXml_test"<<std::endl;
 try {
  Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
 fs::path fitsFileSphere {"testFile.fits"};
 const fs::path filename = test_path / "GalCountSphere.xml";
    std::string product_type = "DpdTwoDMassConvergenceSphere";
    int NResamples = 0;
    auto product = initProduct<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere,
    							    pro::le3::wl::twodmass::twoDMassCollectConvergenceSphere,
			                        int>(product_type, filename, NResamples);
  dm.createSphereGalCountXml(product, fitsFileSphere);
  writeProduct<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere>(product, filename);
  BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createSphereGalCountXml");
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(createSNRSphereOutputXml_test ) {
std::cout << "-- DmOutput:createSNRSphereOutputXml_test"<<std::endl;
 try {
   Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
 fs::path fitsFileSphere {"testFile.fits"};
 const fs::path filename = test_path / "SNRSphere.xml";
    std::string product_type = "DpdTwoDMassConvergenceSphere";
    int NResamples = 2;
    auto product = initProduct<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere,
    							    pro::le3::wl::twodmass::twoDMassCollectConvergenceSphere,
			                        int>(product_type, filename, NResamples);
  dm.createSNRSphereOutputXml(product, fitsFileSphere);
  writeProduct<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere>(product, filename);
  BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createSNRSphereOutputXml");
  }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(createPeakOutputXml_test ) {
std::cout << "-- DmOutput: createPeakOutputXml_test"<<std::endl;
 try {
 Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
 fs::path fitsFilePC {"testFile.fits"};
 const fs::path filename = test_path / "PeakCount.xml";
  dm.createPeakOutputXml(filename, fitsFilePC);
  BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createPeakOutputXml");
  }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(createNoisedPatchtoSphereXml_test ) {
std::cout << "-- DmOutput: createNoisedPatchtoSphereXml_test"<<std::endl;
 try {
  Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
 fs::path fitsFile {"testFile.fits"};
 const fs::path filename = test_path / "NoisedPatchtoSphere.xml";
 std::string product_type = "DpdTwoDMassConvergencePatchesToSphere";
 int NResamples = 0;
   auto product = initProduct<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere,
    							    pro::le3::wl::twodmass::twoDMassCollectConvergencePatchesToSphere,
			                        int>(product_type, filename, NResamples);
  dm.createNoisedPatchtoSphereXml(product, fitsFile);
  writeProduct<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>(product, filename);
  BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createNoisedPatchtoSphereXml");
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(createDenoisedPatchtoSphereXml_test ) {
std::cout << "-- DmOutput: createDenoisedPatchtoSphereXml_test"<<std::endl;
 try {
  Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
 fs::path fitsFile {"testFile.fits"};
 const fs::path filename = test_path / "DeNoisedPatchtoSphere.xml";
 std::string product_type = "DpdTwoDMassConvergencePatchesToSphere";
 int NResamples = 0;
   auto product = initProduct<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere,
    							    pro::le3::wl::twodmass::twoDMassCollectConvergencePatchesToSphere,
			                        int>(product_type, filename, NResamples);
  dm.createDenoisedPatchtoSphereXml(product, fitsFile);
  writeProduct<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>(product, filename);
  BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createDenoisedPatchtoSphereXml");
  }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(createPatchtoSphereGalCountXml_test ) {
std::cout << "-- DmOutput: createPatchtoSphereGalCountXml_test"<<std::endl;
 try {
  Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
 const fs::path filename = test_path / "GalCountPatchtoSphere.xml";
 fs::path fitsFile {"testFile.fits"};
 std::string product_type = "DpdTwoDMassConvergencePatchesToSphere";
 int NResamples = 0;
   auto product = initProduct<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere,
    							    pro::le3::wl::twodmass::twoDMassCollectConvergencePatchesToSphere,
			                        int>(product_type, filename, NResamples);
  dm.createPatchtoSphereGalCountXml(product, fitsFile);
  writeProduct<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>(product, filename);
  BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createPatchtoSphereGalCountXml");
  }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(createSNRPatchtoSphereOutputXml_test ) {
std::cout << "-- DmOutput: createSNRPatchtoSphereOutputXml_test"<<std::endl;
 try {
  Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
 fs::path fitsFile {"testFile.fits"};
 const fs::path filename = test_path / "SNRPatchtoSphere.xml";
 std::string product_type = "DpdTwoDMassConvergencePatchesToSphere";
 int NResamples = 2;
   auto product = initProduct<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere,
    							    pro::le3::wl::twodmass::twoDMassCollectConvergencePatchesToSphere,
			                        int>(product_type, filename, NResamples);
  dm.createSNRPatchtoSphereOutputXml(product, fitsFile);

  writeProduct<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>(product, filename);
  BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createSNRPatchtoSphereOutputXml");
  }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(createPatchtoSphereProjCenterPosXml_test ) {
std::cout << "-- DmOutput: createPatchtoSphereProjCenterPosXml_test"<<std::endl;
 try {
  Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
 fs::path fitsFile {"testFile.fits"};
 const fs::path filename = test_path / "ProjCenterPosPatchtoSphere.xml";
 std::string product_type = "DpdTwoDMassConvergencePatchesToSphere";
 int NResamples = 0;
   auto product = initProduct<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere,
    							    pro::le3::wl::twodmass::twoDMassCollectConvergencePatchesToSphere,
			                        int>(product_type, filename, NResamples);
  dm.createPatchtoSphereProjCenterPosXml(product, fitsFile);
  writeProduct<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>(product, filename);
  BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createPatchtoSphereProjCenterPosXml");
  }
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(createPatchtoSphereXml_test ) {
std::cout << "-- DmOutput: createPatchtoSphereXml_test"<<std::endl;
 try {
  Euclid::WeakLensing::TwoDMass::DmOutput dm;
 // block creation for local variables
 TempDir td;
 test_path = td.path();
 fs::path fitsFile1 {"testFile1.fits"};
 fs::path fitsFile2 {"testFile2.fits"};
 fs::path fitsFile3 {"testFile3.fits"};
 fs::path fitsFile4 {"testFile4.fits"};
 fs::path fitsFile5 {"testFile5.fits"};
 const fs::path filename = test_path / "PatchtoSphere.xml";
 std::string product_type = "DpdTwoDMassConvergencePatchesToSphere";
 int NResamples = 2;
   auto product = initProduct<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere,
    							    pro::le3::wl::twodmass::twoDMassCollectConvergencePatchesToSphere,
			                        int>(product_type, filename, NResamples);
  dm.createNoisedPatchtoSphereXml(product, fitsFile1);
  dm.createDenoisedPatchtoSphereXml(product, fitsFile2);
  dm.createPatchtoSphereGalCountXml(product, fitsFile3);
  dm.createSNRPatchtoSphereOutputXml(product, fitsFile4);
  dm.createPatchtoSphereProjCenterPosXml(product, fitsFile5);
  writeProduct<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>(product, filename);
  BOOST_CHECK(boost::filesystem::is_regular_file(filename));
 } catch (...) {
    BOOST_FAIL("Exception in createNoisedPatchtoSphereXml");
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
