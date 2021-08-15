/**
 * @file tests/src/MassAperturePeakCount_test.cpp
 * @date 10/07/20
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
#include <boost/filesystem.hpp>
#include "ElementsKernel/Auxiliary.h"
#include "ElementsKernel/Temporary.h"
#include "ElementsServices/DataSync.h"

#include "LE3_2D_MASS_WL_PEAK_COUNT/MassAperturePeakCount.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GetMap.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_PEAK_COUNT/PeakParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "LE3_2D_MASS_WL_UTILITIES/NoisyCatalogData.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CoordinateBound.h"

#include <ios>
#include <sstream>
#include <iostream>

using namespace Euclid::WeakLensing::TwoDMass;
using namespace ElementsServices::DataSync;
using std::string;
using LE3_2D_MASS_WL_CARTESIAN::CoordinateBound;

 // handle on created path names
 boost::filesystem::path test_path;

//-----------------------------------------------------------------------------

struct MassApertureDataSyncFixture {
  DataSync sync; // This is the synchronizer
  // These are just shortcuts
  path testPatchParamFile, testPeakParamFile, testCatFile;
  MassApertureDataSyncFixture () :
      sync(
          // Here is the connection configuration file
          "LE3_2D_MASS_WL_PEAK_COUNT/datasync_webdav.conf",
           // Here is the dependency configuration file
          "LE3_2D_MASS_WL_PEAK_COUNT/test_file_list.txt"),
         testPatchParamFile(sync.absolutePath("Cparam_test.xml")),
         testPeakParamFile(sync.absolutePath("PMassAp_test.xml")),
         testCatFile(sync.absolutePath("data/lensmcCat.fits"))
 {
    sync.download();
  }
};

//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_SUITE (MassAperturePeakCount_test)
BOOST_FIXTURE_TEST_SUITE (MassAperturePeakCount_test, MassApertureDataSyncFixture)
//-----------------------------------------------------------------------------


BOOST_AUTO_TEST_CASE( OverallMassAperturePeakCount_test) {

 std::cout <<"OverallMassAperturePeakCount_test" << std::endl;

  LE3_2D_MASS_WL_PEAK_COUNT::PeakParam params;
   // Case 1: check file is of XML format
   if (true == checkFileType(testPeakParamFile, Euclid::WeakLensing::TwoDMass::signXML)) {
    // read XML, get parameter file type
    if (true == fileHasField(testPeakParamFile, "DpdTwoDMassParamsPeakCatalogMassAperture")) {
      std::cout << "Parameter file found for Aperture Mass Peak count" << std::endl;
      params.readPeakMassApertXML(testPeakParamFile.native());
    } else {
         std::cout << "Parameter file not found for Aperture Mass Peak count ...." << std::endl;
     }
    }

  LE3_2D_MASS_WL_CARTESIAN::CartesianParam Patchparams;
   // Case 1: check file is of XML format
   if (true == checkFileType(testPatchParamFile, Euclid::WeakLensing::TwoDMass::signXML)) {
    // read XML, get parameter file type
    if (true == fileHasField(testPatchParamFile, "DpdTwoDMassParamsConvergencePatch")) {
      Patchparams.ReadConvPatchXMLFile (testPatchParamFile.native());
    } else {
      std::cout<<"Parameter file for Patch is not in XML format.."<< std::endl;
    }
   }

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Object to fetch Shear Map
////////////////////////////////////////////////////////////////////////////////////////////////////////
Euclid::WeakLensing::TwoDMass::CartesianKS::CartesianAlgoKS CartesainAlgo(Patchparams);

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Reading input shear catalog
////////////////////////////////////////////////////////////////////////////////////////////////////////
  CatalogData readCat;
  std::cout<< "Reading Input Shear Catalog" << std::endl;
  std::vector<std::vector<double> > catalogData;
  //catalogData = readCat.readCatalog(testCatFile);
  readCat.readCatalog(testCatFile.native(), catalogData);
  double zMin, zMax;
  vecMinMax(catalogData[5], &zMin, &zMax);
  std::cout<< "Reading is completed.."<< std::endl;
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Creating Patch of shear Map from Input shear catalog
////////////////////////////////////////////////////////////////////////////////////////////////////////
 std::vector<double> centerX = Patchparams.getMapCenterX();
 std::vector<double> centerY = Patchparams.getMapCenterY();
  Elements::TempDir one;
  test_path = one.path();
 for (size_t it=0; it<centerX.size(); it++){
  double ramin, decmin;
  ramin = Patchparams.getRaMin(centerX[it]);
  decmin = Patchparams.getDecMin(centerY[it]);
  CoordinateBound m_CB(ramin, Patchparams.getRaMax(ramin), decmin,
                     Patchparams.getDecMax(decmin), zMin, zMax);
   std::cout << "creating shear Map"<< std::endl;
   boost::filesystem::path shearMap= "Test_shearMapFile.fits";
   CartesainAlgo.extractShearMap ((test_path/shearMap).native(), catalogData, m_CB);

////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Randomising shear catalog
////////////////////////////////////////////////////////////////////////////////////////////////////////
  NoisyCatalogData random;
  std::vector<std::vector<double> > NoisData;
  NoisData = random.create_noisy_data(catalogData);
  catalogData.clear();
////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Creating Patch of shear Map from randomised shear
////////////////////////////////////////////////////////////////////////////////////////////////////////
   std::cout << "creating Noisy shear Map"<< std::endl;
   boost::filesystem::path NoisyshearMap= "Test_NoisyshearMap.fits";
   CartesainAlgo.extractShearMap ((test_path/NoisyshearMap).native(), NoisData, m_CB);
   NoisData.clear();

////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Output filename
////////////////////////////////////////////////////////////////////////////////////////////////////////
   boost::filesystem::path outputPeakMACatalog= "Test_PeakCatalog.fits";

    LE3_2D_MASS_WL_CARTESIAN::ShearMap *m_shear;
    m_shear = new LE3_2D_MASS_WL_CARTESIAN::ShearMap((test_path/shearMap).native());

    LE3_2D_MASS_WL_CARTESIAN::ShearMap *m_noisedShear;
    m_noisedShear = new LE3_2D_MASS_WL_CARTESIAN::ShearMap((test_path/NoisyshearMap).native());

    // Then create a PeakCountAlgo object
    LE3_2D_MASS_WL_PEAK_COUNT::MassAperturePeakCount myPeakCounter(*m_shear, *m_noisedShear, params);
    myPeakCounter.saveMAPeakCatalog((test_path/outputPeakMACatalog).string());

    delete m_shear;
    m_shear = nullptr;

    delete m_noisedShear;
    m_noisedShear = nullptr;

    BOOST_CHECK(true);
 }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
