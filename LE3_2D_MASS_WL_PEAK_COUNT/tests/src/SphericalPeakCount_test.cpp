/**
 * @file tests/src/SphericalPeakCount_test.cpp
 * @date 06/18/21
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

#include "LE3_2D_MASS_WL_PEAK_COUNT/SphericalPeakCount.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include <boost/filesystem.hpp>
#include "ElementsKernel/Auxiliary.h"
#include "ElementsKernel/Temporary.h"
#include "ElementsServices/DataSync.h"
#include "LE3_2D_MASS_WL_PEAK_COUNT/PeakParam.h"

#include <ios>
#include <sstream>
#include <iostream>

using namespace Euclid::WeakLensing::TwoDMass;
using namespace ElementsServices::DataSync;
using namespace LE3_2D_MASS_WL_SPHERICAL;
using namespace LE3_2D_MASS_WL_PEAK_COUNT;
using LE3_2D_MASS_WL_PEAK_COUNT::PeakParam;

using std::string;
using boost::filesystem::path;

 // handle on created path names
 boost::filesystem::path test_path;
//-----------------------------------------------------------------------------

struct SphericalPeakDataSyncFixture {
  DataSync sync; // This is the synchronizer
  // These are just shortcuts
  path PtestParamFile, inConvergenceMap;
  SphericalPeakDataSyncFixture () :
      sync(
          // Here is the connection configuration file
          "LE3_2D_MASS_WL_PEAK_COUNT/datasync_webdav.conf",
           // Here is the dependency configuration file
          "LE3_2D_MASS_WL_PEAK_COUNT/test_file_list.txt"),
         PtestParamFile(sync.absolutePath("PConvparam_test.xml")),
         //inConvergenceMap(sync.absolutePath("data/SphKappa.fits"))
         inConvergenceMap(sync.absolutePath("data/SphereKE.fits"))
 {
    sync.download();
  }
};

//-----------------------------------------------------------------------------

//BOOST_AUTO_TEST_SUITE (SphericalPeakCount_test)
BOOST_FIXTURE_TEST_SUITE (SphericalPeakCount_test, SphericalPeakDataSyncFixture)
//-----------------------------------------------------------------------------


BOOST_AUTO_TEST_CASE( overall_test ) {
  Elements::TempDir one;
  test_path = one.path();
  LE3_2D_MASS_WL_PEAK_COUNT::PeakParam Pparams;
  if (true == fileHasField(PtestParamFile, "DpdTwoDMassParamsPeakCatalogConvergence")) {
   std::cout<<"Parameter file is for input Convergence Maps.."<<std::endl;
   Pparams.readPeakConvXML(PtestParamFile.native());
  }
   Healpix_Map< double > mapE;
   read_Healpix_map_from_fits (inConvergenceMap.native(), mapE);
   boost::filesystem::path outputPeakCatalog= "Test_PeakCatalog_conv.fits";

   SphericalPeakCount peakCounter(mapE, Pparams);
   peakCounter.savePeakCatalog_hp((test_path/outputPeakCatalog).native());

    BOOST_CHECK(boost::filesystem::is_regular_file(test_path/outputPeakCatalog));
}


//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
