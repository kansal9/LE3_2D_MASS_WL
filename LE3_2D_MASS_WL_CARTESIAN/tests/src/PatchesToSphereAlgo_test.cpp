/**
 * @file tests/src/PatchesToSphereAlgo_test.cpp
 * @date 08/17/21
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
#include <boost/filesystem/fstream.hpp>

#include "ElementsKernel/Temporary.h"
#include "ElementsServices/DataSync.h"

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "LE3_2D_MASS_WL_CARTESIAN/PatchesToSphereAlgo.h"

#include <iostream>
#include <cstdio>

using namespace LE3_2D_MASS_WL_CARTESIAN;
using namespace LE3_2D_MASS_WL_UTILITIES;
using namespace ElementsServices::DataSync;
namespace fs = boost::filesystem;

// handle on created path names
boost::filesystem::path test_path;

struct PatchesToSphereDataSyncFixture
{
    DataSync sync; // This is the synchronizer
    // These are just shortcuts
    path testParamFile, testFile;
    PatchesToSphereDataSyncFixture() :
            sync(	// Here is the connection configuration file
                    "LE3_2D_MASS_WL_CARTESIAN/datasync_webdav.conf",
                    // Here is the dependency configuration file
                    "LE3_2D_MASS_WL_CARTESIAN/test_file_list.txt"), testParamFile(
                    sync.absolutePath("ParamsConvPatchesToSphere.xml")), testFile(
                    sync.absolutePath("ShearMapsList_Test.json"))
    {
        sync.download();
    }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (PatchesToSphereAlgo_test, PatchesToSphereDataSyncFixture)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( patchesToSphereAlgoOverall_test )
{
    std::cout << "PatchesToSphereAlgoOverall_test" << std::endl;
    CartesianParam params;
    CatalogData dummy;

    auto workdir = testFile.branch_path();
    auto datadir = workdir / "data";

    CartesianAlgoKS CartesainAlgo(params);
    if (fileHasField(testParamFile,
            "DpdTwoDMassParamsConvergencePatchesToSphere"))
    {
        std::cout << "Parameter file is for Convergence Sphere.." << std::endl;
        readParameterFile(testParamFile, params, dummy);
    }

    PatchesToSphereAlgo p2s(params);
    fs::path outputPatchesCenters = fs::path("outputPatchesCenters.fits");
    fs::path outConvergenceMaps = fs::path(
            "NoisySphericalConvergenceMapsList.json");

    std::vector<std::vector<double> > data;
    fs::path inputShearMaps = testFile.filename();
    std::cout << workdir << std::endl;
    std::cout << inputShearMaps << std::endl;
    p2s.getNoisyConvergenceMaps(workdir, inputShearMaps, outConvergenceMaps);

    p2s.getPatchData(outConvergenceMaps, workdir,
            (datadir / outputPatchesCenters).string(), data);

    Healpix_Map<double> mapE, mapB, mapGalCount;
    p2s.getHealpixFormatMap(data, mapE, mapB, mapGalCount);

    // save maps in temp dir
    Elements::TempDir one;
    datadir = one.path() / "data";
    fs::create_directories(datadir);

    fs::path outputHealpixConvergence = fs::path(
            "EUC_LE3_WL_PatchesToSphere_NoisyConvergenceMap.fits");
    fs::path GalCountMap = fs::path(
            "EUC_LE3_WL_PatchesToSphere_GalCount_NSide.fits");

    p2s.write_Map((datadir / outputHealpixConvergence).native(), mapE,
            "KAPPA_E");
    p2s.write_Map((datadir / outputHealpixConvergence).native(), mapB,
            "KAPPA_B");
    p2s.write_Map((datadir / GalCountMap).native(), mapGalCount, "GALCOUNT");

    BOOST_CHECK(fs::exists((datadir / outputHealpixConvergence).native()));
    BOOST_CHECK(fs::exists((datadir / GalCountMap).native()));
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()

