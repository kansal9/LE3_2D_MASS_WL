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
 * @file tests/src/CartesianAlgoKS_test.cpp
 * @date 10/21/19
 * @author user
 */

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

#include "ElementsKernel/Temporary.h"
#include "ElementsServices/DataSync.h"

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "LE3_2D_MASS_WL_UTILITIES/PatchDef.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include <ios>
#include <sstream>
#include <iostream>
#include <exception>

using namespace LE3_2D_MASS_WL_CARTESIAN;
using namespace LE3_2D_MASS_WL_UTILITIES;
using namespace ElementsServices::DataSync;

namespace fs = boost::filesystem;

// handle on created path names
fs::path test_path;
//-----------------------------------------------------------------------------

struct CartesianAlgoKSDataSyncFixture
{
    DataSync sync;
    path testParamFile, testCatFile;
    CartesianAlgoKSDataSyncFixture() :
            sync("LE3_2D_MASS_WL_CARTESIAN/datasync_webdav.conf",
                 "LE3_2D_MASS_WL_CARTESIAN/test_file_list.txt"),
            testParamFile(sync.absolutePath("Cparam_test_CartesianAlgoKS.xml")),
            testCatFile(sync.absolutePath("data/lensmcCat.fits"))
    {
        sync.download();
    }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (CartesianAlgoKS_test, CartesianAlgoKSDataSyncFixture)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( overallCartesianAlgo_test )
{
    std::cout << "OverallCartesianAlgo_test" << std::endl;

    // read catalog
    CatalogData cat;
    fs::path workdir = testCatFile.parent_path().parent_path();
    fs::path cataloFileName = testCatFile.filename();
    cat.getCatalogData(workdir, cataloFileName);

    // get redshift extrema
    double zMin, zMax;
    vecMinMax(cat["z"], zMin, zMax);

    // directories
    Elements::TempDir one;
    test_path = one.path();
    auto datadir = test_path / "data";
    fs::create_directories(datadir);

    CartesianParam params;
    // Case 1: check file is of XML format
    if (checkFileType(testParamFile, signXML))
    {
        // read XML, get parameter file type
        // check if it's to get a patch
        if (fileHasField(testParamFile, "DpdTwoDMassParamsConvergencePatch"))
        {
            params.readConvPatchXMLFile(testParamFile.native(), cat);
        }
    }

    CartesianAlgoKS CartesainAlgo(params);
    for (size_t it = 0; it < params.getNPatches(); it++)
    {
        std::cout << "Patch number: " << it << std::endl;

        const PatchDef& patch = params.getPatches()[it];
        double ramin = patch.getRaMin();
        double decmin = patch.getDecMin();
        double ramax = patch.getRaMax();
        double decmax = patch.getDecMax();
        // PatchDef CB(ramin, ramax, decmin, decmax, zMin, zMax);

        std::cout << "ramin: " << ramin * rad2deg << std::endl;
        std::cout << "decmin: " << decmin * rad2deg << std::endl;
        std::cout << "ramax: " << ramax * rad2deg << std::endl;
        std::cout << "decmax: " << decmax * rad2deg << std::endl;


        std::cout << "ExtractShearMap_test" << std::endl;
        fs::path shearMapPath = "ExtractedShearMap.fits";
        CartesainAlgo.extractShearMap(datadir / shearMapPath, cat, patch);

        BOOST_CHECK(fs::exists((datadir / shearMapPath).native()));

        ShearMap inputShearMap(datadir / shearMapPath);
        ShearMap outputShearMap(inputShearMap, false);
        ConvergenceMap convergenceMap(inputShearMap, false);

        // copy ngal to output maps
        convergenceMap.singleAxisCopy(inputShearMap, 2);
        outputShearMap.singleAxisCopy(inputShearMap, 2);

        std::cout << "PerformKSMassMapping_test" << std::endl;
        CartesainAlgo.performMassMapping(inputShearMap, convergenceMap,
                test_path);

        std::cout << "PerformInverseKSMassMapping_test" << std::endl;
        CartesainAlgo.performInverseKSMassMapping(convergenceMap,
                outputShearMap);

        for (size_t i = 0; i < inputShearMap.getXdim(); i++)
        {
            for (size_t j = 0; j < inputShearMap.getYdim(); j++)
            {
                // check number of galaxies
                BOOST_CHECK(
                        inputShearMap.getBinValue(i, j, 2)
                                == outputShearMap.getBinValue(i, j, 2));

                // check gamma1 and gamma2
                if (inputShearMap.getBinValue(i, j, 2) > 0)
                {
                    BOOST_CHECK(
                            fabs(
                                    inputShearMap.getBinValue(i, j, 0)
                                            - outputShearMap.getBinValue(i, j,
                                                    0)) < 1e-6);
                    BOOST_CHECK(
                            fabs(
                                    inputShearMap.getBinValue(i, j, 1)
                                            - outputShearMap.getBinValue(i, j,
                                                    1)) < 1e-6);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
