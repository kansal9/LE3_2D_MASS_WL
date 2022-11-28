/**
 * @file tests/src/CartesianProcessor_test.cpp
 * @date 02/15/22
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
#include "ElementsServices/DataSync.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianProcessor.h"

using namespace ElementsServices::DataSync;
using LE3_2D_MASS_WL_CARTESIAN::CartesianProcessor;

struct CartesianProcessorFixture
{
    DataSync sync;
    path testParamConvergencePatch, testCatFile;
    std::map<std::string, std::string> args_string;
    std::map<std::string, variable_value> args_values;
    double ra_center, dec_center, size, zmin, zmax;
    CartesianProcessor cp;

    CartesianProcessorFixture() :
            sync("LE3_2D_MASS_WL_CARTESIAN/datasync_webdav.conf",
                 "LE3_2D_MASS_WL_CARTESIAN/test_file_list.txt"),
            testParamConvergencePatch(sync.absolutePath("Cparam_test_1patch.xml")),
            testCatFile(sync.absolutePath("InputLE2Catalog.xml"))
    {
        sync.download();
        ra_center = 20;
        dec_center = 20;
        size = 10;
        zmin = 0.1;
        zmax = 3;
        cp.m_shear_cat.fillTest(100, "LENSMC", ra_center - size,
                                            ra_center + size,
                                            dec_center - size,
                                            dec_center + size,
                                            zmin, zmax);
    }
};

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE (CartesianProcessor_test, CartesianProcessorFixture)

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( cpParseOptionsTest, CartesianProcessorFixture ) {
    try
    {
        cp.setOption("paramfile", testParamConvergencePatch.native());
        cp.setOption("workdir", testCatFile.parent_path().native());
        cp.setOption("shear", testCatFile.filename().native());
        cp.parseOptions();
        BOOST_CHECK(true);
    }
    catch(std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        BOOST_THROW_EXCEPTION(e);
    }
}

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( cpProcessPatchTestNoGal, CartesianProcessorFixture ) {
    try
    {
        cp.m_shear_cat.fillTest(0);
        cp.m_params.setParaFileType(parameterType::DpdTwoDMassParamsConvergencePatch);
        // fill basic options for patch processing
        std::vector<PatchDef> patches;
        patches.emplace_back(ra_center * deg2rad, dec_center * deg2rad,
                             size * deg2rad, 0.586 / 60 * deg2rad);
        cp.m_params.setPatches(patches);
        cp.m_params.setNPatches(patches.size());
        cp.m_params.setZMin(std::vector<double>{zmin});
        cp.m_params.setZMax(std::vector<double>{zmax});
        cp.m_params.setNbins(1);
        // process
        cp.process();
        BOOST_CHECK(true);
    }
    catch(std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        BOOST_THROW_EXCEPTION(e);
    }
}

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( cpProcessPatchTest, CartesianProcessorFixture ) {
    try
    {
        cp.m_params.setParaFileType(parameterType::DpdTwoDMassParamsConvergencePatch);
        // fill basic options for patch processing
        std::vector<PatchDef> patches;
        patches.emplace_back(ra_center * deg2rad, dec_center * deg2rad,
                             size * deg2rad, 0.586 / 60 * deg2rad);
        cp.m_params.setPatches(patches);
        cp.m_params.setNPatches(patches.size());
        cp.m_params.setZMin(std::vector<double>{zmin});
        cp.m_params.setZMax(std::vector<double>{zmax});
        cp.m_params.setNbins(1);
        cp.m_params.setRsCorrection(true);
        cp.m_params.setNResamples(2);
        cp.m_params.setNInpaint(2);
        // process
        cp.process();
        BOOST_CHECK(true);
    }
    catch(std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        BOOST_THROW_EXCEPTION(e);
    }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()


