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
    CartesianProcessorFixture() :
            sync("LE3_2D_MASS_WL_CARTESIAN/datasync_webdav.conf",
                 "LE3_2D_MASS_WL_CARTESIAN/test_file_list.txt"),
            testParamConvergencePatch(sync.absolutePath("Cparam_test_1patch.xml")),
            testCatFile(sync.absolutePath("InputLE2Catalog.xml"))
    {
        sync.download();
    }
};

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE (CartesianProcessor_test, CartesianProcessorFixture)

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_CASE( cpConvergencePatchTest, CartesianProcessorFixture ) {
    try
    {
        CartesianProcessor cp;
        cp.setOption("paramfile", testParamConvergencePatch.native());
        cp.setOption("workdir", testCatFile.parent_path().native());
        cp.setOption("shear", testCatFile.filename().native());
        cp.parseOptions();
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


