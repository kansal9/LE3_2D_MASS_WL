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
    std::map<std::string, std::string> args;
    CartesianProcessorFixture() :
            sync("LE3_2D_MASS_WL_CARTESIAN/datasync_webdav.conf",
                 "LE3_2D_MASS_WL_CARTESIAN/test_file_list.txt"),
            testParamConvergencePatch(sync.absolutePath("Cparam_test.xml")),
            testCatFile(sync.absolutePath("InputLE2Catalog.xml"))
    {
        sync.download();
    }
};

//-----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE (CartesianProcessor_test, CartesianProcessorFixture)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( convergencePatchTest ) {

    args["paramfile"] = testParamConvergencePatch.native();
    args["workdir"] = testCatFile.parent_path().native();
    args["shear"] = testCatFile.filename().native();

    for(auto iter = args.begin(); iter != args.end(); ++iter)
    {
        std::cout << iter->first << ": " << iter->second << std::endl;
    }

    try
    {
        CartesianProcessor cp(args);
        cp.process();
    }
    catch(ExitCode& exitCode)
    {
        int code = static_cast<int>(exitCode);
        std::cout << "Exit code: " << code << std::endl;
    }

}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()


