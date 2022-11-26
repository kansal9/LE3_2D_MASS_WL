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
#include <boost/filesystem.hpp>

#include "ElementsServices/DataSync.h"

#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"

#include <ios>
#include <sstream>
#include <iostream>

using LE3_2D_MASS_WL_SPHERICAL::SphericalParam;

using namespace ElementsServices::DataSync;

using boost::filesystem::path;

struct SphericalParamTestEnv
{
    DataSync sync;
    path xmlFileName;
    SphericalParam param;

    SphericalParamTestEnv() :
        sync("LE3_2D_MASS_WL_SPHERICAL/datasync_webdav.conf",
             "LE3_2D_MASS_WL_SPHERICAL/test_file_list.txt"),
        xmlFileName(sync.absolutePath("Sparam_test.xml"))
    {
        sync.download();
        CatalogData dummy;
        param.readConvSphereXMLFile(xmlFileName, dummy);
    }
};
//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (SphericalParam_test, SphericalParamTestEnv)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( overall_test )
{
    std::cout << " --> SphericalParam: Overall_test" << std::endl;

    BOOST_CHECK(param.getRsCorrection());
    BOOST_CHECK_EQUAL(param.getRsThreshold(), 5);
    BOOST_CHECK_CLOSE(param.getRsGaussStd(), 0, 1.0e-07);

    BOOST_CHECK_EQUAL(param.getDenoisingAlgo(), std::string("GaussFilter"));
    BOOST_CHECK_CLOSE(param.getGaussStd(), 0, 1.0e-07);

    BOOST_CHECK_EQUAL(param.getNResamples(), 0);

    BOOST_CHECK_EQUAL(param.getNside(), 16);

    BOOST_CHECK_EQUAL(param.getNbins(), 2);
    BOOST_CHECK_CLOSE(param.getZMin(0), 0, 0.1);
    BOOST_CHECK_CLOSE(param.getZMax(0), 10, 0.1);

    BOOST_CHECK(param.isBalancedBins());

    BOOST_CHECK_EQUAL(param.getNInpaint(), 10);
    BOOST_CHECK(param.isEqualVarPerScale());
    BOOST_CHECK(param.isForceBMode());
    BOOST_CHECK_EQUAL(param.getNInpScale(), 2);
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
