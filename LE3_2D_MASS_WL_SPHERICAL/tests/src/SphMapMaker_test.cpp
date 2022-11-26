
/**
 * @file tests/src/Sph_map_maker_test.cpp
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
#define BOOST_TEST_MODULE Sph_map_maker_test

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "ElementsServices/DataSync.h"
#include "ElementsKernel/Temporary.h"

#include "LE3_2D_MASS_WL_SPHERICAL/SphMapMaker.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalUtils.h"

using LE3_2D_MASS_WL_SPHERICAL::SphMapMaker;
using LE3_2D_MASS_WL_SPHERICAL::SphericalParam;

using namespace Euclid::WeakLensing::TwoDMass::Spherical;
using namespace Euclid::WeakLensing::TwoDMass;
using namespace ElementsServices::DataSync;

boost::filesystem::path test_path;
//-----------------------------------------------------------------------------

struct SphMapMakerFixture
{
    DataSync sync;
    path testParamFile;
    SphericalParam params;
    CatalogData cat;
    SphMapMakerFixture() :
        sync("LE3_2D_MASS_WL_SPHERICAL/datasync_webdav.conf",
                     "LE3_2D_MASS_WL_SPHERICAL/test_file_list.txt"),
        testParamFile(sync.absolutePath("Sparam_test.xml"))
    {
        sync.download();
        cat.fillTest(100);
        readSphericalParameterFile(testParamFile, params, cat);
    }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (Sph_map_maker_test, SphMapMakerFixture)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( sphMapMaker_test )
{
    std::cout << "-- SphMapMaker_Test" << std::endl;
    SphMapMaker mapMaker(params);
    Healpix_Map<double> g1_hmap;
    Healpix_Map<double> g2_hmap;
    Healpix_Map<double> ngal_hp;
    mapMaker.create_ShearMap(cat, g1_hmap, g2_hmap, ngal_hp);
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
