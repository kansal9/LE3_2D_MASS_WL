/**
 * @file tests/src/SphericalInpainting_test.cpp
 * @date 02/01/21
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

#include "ElementsKernel/Logging.h"
#include "ElementsServices/DataSync.h"

#include "LE3_2D_MASS_WL_SPHERICAL/SphericalInpainting.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalUtils.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphMassMapping.h"

#include <ios>
#include <sstream>
#include <iostream>

using LE3_2D_MASS_WL_SPHERICAL::SphMassMapping;
using LE3_2D_MASS_WL_SPHERICAL::SphericalParam;
using LE3_2D_MASS_WL_SPHERICAL::SphericalInpainting;
using Euclid::WeakLensing::TwoDMass::Spherical::SphericalIO;

using namespace Euclid::WeakLensing::TwoDMass::Spherical;
using namespace Euclid::WeakLensing::TwoDMass;
using namespace ElementsServices::DataSync;
using namespace LE3_2D_MASS_WL_UTILITIES;

using boost::filesystem::path;

static Elements::Logging logger = Elements::Logging::getLogger(
        "SphericalInpainting_test");

//-----------------------------------------------------------------------------

struct SphericalInpaintingDataSyncFixture
{
    DataSync sync;
    path testParamFile, gamma;
    SphericalParam params;
    SphericalInpaintingDataSyncFixture() :
            sync("LE3_2D_MASS_WL_SPHERICAL/datasync_webdav.conf",
                 "LE3_2D_MASS_WL_SPHERICAL/test_file_list.txt"),
            testParamFile(sync.absolutePath("Sparam_test.xml")),
            gamma(sync.absolutePath("data/Gamma_32_mask.fits"))
    {
        sync.download();
        CatalogData dummy;
        readSphericalParameterFile(testParamFile, params, dummy);
    }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (SphericalInpainting_test, SphericalInpaintingDataSyncFixture)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( overall_SphericalInpainting_test)
{
    SphericalIO SphericalIO(params);
    auto mapPair = readHealpixMap(gamma.native());

    //get mask
    int nside = mapPair.first.Nside();
    int m_npix = mapPair.first.Npix();
    Healpix_Map<int> Mask;
    Mask.SetNside(nside, RING);
    Mask.fill(1);
    for (int it = 0; it < m_npix; it++)
    {
        if (fabs(mapPair.first[it]) == 0. && fabs(mapPair.second[it]) == 0.)
        {
            Mask[it] = 0;
        }
    }

    // perform direct mass mapping
    SphMassMapping massmapping;
    auto kappaPair = massmapping.create_SheartoConvMap(mapPair.first, mapPair.second);

    //perform Inpainting
    SphericalInpainting Inpainting(mapPair, kappaPair, params);
    kappaPair = Inpainting.performInpainting();

    for (int it = 0; it < m_npix; it++)
    {
        if (Mask[it] == 0)
        {
            BOOST_CHECK(kappaPair.first[it] != 0.);
        }
    }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
