
/**
 * @file tests/src/Sph_mass_mapping_test.cpp
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
#define BOOST_TEST_MODULE Sph_mass_mapping_test

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "ElementsServices/DataSync.h"
#include "ElementsKernel/Temporary.h"

#include "LE3_2D_MASS_WL_SPHERICAL/SphMassMapping.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphMapMaker.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalUtils.h"

#include <iostream>

using LE3_2D_MASS_WL_SPHERICAL::SphMapMaker;
using LE3_2D_MASS_WL_SPHERICAL::SphMassMapping;
using LE3_2D_MASS_WL_SPHERICAL::SphericalParam;
using Euclid::WeakLensing::TwoDMass::Spherical::SphericalIO;

using namespace ElementsServices::DataSync;
using namespace Euclid::WeakLensing::TwoDMass::Spherical;
using namespace Euclid::WeakLensing::TwoDMass;

boost::filesystem::path test_path;

//-----------------------------------------------------------------------------
struct SphMassMappingFixture
{
    int nside;
    DataSync sync;
    path testParamFile;
    SphericalParam params;
    CatalogData cat;
    SphMassMappingFixture() : nside(0),
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
BOOST_FIXTURE_TEST_SUITE (Sph_mass_mapping_test, SphMassMappingFixture)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( massMappingHealpix_test )
{
    std::cout << "-- SphMassMapping: Test" << std::endl;
    SphMapMaker mapMaker(params);
    SphericalIO SphericalIO(params);
    Healpix_Map<double> g1_hmap;
    Healpix_Map<double> g2_hmap;
    Healpix_Map<double> ngal_hp;
    mapMaker.create_ShearMap(cat, g1_hmap, g2_hmap, ngal_hp);
    using Elements::TempDir;
    // block creation for local variables
    TempDir six;
    test_path = six.path();

    boost::filesystem::path Shearfilename = test_path / "GammaMap.fits";
    bool status = SphericalIO.write_Map(Shearfilename.native(), g1_hmap,
            "GAMMA1");
    BOOST_REQUIRE(status = true);
    status = SphericalIO.write_Map(Shearfilename.native(), g2_hmap, "GAMMA2");
    BOOST_REQUIRE(status = true);

    if (g1_hmap.Scheme() == RING)
    {
        g1_hmap.swap_scheme();
    }
    SphMassMapping massmapping;
    massmapping.correctScheme(&g1_hmap);
    BOOST_CHECK(g1_hmap.Scheme() == RING);
    if (g1_hmap.Scheme() == RING)
    {
        g1_hmap.swap_scheme();
    }
    if (g2_hmap.Scheme() == RING)
    {
        g2_hmap.swap_scheme();
    }
    std::pair<Healpix_Map<double>, Healpix_Map<double> > convPair =
            massmapping.create_SheartoConvMap(g1_hmap, g2_hmap);

    nside = massmapping.getNside();
    BOOST_CHECK(nside == params.getNside());
    fs::path Convfilename = test_path / "KappaMap.fits";
    status = SphericalIO.write_Map(Convfilename.native(), convPair.second,
            "KAPPA_E");
    BOOST_REQUIRE(status = true);
    status = SphericalIO.write_Map(Convfilename.native(), convPair.first,
            "KAPPA_B");
    BOOST_REQUIRE(status = true);

    std::pair<Healpix_Map<double>, Healpix_Map<double> > shearPair_mm =
            massmapping.create_ConvtoShearMap(convPair.first, convPair.second);

    boost::filesystem::path Shearfile2 = test_path / "GammaMassMap.fits";
    status = SphericalIO.write_Map(Shearfile2.native(), shearPair_mm.second,
            "GAMMA1");
    BOOST_REQUIRE(status = true);
    status = SphericalIO.write_Map(Shearfile2.native(), shearPair_mm.first,
            "GAMMA2");
    BOOST_REQUIRE(status = true);
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( reduceShearHealpix_test )
{
    std::cout << "-- SphMassMapping: reduceShearHealpix_Test" << std::endl;
    SphMapMaker mapMaker(params);
    Healpix_Map<double> g1_hmap;
    Healpix_Map<double> g2_hmap;
    Healpix_Map<double> ngal_hp;
    mapMaker.create_ShearMap(cat, g1_hmap, g2_hmap, ngal_hp);
    SphMassMapping massmapping;
    std::pair<Healpix_Map<double>, Healpix_Map<double> > convPair =
            massmapping.create_SheartoConvMap(g1_hmap, g2_hmap);
    std::pair<Healpix_Map<double>, Healpix_Map<double> > shearPair =
            std::make_pair(g1_hmap, g2_hmap);
    for (int it = 0; it < REDUCESHEARNITER; it++)
    {
        massmapping.computeReducedShear_hp(shearPair, convPair);
    }
}
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END ()
