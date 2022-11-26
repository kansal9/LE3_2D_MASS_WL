/**
 * @file tests/src/NoisyCatalogData_test.cpp
 * @date 10/13/20
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

#include "ElementsKernel/Logging.h"
#include "LE3_2D_MASS_WL_UTILITIES/NoisyCatalogData.h"
#include <boost/test/unit_test.hpp>
#include <functional>
#include <fstream>
#include <string>
#include <vector>

using namespace LE3_2D_MASS_WL_UTILITIES;

static Elements::Logging logger = Elements::Logging::getLogger(
        "Radomization_test");

//-----------------------------------------------------------------------------

struct NoisyCatalogDataFixture
{
    CatalogData cat;
    NoisyCatalogDataFixture()
    {
        cat.fillTest(10);
    }
};

//-----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE (NoisyCatalogData_test, NoisyCatalogDataFixture)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( createNoisyData_test )
{
    logger.info() << "-- NoisyCatalogData: createNoisyData";

    CatalogData noisyCat = NoisyCatalogData::getNoisyCatalog(cat);

    BOOST_CHECK(cat.getNentries() == noisyCat.getNentries());
    BOOST_CHECK(cat["ra"](0) == noisyCat["ra"](0));
    BOOST_CHECK_PREDICATE(std::not_equal_to<double>(),
            (cat["g1"](1))(noisyCat["g1"](1)));
    BOOST_CHECK_PREDICATE(std::not_equal_to<double>(),
            (cat["g2"](1))(noisyCat["g2"](1)));
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
