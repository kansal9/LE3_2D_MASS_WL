/**
 * @file tests/src/PeakCatalog_test.cpp
 * @date 07/04/22
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

#include "LE3_2D_MASS_WL_PEAK_COUNT/PeakCatalog.h"

#include "ElementsKernel/Temporary.h"

using LE3_2D_MASS_WL_PEAK_COUNT::PeakCatalog;

using boost::filesystem::path;

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE (PeakCatalog_test)
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( peakCatalog_test )
{
    PeakCatalog peakCat;
    peakCat.addPeak(1, 2, 3, 4, 5, 6);

    BOOST_CHECK_EQUAL(peakCat.getNentries(), 1);

    Elements::TempFile one;
    boost::filesystem::path test_path = one.path();
    peakCat.writePeakCatalog(test_path.native());

    std::cout << "Peak catalog saved to: " << test_path << std::endl;
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()

