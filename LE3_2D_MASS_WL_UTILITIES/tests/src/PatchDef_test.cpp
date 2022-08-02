/**
 * @file tests/src/PatchDef_test.cpp
 * @date 02/25/20
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

#include "LE3_2D_MASS_WL_UTILITIES/PatchDef.h"

#include <boost/test/unit_test.hpp>

using LE3_2D_MASS_WL_UTILITIES::PatchDef;

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE (PatchDef_test)
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( overall_test )
{
    double raCtr(20);
    double decCtr(30);
    double width(10);
    double size(0.586 / 60);
    PatchDef patch(raCtr, decCtr, width, size);

    BOOST_CHECK_CLOSE(patch.getRaMin(), raCtr - width/2., 0.1);
    BOOST_CHECK_CLOSE(patch.getRaMax(), raCtr + width/2., 0.1);
    BOOST_CHECK_CLOSE(patch.getDecMin(), decCtr - width/2., 0.1);
    BOOST_CHECK_CLOSE(patch.getDecMax(), decCtr + width/2., 0.1);
    BOOST_CHECK_CLOSE(patch.getPatchWidth(), width, 0.1);
    BOOST_CHECK_CLOSE(patch.getPixelSize(), size, 0.1);
    BOOST_CHECK_EQUAL(patch.getXbin(), 1024);
    BOOST_CHECK_EQUAL(patch.getYbin(), 1024);
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
