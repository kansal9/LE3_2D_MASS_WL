/**
 * @file tests/src/GenericParam_test.cpp
 * @date 02/04/22
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

#include "LE3_2D_MASS_WL_UTILITIES/GenericParam.h"

using namespace LE3_2D_MASS_WL_UTILITIES;

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE (GenericParam_test)
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( basic_test ) {
    GenericParam param;
    param.print();
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()


