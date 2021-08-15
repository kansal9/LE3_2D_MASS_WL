/**
 * @file tests/src/CoordinateBound_test.cpp
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

#include <boost/test/unit_test.hpp>

#include "LE3_2D_MASS_WL_CARTESIAN/CoordinateBound.h"
using LE3_2D_MASS_WL_CARTESIAN::CoordinateBound;
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE (CoordinateBound_test)

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( Overall_test ) {

  double raMin(21.3);
  double raMax(10.);
  double decMin(20.2);
  double decMax(12.9);
  double zMin(0.2);
  double zMax(2.2);
  CoordinateBound cb(raMin, raMax, decMin, decMax, zMin, zMax);

  BOOST_CHECK_CLOSE(cb.getRaMin(), raMin, 0.1);
  BOOST_CHECK_CLOSE(cb.getRaMax(), raMax, 0.1);
  BOOST_CHECK_CLOSE(cb.getDecMin(), decMin, 0.1);
  BOOST_CHECK_CLOSE(cb.getDecMax(), decMax, 0.1);
  BOOST_CHECK_CLOSE(cb.getZMin(), zMin, 0.1);
  BOOST_CHECK_CLOSE(cb.getZMax(), zMax, 0.1);

}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
