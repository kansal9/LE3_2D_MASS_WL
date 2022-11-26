/**
 * @file tests/src/Projection_test.cpp
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

#include <boost/test/unit_test.hpp>

#include "ElementsKernel/Logging.h"

#include "LE3_2D_MASS_WL_CARTESIAN/Projection.h"

using LE3_2D_MASS_WL_CARTESIAN::Projection;

static Elements::Logging logger = Elements::Logging::getLogger(
        "Projection_test");

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE (Projection_test)
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( example_test )
{
    logger.info() << "--Projection: example_test";
    //Loop over right ascensions
    for (int i = 0; i < 36; i++)
    {
        //Loop over declinations
        for (int j = 0; j < 18; j++)
        {
            double ra0 = 10. * i;
            double dec0 = 10. * j;
            double ra = ra0;
            double dec = dec0;
            Projection Projection;
            std::pair<double, double> myXY = Projection.getGnomonicProjection(
                    ra, dec, ra0, dec0);
            std::pair<double, double> myBackRaDec =
                    Projection.getInverseGnomonicProjection(myXY.first,
                            myXY.second, ra0, dec0);

            BOOST_CHECK_CLOSE(ra, myBackRaDec.first, 0.0001);
            BOOST_CHECK_CLOSE(dec, myBackRaDec.second, 0.0001);
        }
    }
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( example_test2 )
{
    logger.info() << "--Projection: example_test2";
    double raMin = 15 * M_PI / 180.;
    double raMax = 25 * M_PI / 180.;
    double decMin = 15 * M_PI / 180.;
    double decMax = 25 * M_PI / 180.;
    double ra0 = 0.5 * (raMax + raMin);
    double dec0 = 0.5 * (decMax + decMin);

    Projection Projection;
    std::pair<double, double> myXYMin = Projection.getGnomonicProjection(raMin,
            decMin, ra0, dec0);
    std::pair<double, double> myBackRaDecMin =
            Projection.getInverseGnomonicProjection(myXYMin.first,
                    myXYMin.second, ra0, dec0);
    std::pair<double, double> myXYMax = Projection.getGnomonicProjection(raMax,
            decMax, ra0, dec0);
    std::pair<double, double> myBackRaDecMax =
            Projection.getInverseGnomonicProjection(myXYMax.first,
                    myXYMax.second, ra0, dec0);
    BOOST_CHECK_CLOSE(raMin, myBackRaDecMin.first, 0.0001);
    BOOST_CHECK_CLOSE(decMin, myBackRaDecMin.second, 0.0001);
    BOOST_CHECK_CLOSE(raMax, myBackRaDecMax.first, 0.0001);
    BOOST_CHECK_CLOSE(decMax, myBackRaDecMax.second, 0.0001);
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
