/**
 * @file src/lib/Projection.cpp
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

#include "LE3_2D_MASS_WL_CARTESIAN/Projection.h"
static Elements::Logging logger = Elements::Logging::getLogger("Projection");
namespace LE3_2D_MASS_WL_CARTESIAN {

Projection::Projection(){}

std::pair<double, double> Projection::getGnomonicProjection(double ra, double dec, double ra0, double dec0)
{
  ra *= M_PI/180.;
  dec *= M_PI/180.;
  ra0 *= M_PI/180.;
  dec0 *= M_PI/180.;

  double cosc = sin(dec0)*sin(dec) + cos(dec0)*cos(dec)*cos(ra-ra0);
  double x = 1./cosc*cos(dec)*sin(ra-ra0);
  double y = 1./cosc*(cos(dec0)*sin(dec) - sin(dec0)*cos(dec)*cos(ra-ra0));
 // double x = cos(dec)*sin(ra-ra0)*1./cosc;
 // double y = (cos(dec0)*sin(dec) - sin(dec0)*cos(dec)*cos(ra-ra0))*1./cosc;
  //logger.info()<<"x, y : "<<x<<" "<<y;

  return std::pair<double, double> (x, y);
}

std::pair<double, double> Projection::getInverseGnomonicProjection(double x, double y, double ra0, double dec0)
{
  // Handle the case where x = y = 0 when ra = ra0 and dec = dec0
  if (fabs(x)<0.01 && fabs(y)<0.01)
  {
    return std::pair<double, double> (ra0, dec0);
  }
  ra0 *= M_PI/180.;
  dec0 *= M_PI/180.;

  double rho = sqrt(x*x + y*y);
  double c = atan(rho);

  double z1 = 1./rho*(y * sin(c) * cos(dec0));
  double dec = asin(cos(c)*sin(dec0) + z1);

  double z2 = rho * cos(dec0)* cos(c) - y * sin(dec0) *sin(c);
  double ra = ra0 + atan2((x* sin(c)), z2);
 // double dec = asin(cos(c)*sin(dec0) + y*sin(c)*cos(dec0)/rho);
 //double ra = ra0 + atan( x*sin(c)/(rho*cos(dec0)*cos(c) - y*sin(dec0)*sin(c)) );

  dec *= 180./M_PI;
  ra *= 180./M_PI;

  return std::pair<double, double> (ra, dec);
}
}  // namespace LE3_2D_MASS_WL_CARTESIAN
