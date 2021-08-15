/**
 * @file LE3_2D_MASS_WL_CARTESIAN/Projection.h
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_PROJECTION_H
#define _LE3_2D_MASS_WL_CARTESIAN_PROJECTION_H
#include<utility>
#include<iostream>
#include<cmath>
#include<math.h>
#include "ElementsKernel/Logging.h"
namespace LE3_2D_MASS_WL_CARTESIAN {

/**
 * @class Projection
 * @brief
 *
 */
class Projection {

public:

  /**
   * @brief Constructor
   */
  Projection();
  /**
   * @brief Destructor
   */
  virtual ~Projection() = default;

  /**
   * @brief     Gets the gnomonic projection of given ra and dec
   * @param[in] ra the input right ascension in degrees
   * @param[in] dec the input declination in degrees
   * @param[in] ra0 the input right ascension on which to perform the projection in degrees
   * @param[in] dec0 the input declination on which to perform the projection in degrees
   *
   * @return    a pair containing the projected X and Y values
   *
   * @details   This method returns the gnomonic projection of ra and dec for a given ra0 and dec0 on
   *            which to perform the projection. The return value is a pair of doubles containing X and Y projected.
   */
std::pair<double, double> getGnomonicProjection(double ra, double dec, double ra0, double dec0);

  /**
   * @brief     Gets the inverse gnomonic projection of given X and Y
   * @param[in] x the input x
   * @param[in] y the input y
   * @param[in] ra0 the input right ascension on which to perform the projection in degrees
   * @param[in] dec0 the input declination on which to perform the projection in degrees
   *
   * @return    a pair containing the backprojected ra and dec values in degrees
   *
   * @details   This method returns the gnomonic backprojection of x and y for a given ra0 and dec0 on
   *            which to perform the backprojection. The return value is a pair of doubles containing
   *            the right ascension and declination backprojected.
   */
std::pair<double, double> getInverseGnomonicProjection(double x, double y, double ra0, double dec0);

private:

};  // End of Projection class
}  // namespace LE3_2D_MASS_WL_CARTESIAN
#endif
