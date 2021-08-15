/**
 * @file LE3_2D_MASS_WL_CARTESIAN/CoordinateBound.h
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_COORDINATEBOUND_H
#define _LE3_2D_MASS_WL_CARTESIAN_COORDINATEBOUND_H

namespace LE3_2D_MASS_WL_CARTESIAN {

/**
 * @class CoordinateBound
 * @brief
 *
 */
class CoordinateBound {

public:

  /**
   * @brief Destructor
   */
  virtual ~CoordinateBound() = default;

  /**
   * @brief Constructor
   * @param[in] raMin the minimum right ascension
   * @param[in] raMax the maximum right ascension
   * @param[in] decMin the minimum declination
   * @param[in] decMax the maximum declination
   * @param[in] zMin the minimum redshift
   * @param[in] zMax the maximum redshift
   *
   */
  CoordinateBound(double raMin, double raMax, double decMin, double decMax, double zMin, double zMax);

  /**
   * @brief Returns the value of the minimum right ascension
   */
  double getRaMin() const;

  /**
   * @brief Returns the value of the maximum right ascension
   */
  double getRaMax() const;

  /**
    * @brief Returns the value of the minimum declination
    */
  double getDecMin() const;

  /**
    * @brief Returns the value of the maximum declination
    */
  double getDecMax() const;

  /**
    * @brief Returns the value of the minimum redshift
    */
  double getZMin() const;

  /**
    * @brief Returns the value of the maximum redshift
    */
  double getZMax() const;

private:
  double m_raMin, m_raMax, m_decMin, m_decMax, m_zMin, m_zMax;

};  // End of CoordinateBound class

}  // namespace LE3_2D_MASS_WL_CARTESIAN

#endif
