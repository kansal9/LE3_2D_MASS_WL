/**
 * @file LE3_2D_MASS_WL_CARTESIAN/MassMapping.h
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_MASSMAPPING_H
#define _LE3_2D_MASS_WL_CARTESIAN_MASSMAPPING_H

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GetMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "boost/multi_array.hpp"
#include <utility>

namespace LE3_2D_MASS_WL_CARTESIAN {
/**
 * @class MassMapping
 * @brief
 *
 */
class MassMapping {

public:
  /**
   * @brief Constructor
   */
  MassMapping(LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam);
  /**
   * @brief Destructor
   */
  virtual ~MassMapping() = default;

public:
 /**
 * @brief  Perform Map inversion
 * @param  Input Map in FITS format
 * This method perform mass mapping
 */
 ConvergenceMap* getSheartoConv(const std::string& map);
 /**
 * @brief  Perform inverse Map inversion
 * @param  Input Map in FITS format
 * This method perform mass mapping
 */
 ShearMap* getConvtoShear(const std::string& map);

private:
 ConvergenceMap *m_convMap;
 ShearMap *m_shearMap;
  /**
   *  @brief <m_cartesianParam>, CartesianParam object with catalog parameters
  */
 LE3_2D_MASS_WL_CARTESIAN::CartesianParam m_cartesianParam;

};  // End of MassMapping class

}  // namespace LE3_2D_MASS_WL_CARTESIAN

#endif
