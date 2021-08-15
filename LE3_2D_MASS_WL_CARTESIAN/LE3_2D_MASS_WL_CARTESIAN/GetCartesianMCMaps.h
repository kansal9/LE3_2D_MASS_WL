/**
 * @file LE3_2D_MASS_WL_CARTESIAN/GetCartesianMCMaps.h
 * @date 10/13/20
 * @author vkansal
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_GETCARTESIANMCMAPS_H
#define _LE3_2D_MASS_WL_CARTESIAN_GETCARTESIANMCMAPS_H

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GetMap.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "boost/multi_array.hpp"
#include <utility>
#include <vector>

namespace LE3_2D_MASS_WL_CARTESIAN {

/**
 * @class GetCartesianMCMaps
 * @brief
 *
 */
class GetCartesianMCMaps {

public:

  /**
   * @brief Constructor
   */
  GetCartesianMCMaps (std::vector<std::vector<double> >& inputData,
          LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam, LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& CB);

  /**
   * @brief Destructor
   */
  virtual ~GetCartesianMCMaps() = default;

  /**
   * @brief get the Denoised Shear Map
   * @param[in] None
   * @return the denoised shear Map
   */
  LE3_2D_MASS_WL_CARTESIAN::ShearMap* getDeNoisedShearMap();

  /**
   * @brief get the Noised Shear Map after randomising the input Shear
   * @param[in] None
   * @return the a vector containing N noised shear Map
   */
  std::vector<LE3_2D_MASS_WL_CARTESIAN::ShearMap*> getNoisedShearMaps();

  /**
   * @brief perform addition
   * @param[in] DenoisedShearMap
   * @param[in] NoisedShearMap
   * @param[in] filename name of output FITS ShearMap
   * @return true/false
   */
  bool performAddition(LE3_2D_MASS_WL_CARTESIAN::ShearMap& DenoisedShearMap,
                                     LE3_2D_MASS_WL_CARTESIAN::ShearMap& NoisedShearMap, const std::string& ShearMap);

private:
  /**
   *  @brief <m_inData>, Catalog Data
  */
std::vector<std::vector<double> > m_inData;
  /**
   *  @brief <m_cartesianParam>, CartesianParam object with catalog parameters
  */
LE3_2D_MASS_WL_CARTESIAN::CartesianParam m_cartesianParam;
  /**
   *  @brief <m_CB>, CoordinateBound object
  */
LE3_2D_MASS_WL_CARTESIAN::CoordinateBound m_CB;

};  // End of GetCartesianMCMaps class

}  // namespace LE3_2D_MASS_WL_CARTESIAN

#endif
