/**
 * @file LE3_2D_MASS_WL_CARTESIAN/MapMaker.h
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_MAPMAKER_H
#define _LE3_2D_MASS_WL_CARTESIAN_MAPMAKER_H

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GetMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include <utility>
#include <vector>

namespace LE3_2D_MASS_WL_CARTESIAN {
/**
 * @class MapMaker
 * @brief
 *
 */
class MapMaker {

public:

  /**
   * @brief Destructor
   */
  virtual ~MapMaker() = default;

public:
 /**
  * @brief constructor
  */
 MapMaker(std::vector<std::vector<double> >& inData,
          LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam);

 /**
 * @brief  extract Shear Map
 * @param  Catalog Data
 * @param  respective Parameters from Parameter file
 * This method extract the Shear Map from catalog
 * using parameters like min Ra, min Dec, max Ra, max Dec, Pixels values
 * or using parameters like mapSize, mapCenter, pixelSize
 */

 ShearMap* getShearMap(LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& CB);
 /**
 * @brief  extract Convergence Map
 * @param  Catalog Data
 * @param  respective Parameters from Parameter file
 * This method extract the Convergence Map from catalog
 * using parameters like min Ra, min Dec, max Ra, max Dec, Pixels values
 * or using parameters like mapSize, mapCenter, pixelSize
 */

 ConvergenceMap* getConvMap(LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& CB);

private:

 /**
 * @brief   extract Map (either shear or convergence)
 * @param   mapType i.e. shear or Convergece
 * @param   Catalog Data
 * @param   respective Parameters from Parameter file
 * @return  it returns number of galaxies in the map and array of map
 * This method extract the Map from catalog
 * using parameters like min Ra, min Dec, max Ra, max Dec, Pixels values
 * or using parameters like mapSize, mapCenter, pixelSize
 */
std::pair<long, double*> getMap(const Euclid::WeakLensing::TwoDMass::mapType type,
                                              LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& CB);
  /**
   *  @brief <m_inData>, Catalog Data
  */
std::vector<std::vector<double> > inputData;
  /**
   *  @brief <cartesianParam>, CartesianParam object with catalog parameters
  */
LE3_2D_MASS_WL_CARTESIAN::CartesianParam cartesianParam;

};  // End of MapMaker class

}  // namespace LE3_2D_MASS_WL_CARTESIAN

#endif
