/**
 * @file LE3_2D_MASS_WL_CARTESIAN/ShearMap.h
 * @date 05/19/20
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_SHEARMAP_H
#define _LE3_2D_MASS_WL_CARTESIAN_SHEARMAP_H

#include "LE3_2D_MASS_WL_CARTESIAN/GetMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"

namespace LE3_2D_MASS_WL_CARTESIAN {
class ConvergenceMap;
/**
 * @class ShearMap
 * @brief
 *
 */
class ShearMap:public GetMap {

public:

  /**
   * @brief Destructor
   */
  virtual ~ShearMap() = default;
 /**
 * @brief Constructor of a Shear Map
 * @param[in] array input values to give to the class
 * @param[in] nGalaxies number of galaxies in the map
 *
 * This constructor builds a generic Map of dimensions sizeXaxis, sizeYaxis, sizeZaxis
 * from the input data array
 */
 ShearMap(double* array, LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam, int nbGalaxies);

 /**
 * @brief Constructor of a Shear Map
 * @param[in] array input values to give to the class
 * @param[in] sizeXaxis number of pixels in the X axis
 * @param[in] sizeYaxis number of pixels in the Y axis
 * @param[in] sizeZaxis number of pixels in the Z axis
 * @param[in] nGalaxies number of galaxies in the map
 * This constructor builds a generic Map of dimensions sizeXaxis, sizeYaxis, sizeZaxis
 * from the input data array
 */
 ShearMap(double* array, int sizeXaxis, int sizeYaxis,
          int sizeZaxis, int nbGalaxies = 0);

 /**
 * @brief Constructor of a Shear Map
 * @param[in] array input values to give to the class
 * @param[in] sizeXaxis number of pixels in the X axis
 * @param[in] sizeYaxis number of pixels in the Y axis
 * @param[in] sizeZaxis number of pixels in the Z axis
 * @param[in] CoordinateBound contains ra, dec and z information of map
 * @param[in] nGalaxies number of galaxies in the map
 * This constructor builds a generic Map of dimensions sizeXaxis, sizeYaxis, sizeZaxis
 * from the input data array
 */
 ShearMap(double* array, int sizeXaxis, int sizeYaxis,
          int sizeZaxis, LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& coordBound, int nbGalaxies = 0);

 /**
 * @brief Constructor of a convergence Map
 * @param[in] filename name of a FITS input file
 * This constructor builds a shear Map from provided FITS input file filename.
 */
 ShearMap(const std::string& filename);

 /**
 * @brief Copy constructor of a Shear Map
 * @param[in] copyMap the Map to copy
 * This copy constructor builds a copy of the input Map
 */
 ShearMap(GetMap const& copyMap);

 /**
 * @brief Returns a ConvergenceMap using K&S algorithm
 * @return a convergence Map corresponding to the input Shear Map
 * This method creates and returns a Convergence Map from the given ShearMap
 * using Kaiser & Squires algorithm
 */
 ConvergenceMap getConvMap();

 /**
 * @brief Computes the reduced shear
 * This method computes the reduced shear from the original shear map and the
 * provided convergence map
 **/
 void computeReducedShear(ConvergenceMap& inputConvMap);

private:

};  // End of ShearMap class

}  // namespace LE3_2D_MASS_WL_CARTESIAN


#endif
