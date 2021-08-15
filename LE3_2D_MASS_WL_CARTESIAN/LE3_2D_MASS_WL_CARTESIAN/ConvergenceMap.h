/**
 * @file LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_CONVERGENCEMAP_H
#define _LE3_2D_MASS_WL_CARTESIAN_CONVERGENCEMAP_H

#include "LE3_2D_MASS_WL_CARTESIAN/GetMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"

namespace LE3_2D_MASS_WL_CARTESIAN {
class ShearMap;
/**
 * @class ConvergenceMap
 * @brief
 *
 */
class ConvergenceMap:public GetMap {

public:

  /**
   * @brief Destructor
   */
  virtual ~ConvergenceMap() = default;
 /**
 * @brief Constructor of a ConvergenceMap
 * @param[in] array input values to give to the class
 * @param[in] nGalaxies number of galaxies in the map
 *
 * This constructor builds a generic Map of dimensions sizeXaxis, sizeYaxis, sizeZaxis
 * from the input data array
 */
 ConvergenceMap(double* array, LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam, int nbGalaxies);

 /**
 * @brief Constructor of a Convergence Map
 * @param[in] array input values to give to the class
 * @param[in] sizeXaxis number of pixels in the X axis
 * @param[in] sizeYaxis number of pixels in the Y axis
 * @param[in] sizeZaxis number of pixels in the Z axis
 * @param[in] nGalaxies number of galaxies in the map
 * This constructor builds a generic Map of dimensions sizeXaxis, sizeYaxis, sizeZaxis
 * from the input data array
 */
 ConvergenceMap(double* array, int sizeXaxis, int sizeYaxis,
                int sizeZaxis, int nbGalaxies = 0);

 /**
 * @brief Constructor of a Convergence Map
 * @param[in] array input values to give to the class
 * @param[in] sizeXaxis number of pixels in the X axis
 * @param[in] sizeYaxis number of pixels in the Y axis
 * @param[in] sizeZaxis number of pixels in the Z axis
 * @param[in] CoordinateBound contains ra, dec and z information of map
 * @param[in] nGalaxies number of galaxies in the map
 * This constructor builds a generic Map of dimensions sizeXaxis, sizeYaxis, sizeZaxis
 * from the input data array
 */
 ConvergenceMap(double* array, int sizeXaxis, int sizeYaxis, int sizeZaxis,
              LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& coordBound, int nbGalaxies = 0);

 /**
 * @brief Constructor of a Map
 * @param[in] filename name of a FITS input file
 * This constructor builds a convergence Map from provided FITS input file filename.
 * The provided FITS file should have the image as Primary Header in order to be
 * read properly.
 */
 ConvergenceMap(const std::string& filename);

 /**
 * @brief Copy constructor of a Convergence Map
 * @param[in] copyMap the Map to copy
 * This copy constructor builds a copy of the input Map
 */
 ConvergenceMap(GetMap const& copyMap);

 /**
 * @brief Returns a ShearMap using K&S algorithm
 * @return a ShearMap corresponding to the input ConvergenceMap
 * This method creates and returns a ShearMap from the given ConvergenceMap
 * using Kaiser & Squires algorithm
 */
 ShearMap getShearMap();

private:

};  // End of ConvergenceMap class

}  // namespace LE3_2D_MASS_WL_CARTESIAN


#endif
