/**
 * @file WL_2D_SPHERICAL_MASS_MAPPING/Sph_map_maker.h
 * @date 03/21/19
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

#ifndef _WL_2D_SPHERICAL_MASS_MAPPING_SPH_MAP_MAKER_H
#define _WL_2D_SPHERICAL_MASS_MAPPING_SPH_MAP_MAKER_H

#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"
#include <vector>
#include <memory>

namespace  LE3_2D_MASS_WL_SPHERICAL {

/**
 * @class Sph_map_maker
 * @brief
 *
 */
class Sph_map_maker {

public:
  /**
   * @brief Constructor
   */
  Sph_map_maker(LE3_2D_MASS_WL_SPHERICAL::SphericalParam &SphParam);

  /**
   * @brief Destructor
   */
  virtual ~Sph_map_maker() = default;

  /**
  @brief  This method create Shear Map using input catalog data
  @param  <catData> input catalog data
  @retun  <Shear1, Shear2, GalCountMap> Shear healpix maps
  */
  //std::pair<Healpix_Map<double>, Healpix_Map<double> > create_ShearMap(std::vector<std::vector<double> > &catData);
  std::tuple<Healpix_Map<double>, Healpix_Map<double>, Healpix_Map<double> > create_ShearMap
                                                                (std::vector<std::vector<double> > &catData);
private:

  /**
   *  @brief <hb>, Healpix_Base object
  */
 Healpix_Base hb;

  /**
  @brief  This is the number of pixels in healpix map
  */
 unsigned int npix;

  /**
   *  @brief <m_SphParam>, SphericalParam object with input parameters
  */
 LE3_2D_MASS_WL_SPHERICAL::SphericalParam m_SphParam;

};  // End of Sph_map_maker class
}  // namespace  LE3_2D_MASS_WL_SPHERICAL
#endif
