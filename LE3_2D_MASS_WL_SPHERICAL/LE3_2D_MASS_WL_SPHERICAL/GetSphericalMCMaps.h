/**
 * @file LE3_2D_MASS_WL_SPHERICAL/GetSphericalMCMaps.h
 * @date 10/13/20
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

#ifndef _LE3_2D_MASS_WL_SPHERICAL_GETSPHERICALMCMAPS_H
#define _LE3_2D_MASS_WL_SPHERICAL_GETSPHERICALMCMAPS_H

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_SPHERICAL/Sph_map_maker.h"
#include "LE3_2D_MASS_WL_SPHERICAL/Sph_mass_mapping.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalUtils.h"

#include <utility>
#include <vector>
#include <math.h>

namespace LE3_2D_MASS_WL_SPHERICAL {

/**
 * @class GetSphericalMCMaps
 * @brief
 *
 */
class GetSphericalMCMaps {

public:

  /**
   * @brief Constructor
   */
  GetSphericalMCMaps (std::vector<std::vector<double> >& inputData);

  /**
   * @brief Constructor
   */
  GetSphericalMCMaps (std::vector<std::vector<double> >& inputData,
                     LE3_2D_MASS_WL_SPHERICAL::SphericalParam &SphParam);

  /**
   * @brief Destructor
   */
  virtual ~GetSphericalMCMaps() = default;

  /**
   * @brief get the Denoised Shear Map
   * @param[in] None
   * @return the denoised shear Map
   */
  std::pair<Healpix_Map<double>, Healpix_Map<double> > getDeNoisedShearMap();

  /**
   * @brief get the Noised Shear Map after randomising the input Shear
   * @param[in] None
   * @return the a vector containing N noised shear Map
   */
  std::vector<std::pair<Healpix_Map<double>, Healpix_Map<double> > > getNoisedShearMaps();

  /**
   * @brief perform addition
   * @param[in] DenoisedShearMap
   * @param[in] NoisedShearMap
   * @param[in] filename name of output FITS ShearMap
   * @return true/false
   */
  bool performAddition(std::pair<Healpix_Map<double>, Healpix_Map<double> >& DenoisedShearMap,
               std::pair<Healpix_Map<double>, Healpix_Map<double> >& NoisedShearMap, const std::string& ShearMap);

private:

  /**
   *  @brief <m_inData>, Catalog Data
  */
 std::vector<std::vector<double> > m_inData;
  /**
   *  @brief <m_sphericalParam>, SphericalParam object with catalog parameters
  */
 LE3_2D_MASS_WL_SPHERICAL::SphericalParam m_sphericalParam;

};  // End of GetSphericalMCMaps class

}  // namespace LE3_2D_MASS_WL_SPHERICAL

#endif
