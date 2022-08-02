/**
 * @file WL_2D_SPHERICAL_MASS_MAPPING/SphMapMaker.h
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
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include <vector>
#include <memory>

using namespace Euclid::WeakLensing::TwoDMass;

using LE3_2D_MASS_WL_UTILITIES::CatalogData;

namespace LE3_2D_MASS_WL_SPHERICAL
{

/**
 * @class SphMapMaker
 * @brief
 *
 */
class SphMapMaker
{

public:
    /**
     * @brief Constructor
     */
    SphMapMaker(LE3_2D_MASS_WL_SPHERICAL::SphericalParam &SphParam);

    /**
     * @brief Destructor
     */
    virtual ~SphMapMaker() = default;

    /**
     @brief  This method create Shear Map using input catalog data
     @param  <catData> input catalog data
     @retun  <Shear1, Shear2, GalCountMap> Shear healpix maps
     */
    void create_ShearMap(const CatalogData& catData,
                         Healpix_Map<double>& g1_hmap,
                         Healpix_Map<double>& g2_hmap,
                         Healpix_Map<double>& ngal_hp);
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
    SphericalParam m_SphParam;

};// End of Sph_map_maker class

}// namespace  LE3_2D_MASS_WL_SPHERICAL

#endif
