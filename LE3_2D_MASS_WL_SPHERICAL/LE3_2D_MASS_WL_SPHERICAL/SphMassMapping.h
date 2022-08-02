/**
 * @file LE3_2D_MASS_WL_SPHERICAL/SphMassMapping.h
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

#ifndef _LE3_2D_MASS_WL_SPHERICAL_SPH_MASS_MAPPING_H
#define _LE3_2D_MASS_WL_SPHERICAL_SPH_MASS_MAPPING_H

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"
#include <vector>
#include <memory>

using namespace LE3_2D_MASS_WL_UTILITIES;

namespace LE3_2D_MASS_WL_SPHERICAL
{

/**
 * @class SphMassMapping
 * @brief
 *
 */
class SphMassMapping
{

public:

    /**
     * @brief Constructor
     */
    SphMassMapping();

    /**
     * @brief Destructor
     */
    virtual ~SphMassMapping() = default;

    /**
     @brief  This method create Convergence Map using shear Map
     @param  <g1> Gamma1 healpix map
     @param  <g2> Gamma2 healpix map
     @retun  <kappaB, kappaE> convergence healpix maps
     */
    std::pair<Healpix_Map<double>, Healpix_Map<double> > create_SheartoConvMap(
            Healpix_Map<double>& mapE, Healpix_Map<double>& mapB);

    /**
     @brief  This method create Shear Map using convergence Map
     @param  <kE> kappaE healpix map
     @param  <kB> kappaB healpix map
     @retun  <Gamma1, Gamma2> Shear healpix maps
     */
    std::pair<Healpix_Map<double>, Healpix_Map<double> > create_ConvtoShearMap(
            Healpix_Map<double>& mapE, Healpix_Map<double>& mapB);

    /**
     @brief  This method returns nside of the map
     @param  None
     @retun  <Nside> Nside of input map for inversion
     */
    int getNside();

    /**
     @brief  This method swap the scheme of map from NEST to RING
     @param  <map> healpix map
     */
    void correctScheme(Healpix_Map<double>* map);

    /**
     @brief  Computes the reduced shear
     * This method computes the reduced shear from the original shear map and the
     * provided convergence map
     */
    void computeReducedShear_hp(
            std::pair<Healpix_Map<double>, Healpix_Map<double> >& shearMap,
            std::pair<Healpix_Map<double>, Healpix_Map<double> >& convergenceMap);

private:
    /**
     @brief  This method does mass mapping
     @param  <mapType> input healpix map type i.e. either shearMap or convMap
     @param  <inMapE> input E-Mode healpix map
     @param  <inMapB> input B-Mode healpix map
     @retun  <mapE, mapB> healpix maps
     */
    std::pair<Healpix_Map<double>, Healpix_Map<double> > sphMassMap(
            const mapType type,
            Healpix_Map<double>& mapE, Healpix_Map<double>& mapB);
    /**
     @brief  This is Nside of input map for inversion
     */
    int m_nside;
};
// End of Sph_mass_mapping class
}// namespace LE3_2D_MASS_WL_SPHERICAL
#endif
