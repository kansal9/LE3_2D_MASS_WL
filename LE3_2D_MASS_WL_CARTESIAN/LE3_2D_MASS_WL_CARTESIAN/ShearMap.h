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

#include "GenericMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"

namespace LE3_2D_MASS_WL_CARTESIAN
{
class ConvergenceMap;
/**
 * @class ShearMap
 * @brief
 *
 */
class ShearMap: public GenericMap
{

public:

    using GenericMap::GenericMap;
    using GenericMap::operator=;

    /**
     * @brief Copy constructor of a shear map
     * @param[in] copyMap the map to copy
     * @param[in] copyValues if true does copy data from copyMap, if false copy
     *            only metadata
     * This copy constructor builds a copy of the input shear map
     */
    ShearMap(GenericMap const& copyMap, bool copyValues = true);

    /**
     * @brief Returns a ConvergenceMap using K&S algorithm
     * @param[out] a convergence Map corresponding to the input Shear Map
     * This method creates and returns a Convergence Map from the given ShearMap
     * using Kaiser & Squires algorithm
     */
    void getConvMap(ConvergenceMap& outputConvMap) const;

    /**
     * @brief Corrects the reduced shear
     * This method corrects the reduced shear using the original shear map and the
     * provided convergence map
     **/
    void correctReducedShear(const ConvergenceMap& inputConvMap,
            bool useKappaB = false);

private:

};
// End of ShearMap class

}// namespace LE3_2D_MASS_WL_CARTESIAN

#endif
