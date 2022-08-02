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

#include "LE3_2D_MASS_WL_CARTESIAN/GenericMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

namespace LE3_2D_MASS_WL_CARTESIAN
{
class ShearMap;
/**
 * @class ConvergenceMap
 * @brief
 *
 */
class ConvergenceMap: public GenericMap
{

public:

    using GenericMap::GenericMap;
    using GenericMap::operator=;

    /**
     * @brief Copy constructor of a convergence map
     * @param[in] copyMap the map to copy
     * @param[in] copyValues if true does copy data from copyMap, if false copy
     *            only metadata
     * This copy constructor builds a copy of the input convergence map
     */
    ConvergenceMap(GenericMap const& copyMap, bool copyValues = true);

    /**
     * @brief Returns a ShearMap using K&S algorithm
     * @param[out] a ShearMap corresponding to the input ConvergenceMap
     * This method creates and returns a ShearMap from the given ConvergenceMap
     * using Kaiser & Squires algorithm
     */
    void getShearMap(ShearMap& outputShearMap) const;

    /**
     * @brief Applies a gaussian filter on the convergence map
     * @param[in] sigma sigma of the gaussian kernel
     * @param[in] k image on which to apply the filter
     * This method applies a gaussian filter to the Map of width sigma
     * The map is supposed to be a cube with axis 0 standing for real part and
     * axis 1 standing for imaginary part
     */
    void applyGaussianFilter(float sigma, int k = 0);

    /**
     * @brief Applies a gaussian filter on the convergence map
     * @param[in] sigmaX sigma of the gaussian kernel on the X direction
     * @param[in] sigmaY sigma of the gaussian kernel on the Y direction
     * @param[in] k image on which to apply the filter
     * This method applies a gaussian filter to the Map of width
     * sigmaX and sigmaY in X and Y directions respectively.
     * The map is supposed to be a cube with axis 0 standing for real part and
     * axis 1 standing for imaginary part
     */
    void applyGaussianFilter(float sigmax, float sigmay, int k = 0);

    /**
     * @brief Sets the imaginary part kappaB to zero
     * This method sets the imaginary part of the complex convergence (kappaB) to zero
     **/
    void setKappaBToZero();

private:

};
// End of ConvergenceMap class

}// namespace LE3_2D_MASS_WL_CARTESIAN

#endif
