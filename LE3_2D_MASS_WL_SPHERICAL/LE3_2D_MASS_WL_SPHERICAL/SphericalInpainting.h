/**
 * @file LE3_2D_MASS_WL_SPHERICAL/SphericalInpainting.h
 * @date 02/01/21
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

#ifndef _LE3_2D_MASS_WL_SPHERICAL_SPHERICALINPAINTING_H
#define _LE3_2D_MASS_WL_SPHERICAL_SPHERICALINPAINTING_H

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalIO.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalUtils.h"

#include "ElementsKernel/Logging.h"
#include <vector>
#include <memory>

#include <cmath>
#include <complex>
#include <fftw3.h>
#include <iostream>
#include <fstream>
#include "SphMassMapping.h"

namespace LE3_2D_MASS_WL_SPHERICAL
{

/**
 * @class SphericalInpainting
 * @brief
 *
 */
class SphericalInpainting
{

public:

    /**
     * @brief Constructor of an inpainting object
     * @param[in] shear the input shear map
     * @param[in] convergence the input convergence map
     * @param[in] SphericalParam the SphericalParam object with input parameters required for inpainting
     */

    SphericalInpainting(
            std::pair<Healpix_Map<double>, Healpix_Map<double> >& shear,
            std::pair<Healpix_Map<double>, Healpix_Map<double> >& convergence,
            SphericalParam &SphericalParam);

    /**
     * @brief Destructor
     */
    virtual ~SphericalInpainting() = default;

    /**
     * @brief Method to perform inpainting on Sphere
     */
    std::pair<Healpix_Map<double>, Healpix_Map<double> > performInpainting();

    /**
     * @brief perform mask inversion
     */
    std::pair<Healpix_Map<double>, Healpix_Map<double> > maskInversion(
            Healpix_Map<double>& kE, Healpix_Map<double>& kB);

    /**
     * @brief apply wavelet constraints on healpix map
     */
    void applyConstraintsOnWavelets(Healpix_Map<double>& map);

    /**
     * @brief perform wavelets on healpix map
     */
    Healpix_Map<double> performWavelet(Healpix_Map<double>& map);

private:
    /**
     * @brief Minimum Threshold
     */
    double m_minThreshold;

    /**
     * @brief Maximum Threshold
     */
    double m_maxThreshold;

    /**
     * @brief LMax
     */
    int m_nlmax;

    /**
     * @brief Resolution of Map (Nside)
     */
    int m_nside;

    /**
     * @brief number of pixels (NPIX) of the map
     */
    int m_npix;

    /**
     * @brief resolution order of the map
     */
    int m_order;

    /**
     * @brief number of scales
     */
    int m_nbScales;

    /**
     *  @brief <m_SphParam>, SphericalParam object with input parameters
     */
    SphericalParam m_SphParam;

    /**
     * @brief E-Mode Shear Map
     */
    Healpix_Map<double> ShearE;

    /**
     * @brief B-Mode Shear Map
     */
    Healpix_Map<double> ShearB;

    /**
     * @brief E-Mode Convergence Map
     */
    Healpix_Map<double> KappaE;

    /**
     * @brief B-Mode Convergence Map
     */
    Healpix_Map<double> KappaB;

    /**
     * @brief Mask Values inside a map
     */
    Healpix_Map<int> Mask;

}; // End of SphericalInpainting class

}// namespace LE3_2D_MASS_WL_SPHERICAL

#endif
