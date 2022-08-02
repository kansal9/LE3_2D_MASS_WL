/**
 * @file LE3_2D_MASS_WL_CARTESIAN/InpaintingAlgo.h
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_INPAINTINGALGO_H
#define _LE3_2D_MASS_WL_CARTESIAN_INPAINTINGALGO_H
#include "ElementsKernel/Logging.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MatrixProcess.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GenericMap.h"
#include "boost/multi_array.hpp"
#include <utility>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace LE3_2D_MASS_WL_CARTESIAN
{

/**
 * @class InpaintingAlgo
 * @brief
 *
 */
class InpaintingAlgo
{

public:

    /**
     * @brief Constructor of an inpainting object
     * @param[in] shearMap reference on the input shear map
     * @param[in] mask reference on the mask to be filled
     */

    InpaintingAlgo(const ShearMap &shearMap, GenericMap& mask,
            CartesianParam &cartesianParam);

    /**
     * @brief Destructor
     */
    virtual ~InpaintingAlgo() = default;

    /**
     * @brief Method to perform inpainting
     */
    void performInPaintingAlgo(ConvergenceMap& convergenceMap);

    /**
     * @brief Method to get the image and mask sigma
     */
    void getMaskImageSigma(Matrix& kScale, double& maskSigma, double& imSigma,
            double& maskCount, double& imCount);

    /**
     * @brief Method to enforce the powerspectrum by multiplying
     *        each wavelet coefficient by the factor imSigma/maskSigma
     */
    void multiplyWaveletCoeff(Matrix& kScale, double& maskSigma,
            double& imSigma, double& maskCount, double& imCount);

private:
    const ShearMap& m_shearMap;
    GenericMap& m_maskValues;

    int nbScales;

    float m_minThreshold;
    float m_maxThreshold;

    /**
     *  @brief <m_cartesianParam>, CartesianParam object with catalog parameters
     */
    CartesianParam m_cartesianParam;

    /**
     *  @brief <m_MP>, MatrixProcess object to perfom opertions on image matrix
     */
    MatrixProcess m_MP;

    void applyBoundariesOnWavelets(GenericMap& input, int k = 0);

    void performInversionMask(ConvergenceMap& convergenceMap, bool bModeZeros);

};
// End of InpaintingAlgo class

}// namespace LE3_2D_MASS_WL_CARTESIAN

#endif
