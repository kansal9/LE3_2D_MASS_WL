/**
 * @file LE3_2D_MASS_WL_CARTESIAN/MatrixProcess.h
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_MATRIXPROCESS_H
#define _LE3_2D_MASS_WL_CARTESIAN_MATRIXPROCESS_H
#include "LE3_2D_MASS_WL_CARTESIAN/GenericMap.h"
#include <vector>
#include <fftw3.h>

using LE3_2D_MASS_WL_UTILITIES::Matrix;

namespace LE3_2D_MASS_WL_CARTESIAN
{

/**
 * @class MatrixProcess
 * @brief
 *
 */
class MatrixProcess
{

public:

    /**
     * @brief Constructor of an Image matrix object
     * @param[in] sizeXaxis the number of pixels in the X axis
     * @param[in] sizeYaxis the number of pixels in the Y axis
     */
    MatrixProcess(unsigned int sizeXaxis, unsigned int sizeYaxis);

    /**
     * @brief performs DCT on an input matrix
     * @param[in] data array on which to perform the DCT
     * @param[out] data array that contains DCT image
     */
    void performDCT(double* input, double* output);

    /**
     * @brief performs IDCT on an input DCT matrix
     * @param[in] input the input matrix on which to perform the IDCT
     * @param[out] the matrix from the DCT matrix
     */
    void performIDCT(double* input, double* output);

    /**
     * @brief performs DFT on an input matrix
     * @param[in] data array on which to perform the DFT
     * @param[out] data array that contains DFT image
     */
    void performDFT(double* input, fftw_complex* output, double dftFactor = 1);

    /**
     * @brief performs IDFT on an input IDFT matrix
     * @param[in] input the input matrix on which to perform the IDFT
     * @param[out] the matrix from the IDFT matrix
     */
    void performIDFT(fftw_complex* input, double* output, double dftFactor = 1);

    /**
     * @brief performs b spline transformation on an input Image for a given number of scales
     * @param[in] input the input map on which to perform the transform
     * @param[in] nbScales the number of scales
     */
    void transformBspline(Matrix& input, std::vector<Matrix>& scalesValues,
            unsigned int nbScales);

    /**
     * @brief applies b spline transformation for a given step for holes (algorithm a trous)
     * @param[in] input the input matrix
     * @param[out] output the output matrix
     * @param[in] stepTrou the step defining the hole size for the algorithm
     */
    void smoothBspline(Matrix& input, Matrix& output, unsigned int stepTrou);

    /**
     * @brief reconstructs an image from the array of matrices for each scales
     * @param[out] map that will contain the result
     */
    void reconsBspline(std::vector<Matrix>& scalesValues, Matrix& imageOut);

private:
    /**
     *  @brief <m_sizeXaxis>, Xaxis of matrix
     */
    int m_sizeXaxis;
    /**
     *  @brief <m_sizeYaxis>, Yaxis of matrix
     */
    int m_sizeYaxis;
};
// End of MatrixProcess class

}// namespace LE3_2D_MASS_WL_CARTESIAN

#endif
