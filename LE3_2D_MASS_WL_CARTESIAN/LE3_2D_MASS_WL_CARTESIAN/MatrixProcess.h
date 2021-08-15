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
#include "LE3_2D_MASS_WL_CARTESIAN/Matrix.h"
#include<vector>
namespace LE3_2D_MASS_WL_CARTESIAN {

/**
 * @class MatrixProcess
 * @brief
 *
 */
class MatrixProcess {

public:

  /**
   * @brief Destructor
   */
  virtual ~MatrixProcess() = default;
 /**
  * @brief Constructor of an Image matrix object
  * @param[in] sizeXaxis the number of pixels in the X axis
  * @param[in] sizeYaxis the number of pixels in the Y axis
 */
  MatrixProcess(unsigned int sizeXaxis, unsigned int sizeYaxis);
  /**
   * @brief performs DCT on an input matrix
   * @param[in] input the input matrix on which to perform the DCT
   * @return the DCT image
   */
  Matrix performDCT(Matrix input);

  /**
   * @brief performs IDCT on an input DCT matrix
   * @param[in] input the input matrix on which to perform the IDCT
   * @return the matrix from the DCT matrix
   */
  Matrix performIDCT(Matrix input);

  /**
   * @brief performs (I)DCT on an input image
   * @param[in] input the input Image on which to perform the DCT
   * @param[in] blockSizeX the number of pixels to use per block on X axis
   * @param[in] blockSizeY the number of pixels to use per block on Y axis
   * @param[in] forward set to true to perform DCT, false to perform IDCT
   * @return the (I)DCT image
   */
  Matrix performDCT(Matrix input, unsigned int blockSizeX, unsigned int blockSizeY, bool forward);
  /**
   * @brief performs b spline transformation on an input Image for a given number of scales
   * @param[in] input the input Image on which to perform the transform
   * @param[in] nbScales the number of scales
   * @return a vector containing the Image of the b spline transformations for each scale
   */
  std::vector<Matrix> transformBspline(Matrix input, unsigned int nbScales);

  /**
   * @brief applies b spline transformation on a given input image and a given step for holes (algorithm a trous)
   * @param[in] input the input Image on which to perform the transform
   * @param[in] stepTrou the step defining the hole size for the algorithm
   * @return an image after transformation
   */
  Matrix smoothBspline(Matrix input, unsigned int stepTrou);

  /**
   * @brief reconstructs an image from the vector of images for each scales (returned by the transformBspline method)
   * @param[in] band the input vector of Images for each scales
   * @return an image after reconstruction
   */
  Matrix reconsBspline(std::vector<Matrix> band);
private:
  /**
   *  @brief <m_sizeXaxis>, Xaxis of matrix
  */
 unsigned int m_sizeXaxis;
  /**
   *  @brief <m_sizeYaxis>, Yaxis of matrix
  */
 unsigned int m_sizeYaxis;


};  // End of MatrixProcess class

}  // namespace LE3_2D_MASS_WL_CARTESIAN


#endif
