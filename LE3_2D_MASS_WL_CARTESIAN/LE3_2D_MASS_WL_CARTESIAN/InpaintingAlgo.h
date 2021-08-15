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
#include "LE3_2D_MASS_WL_CARTESIAN/GetMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "boost/multi_array.hpp"
#include <utility>
#include <algorithm>
#include <iostream>
#include <iomanip>

using LE3_2D_MASS_WL_CARTESIAN::MatrixProcess;
using LE3_2D_MASS_WL_CARTESIAN::GetMap;

namespace LE3_2D_MASS_WL_CARTESIAN {

/**
 * @class InpaintingAlgo
 * @brief
 *
 */
class InpaintingAlgo {

public:

  /**
   * @brief Destructor
   */
 virtual ~InpaintingAlgo();

  /**
   * @brief Constructor of an inpainting object
   * @param[in] shearMap the input shear map
   * @param[in] convMap the input convergence map
   */

 InpaintingAlgo(ShearMap &shearMap, ConvergenceMap &convMap,
                 LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam);

  /**
   * @brief Copy constructor of an InPainting object
   * @param[in] copy the InPainting object to copy
   */
 InpaintingAlgo(const InpaintingAlgo& copy);

  /**
   * @brief Method to perform inpainting
   */
 ConvergenceMap* performInPaintingAlgo();

  /**
   * @brief Method to perform inpainting for given block size of the image
   * @param[in] blockSizeX size of the block on X dimensions
   * @param[in] blockSizeY size of the block on Y dimensions
   */
 ConvergenceMap* performInPaintingAlgo(unsigned int blockSizeX, unsigned int blockSizeY);

  /**
   * @brief Method to get the image and mask sigma
   */
 void getMaskImageSigma (LE3_2D_MASS_WL_CARTESIAN::Matrix& Image, double& maskSigma, double& imSigma,
                         double& maskCount, double& imCount);

  /**
   * @brief Method to enforce the powerspectrum by multiplying 
   *        each wavelet coefficient by the factor imSigma/maskSigma
   */
 void multiplyWaveletCoeff (LE3_2D_MASS_WL_CARTESIAN::Matrix& Image, double& maskSigma, double& imSigma,
                         double& maskCount, double& imCount);

private:
 ShearMap shearMap;
 ConvergenceMap convMap;

 unsigned int Xaxis, Yaxis, Zaxis;
 int nbScales;
 boost::multi_array<int, 3> *m_maskValues;

 float m_minThreshold;
 float m_maxThreshold;
  /**
   *  @brief <m_cartesianParam>, CartesianParam object with catalog parameters
  */
 LE3_2D_MASS_WL_CARTESIAN::CartesianParam m_cartesianParam;
  /**
   *  @brief <m_MP>, MatrixProcess object to perfom opertions on image matrix
  */
 MatrixProcess m_MP;
 Matrix applyBoundariesOnWavelets(Matrix input);

 ConvergenceMap performInversionMask(Matrix kappaE, Matrix kappaB, bool bModeZeros);

};  // End of InpaintingAlgo class

}  // namespace LE3_2D_MASS_WL_CARTESIAN

#endif
