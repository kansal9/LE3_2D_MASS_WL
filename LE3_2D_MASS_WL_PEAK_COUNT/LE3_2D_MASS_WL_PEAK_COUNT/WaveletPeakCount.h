/**
 * @file LE3_2D_MASS_WL_PEAK_COUNT/WaveletPeakCount.h
 * @date 10/07/20
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

#ifndef _LE3_2D_MASS_WL_PEAK_COUNT_WAVELETPEAKCOUNT_H
#define _LE3_2D_MASS_WL_PEAK_COUNT_WAVELETPEAKCOUNT_H

#include "LE3_2D_MASS_WL_CARTESIAN/MatrixProcess.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GetMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/Projection.h"
#include "LE3_2D_MASS_WL_PEAK_COUNT/PeakParam.h"

#include "boost/multi_array.hpp"
#include <utility>
#include <algorithm>
#include <iostream>
#include <iomanip>

using LE3_2D_MASS_WL_CARTESIAN::GetMap;
using LE3_2D_MASS_WL_CARTESIAN::MatrixProcess;
using LE3_2D_MASS_WL_CARTESIAN::Projection;

namespace LE3_2D_MASS_WL_PEAK_COUNT {

/**
 * @class WaveletPeakCount
 * @brief this class will use to get the peak catalog using wavelets
 *
 */
class WaveletPeakCount {

public:

  /**
   * @brief Destructor
   */
  virtual ~WaveletPeakCount() = default;

  /**
   * @brief Constructor
   * @param[in] convMap the convergence map on which to compute the peak counting
   * @param[in] Peak Parameters for SNR estimation
   */
  WaveletPeakCount(LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap &convMap, LE3_2D_MASS_WL_PEAK_COUNT::PeakParam &PeakParam);

  /**
   * @brief saving peaks as a FITS catalog
   * @param[in] filename the FITS catalog output filename
   */
  void savePeakCatalog(const std::string& filename);

  /**
   * @brief method that returns the peaks from a vector of bands
   * @param[in] myBand the vector containing all the scales of the wavelet decomposition of the convergence map
   *
   * @return a vector containing all the peak information: right ascension, declination, redshift, scale and SNR
   */
 std::vector<std::vector<double> > getPeaks(std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> &myBand, double globalNoise);

  /**
   * @brief method that returns the SNR image from an input image and a noise value
   * @param[in] inputKappaImage the input image on which to measure the SNR
   * @param[in] globalNoise the input noise to be used to measure the SNR
   *
   * @return a SNR image
   */
  LE3_2D_MASS_WL_CARTESIAN::Matrix getSNRimage(LE3_2D_MASS_WL_CARTESIAN::Matrix inputKappaImage, double globalNoise);

private:

  LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap m_convMap;
  /**
   *  @brief <m_MP>, MatrixProcess object to perfom opertions on image matrix
  */
  LE3_2D_MASS_WL_CARTESIAN::MatrixProcess m_MP;
  /**
   *  @brief <m_peakParam>, PeakParam object
  */
  LE3_2D_MASS_WL_PEAK_COUNT::PeakParam m_peakParam;

  unsigned int m_sizeXaxis;
  unsigned int m_sizeYaxis;
  unsigned int m_sizeZaxis;

  double raMin, raMax, decMin, decMax, m_zMin, m_zMax;

};  // End of WaveletPeakCount class

  /**
   * @brief      writing FITS file
   * @param[in]  <filename> the output filename
   * @param[in]  <Data> that has to be written in FITS file
   */
  void writePeakCatalog (const std::string& filename, std::vector<std::vector<double> > Data);

}  // namespace LE3_2D_MASS_WL_PEAK_COUNT


#endif
