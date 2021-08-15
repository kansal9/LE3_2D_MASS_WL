/**
 * @file LE3_2D_MASS_WL_PEAK_COUNT/SphericalPeakCount.h
 * @date 06/18/21
 * @author vanshika kansal
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

#ifndef _LE3_2D_MASS_WL_PEAK_COUNT_SPHERICALPEAKCOUNT_H
#define _LE3_2D_MASS_WL_PEAK_COUNT_SPHERICALPEAKCOUNT_H

#include "LE3_2D_MASS_WL_PEAK_COUNT/PeakParam.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_SPHERICAL/SphericalUtils.h"
#include "LE3_2D_MASS_WL_PEAK_COUNT/WaveletPeakCount.h"

#include "ElementsKernel/Logging.h"
#include <vector>
#include <memory>

#include <utility>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace LE3_2D_MASS_WL_PEAK_COUNT {

/**
 * @class SphericalPeakCount
 * @brief
 *
 */
class SphericalPeakCount {

public:

  /**
   * @brief Destructor
   */
  virtual ~SphericalPeakCount() = default;

  /**
   * @brief Constructor
   * @param[in] convergenceE the E-Mode convergence map on which to compute the peak counting
   * @param[in] Peak Parameters for SNR estimation
   */
  SphericalPeakCount(Healpix_Map<double>& convergenceE, LE3_2D_MASS_WL_PEAK_COUNT::PeakParam &PeakParam);

  /**
   * @brief saving peaks as a FITS catalog
   * @param[in] filename the FITS catalog output filename
   */
  void savePeakCatalog_hp(const std::string& filename);

  /**
   * @brief method that returns the SNR image from an input image and a noise value
   * @param[in] inputKappa the input image on which to measure the SNR
   * @param[in] globalNoise the input noise to be used to measure the SNR
   *
   * @return a SNR image
   */
  Healpix_Map<double> getSNRimage_hp(Healpix_Map<double>& inputKappa, double globalNoise);

  /**
   * @brief method that returns the peaks from a vector of bands
   * @param[in] myBand the vector containing all the scales of the wavelet decomposition of the convergence map
   *
   * @return a vector containing all the peak information: right ascension, declination, redshift, scale and SNR
   */
  std::vector<std::vector<double> > getPeaks_hp(std::vector<Healpix_Map<double> >& myBand, double globalNoise);

private:
  /**
   * @brief E-Mode Convergence Map
   */
  Healpix_Map<double> KappaE;

  /**
   *  @brief <m_peakParam>, PeakParam object
  */
  LE3_2D_MASS_WL_PEAK_COUNT::PeakParam m_peakParam;

  /**
  @brief  This is Nside of input map
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
   * @brief LMax
   */
  int m_nlmax;

  /**
   * @brief number of scales
   */
  int m_nbScales;

};  // End of SphericalPeakCount class

}  // namespace LE3_2D_MASS_WL_PEAK_COUNT

#endif
