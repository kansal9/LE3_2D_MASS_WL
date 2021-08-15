/**
 * @file LE3_2D_MASS_WL_PEAK_COUNT/MassAperturePeakCount.h
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

#ifndef _LE3_2D_MASS_WL_PEAK_COUNT_MASSAPERTUREPEAKCOUNT_H
#define _LE3_2D_MASS_WL_PEAK_COUNT_MASSAPERTUREPEAKCOUNT_H

#include "LE3_2D_MASS_WL_PEAK_COUNT/PeakParam.h"
#include "LE3_2D_MASS_WL_PEAK_COUNT/WaveletPeakCount.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MatrixProcess.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GetMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/Projection.h"

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
 * @class MassAperturePeakCount
 * @brief
 *
 */
class MassAperturePeakCount {

public:

  /**
   * @brief Destructor
   */
  virtual ~MassAperturePeakCount() = default;

  /**
   * @brief Constructor
   * @param[in] shearMap the Shear map on which to compute the peak counting
   * @param[in] noisedShearMap the Noised Shear map needed for SNR estimation
   * @param[in] Peak Parameters for SNR estimation
   */
  MassAperturePeakCount(LE3_2D_MASS_WL_CARTESIAN::ShearMap &shearMap,
        LE3_2D_MASS_WL_CARTESIAN::ShearMap &noisedShearMap, LE3_2D_MASS_WL_PEAK_COUNT::PeakParam &PeakParam);
  /*MassAperturePeakCount(LE3_2D_MASS_WL_CARTESIAN::ShearMap &shearMap,
        LE3_2D_MASS_WL_CARTESIAN::ShearMap &noisedShearMap, LE3_2D_MASS_WL_PEAK_COUNT::PeakParam &PeakParam,
         LE3_2D_MASS_WL_CARTESIAN::CartesianParam& CParam);*/

  /**
   * @brief returns mass aperture map wrt input shear map
   * @param[in] shearMap the shear map to get mass aperture map
   * @param[in] Peak Parameters for SNR estimation
   */
  std::vector<LE3_2D_MASS_WL_CARTESIAN::Matrix> getMassApertureMap(LE3_2D_MASS_WL_CARTESIAN::ShearMap &inShearMap);

  /**
   * @brief      returns mass aperture map wrt input shear map
   * @param[in]  E mode shearMap
   * @param[in]  B mode shearMap
   * @param[in]  <iter> index of element from the list of radius 
   */
  LE3_2D_MASS_WL_CARTESIAN::Matrix createMassApertureMap(LE3_2D_MASS_WL_CARTESIAN::Matrix &shearE,
                                        LE3_2D_MASS_WL_CARTESIAN::Matrix &shearB, unsigned int iter);
  /**
   * @brief saving peaks as a FITS catalog
   * @param[in] filename the FITS catalog output filename
   */
  void saveMAPeakCatalog(const std::string& filename);

private:
  LE3_2D_MASS_WL_CARTESIAN::ShearMap m_shearMap;
  LE3_2D_MASS_WL_CARTESIAN::ShearMap m_noisedShearMap;
  /**
   *  @brief <m_MP>, MatrixProcess object to perfom opertions on image matrix
  */
 // LE3_2D_MASS_WL_CARTESIAN::MatrixProcess m_MP;
  /**
   *  @brief <m_peakParam>, PeakParam object
  */
  LE3_2D_MASS_WL_PEAK_COUNT::PeakParam m_peakParam;
  //LE3_2D_MASS_WL_CARTESIAN::CartesianParam m_cartesianParam;
  std::vector<double> radius;
  unsigned int m_sizeXaxis, m_sizeYaxis, m_sizeZaxis;
  double m_raMin, m_raMax, m_decMin, m_decMax, m_zMin, m_zMax;

};  // End of MassAperturePeakCount class

}  // namespace LE3_2D_MASS_WL_PEAK_COUNT


#endif
