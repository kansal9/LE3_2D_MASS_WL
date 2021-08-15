/**
 * @file LE3_2D_MASS_WL_SPHERICAL/SphericalUtils.h
 * @date 06/17/21
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

#ifndef _SPHERICALUTILS_H
#define _SPHERICALUTILS_H

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "ElementsKernel/Logging.h"

#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"

namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {
   namespace Spherical {

  /**
  @brief  This method read input Spherical Parameter file
  @param  <ParamFile> name of the input parameter file in xml format
  @retun  <params> SphericalParam object to get input parameters
  */
  void readSphericalParameterFile (const boost::filesystem::path& ParamFile,
                                   LE3_2D_MASS_WL_SPHERICAL::SphericalParam &params);
  /**
  * @brief This method applies a gaussian filter
  * @param <sigma> sigma of the gaussian kernel
  * @param <map> map in healpix format
  * This method applies a gaussian filter to the Map of width sigma
  */
  void applyGaussianFilter_hp(Healpix_Map<double>& map, double sigma);

  /**
   * @brief perform B-spline scaling function of order 3
   */
  double b3_spline (double val);

  /**
   * @brief sets all values below a threshold to zero
   */
  void applyThreshold (Alm<xcomplex<double> >& alphaT, double threshVal, int m_nlmax);

  /**
   * @brief returns the maximum absolute value of alm
   */
  double max_abs_alm(Alm<xcomplex<double> >& alm, int m_nlmax);

  /**
   * @brief returns the standard deviation of input map
   */
  double getStd(Healpix_Map<double>& map);

  /**
   * @brief get the filter values
   */
  void getFilter(double* array, double lc, int m_nlmax);

  /**
   * @brief performs b spline transformation on an input healpix map for a given number of scales
   * @param[in] the input healpix map on which to perform the transform
   * @param[in] the number of scales
   * @return a vector containing the healpix map of the b spline transformations for each scale
   */
  std::vector<Healpix_Map<double> > transformBspline_hp(Healpix_Map<double>& map, int m_nbScales);

  /**
   * @brief reconstructs a healpix map from the vector of healpix maps for each scales
            (returned by the transformBspline method)
   * @param[in] band the input vector of healpix map for each scales
   * @return a healpix map after reconstruction
   */
  Healpix_Map<double> reconsBspline_hp(std::vector<Healpix_Map<double> >& band);

  /**
   * @brief applies b spline transformation on a given input healpix map and a given step for holes (algorithm a trous)
   * @param[in] the input healpix map on which to perform the transform
   * @param[in] scale the step defining the hole size for the algorithm
   * @return a healpix map after transformation
   */
  Healpix_Map<double> smoothBspline_hp(Healpix_Map<double>& map, int scale);

   } /* namespace Spherical */
  } /* namespace TwoDMass */
 } /* namespace WeakLensing */
} /* namespace Euclid */
#endif
