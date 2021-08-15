/**
 * @file LE3_2D_MASS_WL_CARTESIAN/ReconstructMR.h
 * @date 03/12/21
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_RECONSTRUCTMR_H
#define _LE3_2D_MASS_WL_CARTESIAN_RECONSTRUCTMR_H

#include <algorithm>
#include <cmath>
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "ElementsKernel/Logging.h"

namespace LE3_2D_MASS_WL_CARTESIAN {

/**
 * @class ReconstructMR
 * @brief
 *
 */
class ReconstructMR {

public:
  /**
   * @brief Constructor
   */
 ReconstructMR();
  /**
   * @brief Destructor
   */
  ~ReconstructMR();

  /**
   * @brief copy Constructor
   */
  ReconstructMR(const ReconstructMR& copy);

  /**
   * @brief       calculate the derivative of hn in case of Gaussian noise (sigma=1)
   *              hn = x/Sigma^2*erfc(x/(sqrt(2)sigma)) + sqrt(2/PI)/sigma[1-exp(-x^2/(2sigma^2))]
   * @param[in]   the input value
   * @return      derivative of hn
   */
  double grad_hn_sig1(double Val);

  /**
   * @brief       calculate the derivative of hs in case of Gaussian noise (sigma=1)
   *              hs = -x/Sigma^2*erf(x/(sqrt(2)sigma)) + sqrt(2/PI)/sigma[1-exp(-x^2/(2sigma^2))]
   * @param[in]   the input value
   * @return      derivative of hs
   */
  double grad_hs_sig1(double Val);

  /**
   * @brief       calculate the derivative of hs at input Val
   * @param[in]   the input value
   * @return      derivative of hs
   */
   double grad_hs(double Val);

  /**
   * @brief       calculate the derivative of hn at input Val
   * @param[in]   the input value
   * @return      derivative of hn
   */
   double grad_hn(double Val);

  /**
   * @brief       calculate by the dichotomy method the solution wich minimize:
   *              J = hs(y-x) + Alpha hn(x), It is obtained when:  grad_hs(y-x) = Alpha grad_hn(x)
   */
   double filter (double CoefDat, double Alpha, double Sigma=1.);

private:
  /**
   * @brief  internal constant = sqrt(2/PI)
   */
   double C1;

  /**
   * @brief  internal constant = sqrt(2)
   */
   double C2;

  /**
   * @brief  array for the derivative of Hs
   */
   double *TabHsGauss;

  /**
   * @brief  array for the derivative of Hn
   */
   double *TabHnGauss;

  /**
   * @brief  size of the  array
   */
   int Np;

  /**
   * @brief  Step used for the probability integration
   */
   float Step;

};  // End of ReconstructMR class

}  // namespace LE3_2D_MASS_WL_CARTESIAN

#endif
