/**
 * @file LE3_2D_MASS_WL_CARTESIAN/FilterMR.h
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_FILTERMR_H
#define _LE3_2D_MASS_WL_CARTESIAN_FILTERMR_H

#include "LE3_2D_MASS_WL_CARTESIAN/MatrixProcess.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ReconstructMR.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GenericMap.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

#include "ElementsKernel/Logging.h"
#include <utility>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>

namespace LE3_2D_MASS_WL_CARTESIAN
{

/**
 * @class FilterMR
 * @brief
 *
 */
class FilterMR
{

public:

    /**
     * @brief Destructor
     */
    virtual ~FilterMR() = default;

    /**
     * @brief Constructor of class FilterMR
     * @param[in] convMap the input convergence map
     * @param[in] PositiveCons Apply the positivty constraint
     * @param[in] KillLastScale Remove the smoothed plane
     * @param[in] KillIsol suppress isolated pixels
     * @param[in] niter number of loops in reconstruction iterative process (default is 10)
     * @param[in] FirstScale used for detection (default is 0)
     * @param[in] CartesianParam object with input parameters from input parameter file
     */
    FilterMR(LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap &convMap,
            bool PositiveCons, bool KillLastScale, bool KillIsol, int niter,
            int FirstScale,
            LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam);

    //FilterMR(LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap &convMap,
    //           LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam);

    /**
     * @brief Copy constructor of FilterMR object
     * @param[in] copy the FilterMR object to copy
     */
    FilterMR(const FilterMR& copy);

    double fdr_pv(std::vector<double> &PValue, double alpha);

    /**
     * @brief this function calculates the inverse error function
     */
    double myErfInv(double x);
    double xerfc(double F);

    /**
     * @brief return the integrated probility (repartition function) that a noise
     *        produce the wavelet coefficient Val at scale s and at position i,j
     */
    double getProbility(double val, double sigLevel);

    /**
     * @brief Compute FDR threshold
     */
    void get_fdr(std::vector<LE3_2D_MASS_WL_UTILITIES::Matrix>& myBand,
            std::vector<double>& NSigma);

    void mr_support(std::vector<LE3_2D_MASS_WL_UTILITIES::Matrix>& myBand,
            std::vector<double>& NSigma,
            std::vector<LE3_2D_MASS_WL_UTILITIES::Matrix>& mr_myBand);

    /**
     * @brief Method to perform Filtering
     */
    void performFiltering(ConvergenceMap& outConvergenceMap);

    /**
     * @brief Apply the multiscale entropy filter
     */
    void applyfilter(std::vector<LE3_2D_MASS_WL_UTILITIES::Matrix>& band,
            std::vector<double>& NSigma, ConvergenceMap& outConvergenceMap);

    /**
     * @brief Apply the multiscale entropy to correct the wavelet coefficients
     */
    void fixedAlphaFilter(std::vector<LE3_2D_MASS_WL_UTILITIES::Matrix>& band,
            std::vector<double>& TabAlpha, std::vector<double>& NSigma,
            LE3_2D_MASS_WL_CARTESIAN::ReconstructMR & restore);

    /**
     * @brief perform wavelet reconstruction
     */
    void reconstruct_optimized(
            std::vector<LE3_2D_MASS_WL_UTILITIES::Matrix>& band,
            std::vector<LE3_2D_MASS_WL_UTILITIES::Matrix>& mr_band,
            ConvergenceMap& outConvergenceMap);

    /**
     * @brief perform reversible reconstruction
     */
    void reconstruct(std::vector<LE3_2D_MASS_WL_UTILITIES::Matrix>& band,
            LE3_2D_MASS_WL_UTILITIES::Matrix& outConvergenceMap);

    /**
     * @brief removes the isolated pixels from the image
     */
    void kill_isolated(LE3_2D_MASS_WL_UTILITIES::Matrix& image);

    /**
     * @brief this method estimate the new alpha parameters
     */
    void estimateNewAlphaParam(
            std::vector<LE3_2D_MASS_WL_UTILITIES::Matrix>& band,
            std::vector<LE3_2D_MASS_WL_UTILITIES::Matrix>& band_sol,
            std::vector<double>& regulMin, std::vector<double>& regulMax,
            std::vector<double>& TabAlpha, std::vector<double>& TabDelta);

    /**
     * @brief this method provide support to mr_support method
     */
    void set_support(std::vector<LE3_2D_MASS_WL_UTILITIES::Matrix>& myBand,
            std::vector<double>& NSigma, int it);

private:
    /**
     * @brief <convMap> Convergence Map
     */
    LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap m_convMap;

    /**
     * @brief <m_Xaxis, m_Yaxis, m_Zaxis> axis of Map
     */
    int m_Xaxis, m_Yaxis, m_Zaxis;

    /**
     * @brief <nbScales> number of scales
     */
    int m_nbScales;

    /**
     * @brief <m_niter> number of loops in reconstruction iterative process
     */
    int m_niter;

    /**
     * @brief <m_FirstScale> FirstScale used for detection
     */
    int m_FirstScale;

    /**
     * @brief <MaxIter> Maximum number of Multiscale Entropy Iteration
     */
    int m_MaxIter;

    /**
     * @brief <m_SigmaNoise> Standard Deviation of the Map
     */
    double m_SigmaNoise;

    /**
     * @brief <m_FDR> FDR threshold
     */
    double m_FDR;
    double m_avar;

    /**
     * @brief <m_PositveConstraint> Apply the positivty constraint
     */
    bool m_PositveConstraint;

    /**
     * @brief <m_KillLastScale> Remove the smoothed plane
     */
    bool m_KillLastScale;

    /**
     * @brief <m_KillIsol> suppress isolated pixels
     */
    bool m_KillIsol;

    /**
     *  @brief <m_cartesianParam>, CartesianParam object with catalog parameters
     */
    LE3_2D_MASS_WL_CARTESIAN::CartesianParam m_cartesianParam;
    /**
     *  @brief <m_MP>, MatrixProcess object to perfom opertions on image matrix
     */
    LE3_2D_MASS_WL_CARTESIAN::MatrixProcess m_MP;
};
// End of FilterMR class

}// namespace LE3_2D_MASS_WL_CARTESIAN

#endif
