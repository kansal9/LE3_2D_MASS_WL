/**
 * @file LE3_2D_MASS_WL_PEAK_COUNT/CartesianPeakCount.h
 * @date 06/30/22
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

#ifndef _LE3_2D_MASS_WL_PEAK_COUNT_PEAKCOUNTCARTESIAN_H
#define _LE3_2D_MASS_WL_PEAK_COUNT_PEAKCOUNTCARTESIAN_H

#include "LE3_2D_MASS_WL_CARTESIAN/MatrixProcess.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/Projection.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GenericMap.h"
#include "LE3_2D_MASS_WL_PEAK_COUNT/GenericPeakCount.h"
#include "LE3_2D_MASS_WL_PEAK_COUNT/PeakParam.h"
#include "LE3_2D_MASS_WL_UTILITIES/PatchDef.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

using LE3_2D_MASS_WL_CARTESIAN::ShearMap;
using LE3_2D_MASS_WL_CARTESIAN::MatrixProcess;
using LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap;
using LE3_2D_MASS_WL_CARTESIAN::Projection;
using LE3_2D_MASS_WL_UTILITIES::Matrix;
using LE3_2D_MASS_WL_UTILITIES::PatchDef;

namespace LE3_2D_MASS_WL_PEAK_COUNT
{

/**
 * @class PeakCountCartesian
 * @brief
 *
 */

enum estimator {
    APERTURE_MASS = 0,
    WAVELET = 1
};

class CartesianPeakCount : public GenericPeakCount
{

public:

    CartesianPeakCount(PatchDef &patch, double zmin, double zmax,
            unsigned int estimator);

    /**
     * @brief Destructor
     */
    virtual ~CartesianPeakCount() = default;

    void findPeaksAtTheta(Matrix &map, double theta);

    void findPeaksFromShear(ShearMap &shearMap, ShearMap &noisedShearMap,
            std::vector<double> radius);

    void findPeaksFromConvergence(ConvergenceMap &convMap, unsigned int nbScales);

    void getMassApertureMap(const ShearMap &inShearMap, ShearMap &outShearMap,
            std::vector<double> radius) const;

    void createMassApertureMap(const Matrix &shearE, const Matrix &shearB,
            Matrix& outShear, double radius) const;

private:

    Projection m_projection;

    unsigned int m_estimator;

    PatchDef m_patch;

    double m_zmin, m_zmax;

};
// End of PeakCountCartesian class

}// namespace LE3_2D_MASS_WL_PEAK_COUNT

#endif // _LE3_2D_MASS_WL_PEAK_COUNT_PEAKCOUNTCARTESIAN_H
