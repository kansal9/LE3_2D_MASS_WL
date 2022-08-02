/**
 * @file src/lib/CartesianPeakCount.cpp
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

#include "LE3_2D_MASS_WL_PEAK_COUNT/CartesianPeakCount.h"

using namespace LE3_2D_MASS_WL_UTILITIES;

namespace LE3_2D_MASS_WL_PEAK_COUNT
{

CartesianPeakCount::CartesianPeakCount(PatchDef &patch, double zmin,
        double zmax, unsigned int estimator) :
        m_projection(), m_estimator(estimator), m_patch(patch),
        m_zmin(zmin), m_zmax(zmax)
{

}

void CartesianPeakCount::findPeaksAtTheta(Matrix &map, double theta)
{
    std::pair<double, double> radec;
    double ra0 = 0.5 * (m_patch.getRaMin() + m_patch.getRaMax()); // in degree
    double dec0 = 0.5 * (m_patch.getDecMin() + m_patch.getDecMax()); // in degree
    double raRange = (m_patch.getRaMax() - m_patch.getRaMin()) * M_PI / 180.; // in rad
    double decRange = (m_patch.getDecMax() - m_patch.getDecMin()) * M_PI / 180.; // in rad

    for (unsigned int i = 0; i < map.getXdim(); i++)
    {
        for (unsigned int j = 0; j < map.getYdim(); j++)
        {
            if (map.isLocalMax(i, j))
            {
                // compute ra/dec
                double tmpx = (i + 0.5) * raRange / map.getXdim() - 0.5 * raRange;
                double tmpy = (j + 0.5) * decRange / map.getYdim() - 0.5 * decRange;
                radec = m_projection.getInverseGnomonicProjection(tmpx, tmpy, ra0, dec0);
                // fill cat
                m_peakCat.addPeak(radec.first, radec.second, m_zmin, m_zmax,
                        theta, map(i,j));
            }
        }
    }
}

void CartesianPeakCount::findPeaksFromShear(ShearMap &shearMap,
        ShearMap &noisedShearMap, std::vector<double> radius)
{
    // construct shear map with number of radius as 3rd dimension
    ShearMap shearBands(shearMap.getXdim(), shearMap.getYdim(), radius.size());
    ShearMap noisedShearBands(noisedShearMap.getXdim(),
            noisedShearMap.getYdim(), radius.size());

    // fill shearBands and noisedShearBands with mass aperture maps wrt each radius
    getMassApertureMap(shearMap, shearBands, radius);
    getMassApertureMap(noisedShearMap, noisedShearBands, radius);

    // loop over radius
    for (size_t r = 0; r < radius.size(); ++r)
    {
        double stdevNoise = noisedShearMap.getSigma(r);
        shearBands[r] *= 1. / stdevNoise;
        double theta = r;
        findPeaksAtTheta(shearBands[r], theta);
    }
}

void CartesianPeakCount::findPeaksFromConvergence(ConvergenceMap &convMap,
        unsigned int nbScales)
{
    // Get the noise on the input kappa
    double stdevNoise = convMap.getSigma();

    // Perform the kappa map decomposition
    std::vector<Matrix> kappaBands;
    MatrixProcess mp(convMap.getXdim(), convMap.getYdim());
    mp.transformBspline(convMap[0], kappaBands, nbScales);

    // loop over scales
    for (size_t scale = 0; scale < kappaBands.size() - 1; ++scale)
    {
        kappaBands[scale] *= (norm[scale] / stdevNoise);
        double theta = pow(2, scale + 1) * convMap.getPixelSize() * 60.;
        findPeaksAtTheta(kappaBands[scale], theta);
    }
}


void CartesianPeakCount::getMassApertureMap(const ShearMap &inShearMap,
        ShearMap &outShearMap, std::vector<double> radius) const
{
    for (size_t iter = 0; iter < radius.size(); iter++)
    {
        createMassApertureMap(inShearMap[0], inShearMap[1], outShearMap[iter],
                iter);
    }
}

void CartesianPeakCount::createMassApertureMap(const Matrix &shearE,
        const Matrix &shearB, Matrix& outShear, double radius) const
{
    double rmax = radius * 4.;

    for (int i0 = rmax; i0 < shearE.getXdim() - rmax; i0++)
    {
        for (int j0 = rmax; j0 < shearE.getYdim() - rmax; j0++)
        {
            for (int i = 0; i < 2. * rmax; i++)
            {
                for (int j = 0; j < 2. * rmax; j++)
                {
                    int ii, jj, a, b;
                    double a1, a2, phi, kernel, r;
                    ii = i0 + (i - rmax);
                    jj = j0 + (j - rmax);
                    a = (ii - i0);
                    b = (jj - j0);
                    phi = atan2(a, b);
                    r = sqrt(pow(a, 2.) + pow(b, 2.));
                    kernel = 0;
                    if (r < rmax)
                    {
                        kernel = pow(r, 2) / (4. * M_PI * pow(radius, 4))
                                * exp(-pow(r, 2.)/(2. * pow(radius, 2.)));
                    }
                    a1 = shearE(ii, jj) * cos(2. * phi);
                    a2 = shearB(ii, jj) * sin(2. * phi);
                    outShear(i0, j0) += (a1 - a2) * kernel;
                }
            }
        }
    }
}

}  // namespace LE3_2D_MASS_WL_PEAK_COUNT

