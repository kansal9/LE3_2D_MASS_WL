/**
 * @file src/lib/InpaintingAlgo.cpp
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

#include "LE3_2D_MASS_WL_CARTESIAN/InpaintingAlgo.h"

#include <map>
#include <fftw3.h>
#include <cmath>

static Elements::Logging logger = Elements::Logging::getLogger("Inpainting");

namespace LE3_2D_MASS_WL_CARTESIAN
{

InpaintingAlgo::InpaintingAlgo(const ShearMap& shearMap, GenericMap& mask,
        CartesianParam& cartesianParam) :
        m_shearMap(shearMap), m_maskValues(mask), m_cartesianParam(
                cartesianParam), m_MP(shearMap.getXdim(), shearMap.getYdim())
{
    m_minThreshold = 0.0;
    m_maxThreshold = 0.5;
    nbScales = m_cartesianParam.getNInpScale();

    if (nbScales == 0)
    {
        nbScales = int(log(shearMap.getXdim()) / log(2.)) - 2.;
    }

    logger.info() << "axis dim: " << m_shearMap.getXdim() << " "
            << m_shearMap.getYdim() << " " << m_shearMap.getZdim();
    logger.info() << "number of scales: " << nbScales;

    unsigned int count(0);
    unsigned int count2(0);

    // Initialize the mask
    for (size_t j = 0; j < m_shearMap.getYdim(); ++j)
    {
        for (size_t i = 0; i < m_shearMap.getXdim(); ++i)
        {
            if (fabs(shearMap.getBinValue(i, j, 0)) < 0.0000000001
                    && fabs(shearMap.getBinValue(i, j, 1)) < 0.0000000001)
            {
                m_maskValues.setBinValue(i, j, 0, 0.);
                count++;
            }
            else
            {
                m_maskValues.setBinValue(i, j, 0, 1.);
                count2++;
            }
        }
    }

    logger.info() << "sum of the mask: " << mask.getFlux(0);
    logger.info() << "number of zeros: " << count;
    logger.info() << "number of ones: " << count2;
    logger.info() << "over the number of pixels: "
            << shearMap.getXdim() * shearMap.getYdim() * shearMap.getZdim();
}

void InpaintingAlgo::performInPaintingAlgo(ConvergenceMap& convergenceMap)
{
    unsigned int nbIter = m_cartesianParam.getNInpaint();
    bool sigmaBounds = m_cartesianParam.isEqualVarPerScale();
    bool bModeZeros = m_cartesianParam.isForceBMode();

    logger.info() << "number of inpainting iteration: " << nbIter;
    logger.info() << "status of Equal variance per scale is: " << sigmaBounds;
    logger.info() << "status of ForceBMode: " << bModeZeros;

    double maxThreshold(m_maxThreshold);
    double minThreshold(m_minThreshold);

    ConvergenceMap DCTconvergenceMap(convergenceMap.getXdim(),
            convergenceMap.getYdim(), 2);

    logger.info() << "convergenceMap size: " << convergenceMap.getXdim() << " "
            << convergenceMap.getYdim() << " " << convergenceMap.getZdim();
    logger.info() << "DCTconvergenceMap size: " << DCTconvergenceMap.getXdim()
            << " " << DCTconvergenceMap.getYdim() << " "
            << DCTconvergenceMap.getZdim();

    for (size_t iter = 0; iter < nbIter; iter++)
    {
        logger.info() << "iteration " << iter << " beginning";
        logger.info() << "Initializing kappa map";

        // Perform the DCT
        m_MP.performDCT(convergenceMap.getImageAddress(0),
                DCTconvergenceMap.getImageAddress(0));
        m_MP.performDCT(convergenceMap.getImageAddress(1),
                DCTconvergenceMap.getImageAddress(1));

        // Update the threshold value with the max value (kappaE) at first iteration
        // Unused: maxThreshold is fixed at 0.5
        if (iter == 0 && maxThreshold <= 0.)
        {
            maxThreshold = DCTconvergenceMap.getMax(0);
        }

        double lambda = minThreshold
                + (maxThreshold - minThreshold) * (erfc(2.8 * iter / nbIter));

        if (lambda < minThreshold || iter == nbIter - 1)
        {
            lambda = minThreshold;
        }
        logger.info() << "threshold: " << lambda;

        // Keep in memory the 0 values of kappa
        double DCTkappaE0 = DCTconvergenceMap.getBinValue(0, 0, 0);
        double DCTkappaB0 = DCTconvergenceMap.getBinValue(0, 0, 1);

        // Cut all values below the threshold value
        DCTconvergenceMap.applyThreshold(lambda, 0);
        DCTconvergenceMap.applyThreshold(lambda, 1);

        // Restore 0 values
        DCTconvergenceMap.setBinValue(0, 0, 0, DCTkappaE0);
        DCTconvergenceMap.setBinValue(0, 0, 1, DCTkappaB0);

        // Perform the IDCT
        m_MP.performIDCT(DCTconvergenceMap.getImageAddress(0),
                convergenceMap.getImageAddress(0));
        m_MP.performIDCT(DCTconvergenceMap.getImageAddress(1),
                convergenceMap.getImageAddress(1));

        // Apply sigma boundaries (inplace on kappaE)
        if (sigmaBounds)
        {
            applyBoundariesOnWavelets(convergenceMap, 0);
        }

        // Perform the inversion and apply the mask to get the final convergence map
        performInversionMask(convergenceMap, bModeZeros);
        logger.info() << "end of iteration " << iter;

    } // end loop on iterations
}

void InpaintingAlgo::applyBoundariesOnWavelets(GenericMap& input, int k)
{
    std::vector<Matrix> scalesValues;

    m_MP.transformBspline(input[k], scalesValues, nbScales);

    for (int kScale = 0; kScale < nbScales - 1; kScale++)
    {
        double maskSigma = 0.0, imSigma = 0.0, maskCount = 0.0, imCount = 0.0;
        getMaskImageSigma(scalesValues[kScale], maskSigma, imSigma, maskCount,
                imCount);
        multiplyWaveletCoeff(scalesValues[kScale], maskSigma, imSigma,
                maskCount, imCount);
    }
    m_MP.reconsBspline(scalesValues, input[k]);
}

void InpaintingAlgo::getMaskImageSigma(Matrix& kScale, double& maskSigma,
        double& imSigma, double& maskCount, double& imCount)
{

    double maskMean = 0.;
    double maskSquareMean = 0.;
    double imMean = 0.;
    double imSquareMean = 0.;
    for (size_t i = 0; i < m_shearMap.getXdim(); i++)
    {
        for (size_t j = 0; j < m_shearMap.getYdim(); j++)
        {
            double tmp = kScale(i, j);

            if (m_maskValues.getBinValue(i, j, 0) == 0)
            {
                maskMean += tmp;
                maskSquareMean += tmp * tmp;
                maskCount++;
            }
            else
            {
                imMean += tmp;
                imSquareMean += tmp * tmp;
                imCount++;
            }
        }
    }
    maskSigma = sqrt(
            (maskSquareMean / maskCount)
                    - ((maskMean / maskCount) * (maskMean / maskCount)));
    imSigma = sqrt(
            (imSquareMean / imCount)
                    - ((imMean / imCount) * (imMean / imCount)));
}

void InpaintingAlgo::multiplyWaveletCoeff(Matrix& kScale, double& maskSigma,
        double& imSigma, double& maskCount, double& imCount)
{
    if (maskSigma > imSigma * (1 + sqrt(sqrt(2. / (maskCount + 1)))))
    {
        for (size_t i = 0; i < m_shearMap.getXdim(); i++)
        {
            for (size_t j = 0; j < m_shearMap.getYdim(); j++)
            {
                if (m_maskValues.getBinValue(i, j, 0) == 0)
                {
                    if ((maskCount > 9) && (imCount > 9) && (maskSigma > 0))
                    {
                        kScale(i, j) *= imSigma / maskSigma;
                    }
                }
            }
        }
    }
}

void InpaintingAlgo::performInversionMask(ConvergenceMap& convergenceMap,
        bool bModeZeros)
{
    // Force kappaB to zero if required
    if (bModeZeros)
    {
        for (size_t i = 0; i < m_shearMap.getXdim(); i++)
        {
            for (size_t j = 0; j < m_shearMap.getYdim(); j++)
            {
                if (m_maskValues.getBinValue(i, j, 0) == 0)
                {
                    convergenceMap.setBinValue(i, j, 1, 0);
                }
            }
        }
    }

    // Create working shear map
    ShearMap workingShearMap(m_shearMap.getXdim(), m_shearMap.getYdim(), 2);

    // Get the shear map from this convergence map
    convergenceMap.getShearMap(workingShearMap);

    // Apply the mask on the shear map and perform inpainting on it (inplace)
    logger.info() << "applying mask ";

    for (size_t i = 0; i < m_shearMap.getXdim(); i++)
    {
        for (size_t j = 0; j < m_shearMap.getYdim(); j++)
        {
            double m = m_maskValues.getBinValue(i, j, 0);
            double gammaCorr0 = workingShearMap.getBinValue(i, j, 0) * (1 - m)
                    + m_shearMap.getBinValue(i, j, 0) * m;
            double gammaCorr1 = workingShearMap.getBinValue(i, j, 1) * (1 - m)
                    + m_shearMap.getBinValue(i, j, 1) * m;
            workingShearMap.setBinValue(i, j, 0, gammaCorr0);
            workingShearMap.setBinValue(i, j, 1, gammaCorr1);
        }
    }
    workingShearMap.getConvMap(convergenceMap);
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN
