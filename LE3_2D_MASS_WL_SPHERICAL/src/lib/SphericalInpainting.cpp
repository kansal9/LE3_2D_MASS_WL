/**
 * @file src/lib/SphericalInpainting.cpp
 * @date 02/01/21
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

#include "LE3_2D_MASS_WL_SPHERICAL/SphericalInpainting.h"

using namespace Euclid::WeakLensing::TwoDMass::Spherical;

static Elements::Logging logger = Elements::Logging::getLogger(
        "Spherical Inpainting");

namespace LE3_2D_MASS_WL_SPHERICAL
{

SphericalInpainting::SphericalInpainting(
        std::pair<Healpix_Map<double>, Healpix_Map<double>>& shear,
        std::pair<Healpix_Map<double>, Healpix_Map<double>>& convergence,
        SphericalParam &SphericalParam) :
        m_SphParam(SphericalParam), ShearE(shear.first), ShearB(shear.second),
        KappaE(convergence.first), KappaB(convergence.second)
{

    m_nside = ShearE.Nside();
    m_nlmax = 3 * m_nside - 1;
    m_npix = ShearE.Npix();
    m_order = ShearE.Order();
    m_minThreshold = 0.;
    m_maxThreshold = 0.;
    logger.info() << "lmax: " << m_nlmax;
    logger.info() << "npix: " << m_npix;
    m_nbScales = m_SphParam.getNInpScale();
    if (m_nbScales < 2)
    {
        m_nbScales = int(log10(m_nside) / log10(2.));
    }
    logger.info() << "nbScales: " << m_nbScales;
    unsigned int count(0);
    Mask.SetNside(m_nside, RING);
    Mask.fill(1);
    for (int it = 0; it < m_npix; it++)
    {
        if (fabs(ShearE[it]) == 0. && fabs(ShearB[it]) == 0.)
        {
            Mask[it] = 0;
            count++;
        }
    }
    logger.info() << "number of zeros: " << count;
    logger.info() << "over the number of pixels: " << m_npix;

}

std::pair<Healpix_Map<double>, Healpix_Map<double> > SphericalInpainting::maskInversion(
        Healpix_Map<double>& kE, Healpix_Map<double>& kB)
{
    bool bModeZeros = m_SphParam.isForceBMode();
    Healpix_Map<double> mapE;
    Healpix_Map<double> mapB;
    mapE.SetNside(m_nside, RING);
    mapE.fill(0.);
    mapB.SetNside(m_nside, RING);
    mapB.fill(0.);
    for (int it = 0; it < m_npix; it++)
    {
        mapE[it] = kE[it];
        if ((bModeZeros == true) && (Mask[it] == 0))
        {
            mapB[it] = 0.;
        }
        else
        {
            mapB[it] = kB[it];
        }
    }

    SphMassMapping mapping;
    auto gammaPair = mapping.create_ConvtoShearMap(mapE, mapB);
    Healpix_Map<double> gamma1, gamma2;
    gamma1.SetNside(m_nside, RING);
    gamma1.fill(0.);
    gamma2.SetNside(m_nside, RING);
    gamma2.fill(0.);
    for (int it = 0; it < m_npix; it++)
    {
        gamma1[it] = gammaPair.first[it] * (1 - Mask[it])
                + ShearE[it] * Mask[it];
        gamma2[it] = gammaPair.second[it] * (1 - Mask[it])
                + ShearB[it] * Mask[it];
    }
    auto kappaPair = mapping.create_SheartoConvMap(gamma1, gamma2);
    return kappaPair;
}

void SphericalInpainting::applyConstraintsOnWavelets(Healpix_Map<double>& map)
{
    double meanMask = 0.;
    double squareMeanMask = 0.;
    double maskCount = 0.;
    double meanImage = 0.;
    double squareMeanImage = 0.;
    double imageCount = 0.;
    for (int it = 0; it < m_npix; it++)
    {
        if (Mask[it] == 0)
        {
            double tmp = map[it];
            meanMask += tmp;
            squareMeanMask += tmp * tmp;
            maskCount++;
        }
        else
        {
            double tmp = map[it];
            meanImage += tmp;
            squareMeanImage += tmp * tmp;
            imageCount++;
        }
    }
    double maskSigma = sqrt(
            (squareMeanMask / maskCount)
                    - ((meanMask / maskCount) * (meanMask / maskCount)));
    double imSigma = sqrt(
            (squareMeanImage / imageCount)
                    - ((meanImage / imageCount) * (meanImage / imageCount)));
    logger.info() << "imSigma: " << imSigma;
    logger.info() << "maskSigma: " << maskSigma;
    if (maskSigma > imSigma)
    {
        for (int i = 0; i < m_npix; i++)
        {
            if ((Mask[i] == 0))
            {
                if ((maskCount > 9) && (imageCount > 9) && (maskSigma > 0))
                {
                    map[i] *= imSigma / maskSigma;
                }
            }
        }
    }
}

Healpix_Map<double> SphericalInpainting::performWavelet(
        Healpix_Map<double>& map)
{

    arr<double> weight;
    weight.alloc(3 * m_nside - 1);
    weight.fill(1.);

    Alm<xcomplex<double> > map_lm(m_nlmax, m_nlmax);
    map_lm.SetToZero();

    map2alm_iter(map, map_lm, 3, weight);

    Healpix_Map<double> Result;
    Healpix_Map<double> Resi;
    Healpix_Map<double> Band;
    Result.SetNside(m_nside, RING);
    Resi.SetNside(m_nside, RING);
    Band.SetNside(m_nside, RING);
    Result.fill(0.);
    Resi.fill(0.);
    Band.fill(0.);

    for (int it = 0; it < m_npix; it++)
    {
        Result[it] = map[it];
    }

    for (int b = 0; b < m_nbScales; b++)
    {

        double lc = m_nlmax * pow(0.5, b);
        std::vector<double> filter(m_nlmax + 1, 0);
        getFilter(filter, lc, m_nlmax);

        if (b == m_nbScales - 1)
        {
            for (int it = 0; it < m_npix; it++)
            {
                Band[it] = Result[it];
            }
        }
        else
        {
            Alm<xcomplex<double> > Band_lm(m_nlmax, m_nlmax);
            Band_lm.SetToZero();

            for (int l=0; l<=m_nlmax; ++l)
            {
                for (int m=0; m<=l; ++m)
                {
                    Band_lm(l,m) = filter[l] * map_lm(l,m);
                }
            }
            alm2map(Band_lm, Band);

            for (int it = 0; it<m_npix; it++)
            {
                double val = Band[it];
                Band[it] = Result[it] - val;
                Result[it] = val;
            }
            applyConstraintsOnWavelets(Band);
        } // end of if-else
        for (int it = 0; it < m_npix; it++)
        {
            Resi[it] += Band[it];
        }
    } // end of for (nbSclaes)

    for (int it = 0; it < m_npix; it++)
    {
        Result[it] = Resi[it];
    }
    weight.dealloc();
    return Result;
}

std::pair<Healpix_Map<double>, Healpix_Map<double> > SphericalInpainting::performInpainting()
{

    int nbIter = m_SphParam.getNInpaint();
    bool sigmaBounds = m_SphParam.isEqualVarPerScale();
    //logger.info()<<"sigmaBounds: "<<sigmaBounds;
    double maxThreshold(m_maxThreshold);
    double minThreshold(m_minThreshold);

    Healpix_Map<double> outKappaE(m_order, RING);
    Healpix_Map<double> outKappaB(m_order, RING);
    outKappaE.SetNside(m_nside, RING);
    outKappaE.fill(0.);
    outKappaB.SetNside(m_nside, RING);
    outKappaB.fill(0.);

    outKappaE = KappaE;
    outKappaB = KappaB;
    for (int iter = 0; iter < nbIter; iter++)
    {
        logger.info() << "start of iteration: " << iter;

        arr<double> weightE;
        weightE.alloc(3 * m_nside - 1);
        weightE.fill(1.);

        arr<double> weightB;
        weightB.alloc(3 * m_nside - 1);
        weightB.fill(1.);

        Alm<xcomplex<double> > kE_lm( m_nlmax, m_nlmax);
        Alm<xcomplex<double> > kB_lm( m_nlmax, m_nlmax);
        kE_lm.SetToZero();
        kB_lm.SetToZero();

        map2alm_iter(outKappaE, kE_lm, 3, weightE);
        map2alm_iter(outKappaB, kB_lm, 3, weightB);

        //if (maxThreshold <=0 || iter == 0) {
        if (iter == 0)
        {
            maxThreshold = max_abs_alm(kE_lm, m_nlmax);
        }
        double lambda = minThreshold
                + (maxThreshold - minThreshold) * (erfc(2.8 * iter / nbIter));
        //if (lambda < minThreshold || iter==nbIter-1) {
        if (iter == nbIter - 1)
        {
            lambda = minThreshold;
        }

        std::complex<double> kE_lm_0 = kE_lm(0, 0);
        std::complex<double> kB_lm_0 = kB_lm(0, 0);
        logger.info() << "threshold: " << lambda;

        applyThreshold(kE_lm, lambda, m_nlmax);
        applyThreshold(kB_lm, lambda, m_nlmax);

        kE_lm(0, 0) = kE_lm_0;
        kB_lm(0, 0) = kB_lm_0;

        alm2map(kE_lm, outKappaE);
        alm2map(kB_lm, outKappaB);

        if (sigmaBounds)
        {
            outKappaE = performWavelet(outKappaE);
        }

        auto outMapPair = maskInversion(outKappaE, outKappaB);
        //for (int ipix = 0; ipix<m_npix; ipix++) {
        //    outKappaE[ipix] = outMapPair.first[ipix];
        //   outKappaB[ipix] = outMapPair.second[ipix];
        //}
        outKappaE = outMapPair.first;
        outKappaB = outMapPair.second;

        logger.info() << "end of iteration: " << iter;

        weightE.dealloc();
        weightB.dealloc();
    } //end of for loop => max iter for inpainting

    return std::pair<Healpix_Map<double>, Healpix_Map<double> >(outKappaE,
            outKappaB);
}

}  // namespace LE3_2D_MASS_WL_SPHERICAL
