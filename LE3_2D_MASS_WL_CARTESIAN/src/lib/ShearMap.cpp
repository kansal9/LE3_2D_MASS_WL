/**
 * @file src/lib/ShearMap.cpp
 * @date 05/19/20
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

#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"

static Elements::Logging logger = Elements::Logging::getLogger("ShearMap");

using namespace LE3_2D_MASS_WL_UTILITIES;

namespace LE3_2D_MASS_WL_CARTESIAN
{

ShearMap::ShearMap(GenericMap const& copyMap, bool copyValues) :
        GenericMap(copyMap, copyValues)
{
}

void ShearMap::getConvMap(ConvergenceMap& outputConvMap) const
{
    logger.info() << "Will perform getConvMap with sizes: " << m_sizeXaxis
            << " " << m_sizeYaxis;

    double fftFactor = 1.0 / m_sizeXaxis / m_sizeYaxis;

    // Create the complex maps
    fftw_complex *gamma_complex = (fftw_complex *) fftw_malloc(
            sizeof(fftw_complex) * m_sizeXaxis * m_sizeYaxis);
    fftw_complex *fft_gamma_complex = (fftw_complex *) fftw_malloc(
            sizeof(fftw_complex) * m_sizeXaxis * m_sizeYaxis);
    fftw_complex *Psi_complex = (fftw_complex *) fftw_malloc(
            sizeof(fftw_complex) * m_sizeXaxis * m_sizeYaxis);
    fftw_complex *fft_kappa_complex = (fftw_complex *) fftw_malloc(
            sizeof(fftw_complex) * m_sizeXaxis * m_sizeYaxis);
    fftw_complex *kappa_complex = (fftw_complex *) fftw_malloc(
            sizeof(fftw_complex) * m_sizeXaxis * m_sizeYaxis);

    // Create the plans for transformation
    fftw_plan plan_k_backward;
    fftw_plan plan_g_forward;

#pragma omp critical
    {
        // Create the plans for transformation
        plan_g_forward = fftw_plan_dft_2d(m_sizeXaxis, m_sizeYaxis,
                gamma_complex, fft_gamma_complex, FFTW_FORWARD, FFTW_MEASURE);
        plan_k_backward = fftw_plan_dft_2d(m_sizeXaxis, m_sizeYaxis,
                fft_kappa_complex, kappa_complex, FFTW_BACKWARD, FFTW_MEASURE);
    }

    // Fill the complex gamma map with shear map values
    for (int j = 0; j < m_sizeYaxis; j++)
    {
        for (int i = 0; i < m_sizeXaxis; i++)
        {
            gamma_complex[j * m_sizeYaxis + i][0] = getBinValue(i, j, 0);
            gamma_complex[j * m_sizeYaxis + i][1] = getBinValue(i, j, 1);
        }
    }

    // Perform the fourier transform of the complex shear map
    fftw_execute(plan_g_forward);

    // Create the P factor
    for (int j = 0; j < m_sizeYaxis; j++)
    {
        for (int i = 0; i < m_sizeXaxis; i++)
        {
            // f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / n if n is even
            double l1 =
                    double(i) < double(m_sizeXaxis) / 2. ? i : i - m_sizeXaxis;
            double l2 =
                    double(j) < double(m_sizeYaxis) / 2. ? j : j - m_sizeYaxis;
            double denom = l1 * l1 + l2 * l2;
            Psi_complex[j * m_sizeYaxis + i][0] = (l1 * l1 - l2 * l2) / denom;
            Psi_complex[j * m_sizeYaxis + i][1] = -(2. * l1 * l2) / denom;
        }
    }

    Psi_complex[0][0] = 0.0;
    Psi_complex[0][1] = 0.0;

    // Create the complex convergence map in Fourier space
    for (int k = 0; k < m_sizeXaxis * m_sizeYaxis; k++)
    {
        fft_kappa_complex[k][0] = Psi_complex[k][0] * fft_gamma_complex[k][0]
                - Psi_complex[k][1] * fft_gamma_complex[k][1];
        fft_kappa_complex[k][1] = Psi_complex[k][0] * fft_gamma_complex[k][1]
                + Psi_complex[k][1] * fft_gamma_complex[k][0];
    }

    // Perform inverse Fourier transform to get the convergence map
    fftw_execute(plan_k_backward);

    // Fill the convergence map
    for (int j = 0; j < m_sizeYaxis; j++)
    {
        for (int i = 0; i < m_sizeXaxis; i++)
        {
            outputConvMap.setBinValue(i, j, 0,
                    kappa_complex[j * m_sizeYaxis + i][0] * fftFactor);
            outputConvMap.setBinValue(i, j, 1,
                    kappa_complex[j * m_sizeYaxis + i][1] * fftFactor);
        }
    }

    fftw_destroy_plan(plan_g_forward);
    fftw_destroy_plan(plan_k_backward);

    fftw_free(gamma_complex);
    fftw_free(Psi_complex);
    fftw_free(kappa_complex);
    fftw_free(fft_gamma_complex);
    fftw_free(fft_kappa_complex);
}

void ShearMap::correctReducedShear(
        const LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap& inputConvMap,
        bool useKappaB)
{
    // declare complex numbers
    std::complex<double> g(0, 0);
    std::complex<double> kappa(0, 0);
    std::complex<double> gamma(0, 0);

    for (int i = 0; i < m_sizeXaxis; i++)
    {
        for (int j = 0; j < m_sizeYaxis; j++)
        {
            g = std::complex<double>(getBinValue(i, j, 0),
                    getBinValue(i, j, 1));
            if (useKappaB)
            {
                kappa = std::complex<double>(inputConvMap.getBinValue(i, j, 0),
                        inputConvMap.getBinValue(i, j, 1));
            }
            else
            {
                kappa = std::complex<double>(inputConvMap.getBinValue(i, j, 0),
                        0);
            }
            gamma = g * (std::complex<double>(1., 0) - kappa);
            setBinValue(i, j, 0, std::real(gamma));
            setBinValue(i, j, 1, std::imag(g));
        }
    }
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN
