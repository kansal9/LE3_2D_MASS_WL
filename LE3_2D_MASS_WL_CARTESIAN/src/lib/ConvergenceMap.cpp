/**
 * @file src/lib/ConvergenceMap.cpp
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

#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"

static Elements::Logging logger = Elements::Logging::getLogger(
        "ConvergenceMap");

using namespace LE3_2D_MASS_WL_UTILITIES;

namespace LE3_2D_MASS_WL_CARTESIAN
{

ConvergenceMap::ConvergenceMap(GenericMap const& copyMap, bool copyValues) :
        GenericMap(copyMap, copyValues)
{
}

void ConvergenceMap::getShearMap(ShearMap& outputShearMap) const
{
    logger.info() << "Will perform getShearMap with sizes: " << m_sizeXaxis
            << " " << m_sizeYaxis;

    double fftFactor = 1.0 / (m_sizeXaxis * m_sizeYaxis);

    // Create the complex maps
    //logger.info() << "Allocate complex arrays";
    fftw_complex* kappa_complex = (fftw_complex *) fftw_malloc(
            sizeof(fftw_complex) * m_sizeXaxis * m_sizeYaxis);
    fftw_complex* fft_kappa_complex = (fftw_complex *) fftw_malloc(
            sizeof(fftw_complex) * m_sizeXaxis * m_sizeYaxis);
    fftw_complex* Psi_complex = (fftw_complex *) fftw_malloc(
            sizeof(fftw_complex) * m_sizeXaxis * m_sizeYaxis);
    fftw_complex* gamma_complex = (fftw_complex *) fftw_malloc(
            sizeof(fftw_complex) * m_sizeXaxis * m_sizeYaxis);
    fftw_complex* fft_gamma_complex = (fftw_complex *) fftw_malloc(
            sizeof(fftw_complex) * m_sizeXaxis * m_sizeYaxis);

    // Create the plans for transformation
    //logger.info() << "Declare DFT plans";
    fftw_plan plan_k_forward;
    fftw_plan plan_g_backward;

    //logger.info() << "Allocate DFT plans";
#pragma omp critical
    {
        plan_k_forward = fftw_plan_dft_2d(m_sizeXaxis, m_sizeYaxis,
                kappa_complex, fft_kappa_complex, FFTW_FORWARD, FFTW_MEASURE);
        plan_g_backward = fftw_plan_dft_2d(m_sizeXaxis, m_sizeYaxis,
                fft_gamma_complex, gamma_complex, FFTW_BACKWARD, FFTW_MEASURE);
    }

    // Fill the complex kappa map with convergence map values
    //logger.info() << "Copy values from convergenceMap to complex array";
    for (int i = 0; i < m_sizeXaxis; i++)
    {
        for (int j = 0; j < m_sizeYaxis; j++)
        {
            kappa_complex[j * m_sizeXaxis + i][0] = getBinValue(i, j, 0);
            kappa_complex[j * m_sizeXaxis + i][1] = getBinValue(i, j, 1);
        }
    }

    // Perform the fourier transform of the complex convergence map
    //logger.info() << "Perform DFT";
    fftw_execute(plan_k_forward);

    // Create the P factor
    //logger.info() << "Compute P factor";
    for (int i = 0; i < m_sizeXaxis; i++)
    {
        for (int j = 0; j < m_sizeYaxis; j++)
        {
            // f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / n if n is even
            double l1 =
                    double(i) < double(m_sizeXaxis) / 2. ? i : i - m_sizeXaxis;
            double l2 =
                    double(j) < double(m_sizeYaxis) / 2. ? j : j - m_sizeYaxis;
            double denom = l1 * l1 + l2 * l2;
            Psi_complex[j * m_sizeXaxis + i][0] = (l1 * l1 - l2 * l2) / denom;
            Psi_complex[j * m_sizeXaxis + i][1] = (2. * l1 * l2) / denom;
        }
    }
    Psi_complex[0][0] = 0.0;
    Psi_complex[0][1] = 0.0;

    // Create the complex shear map in Fourier space
    //logger.info() << "Multiplication in Fourier space";
    for (int k = 0; k < m_sizeXaxis * m_sizeYaxis; k++)
    {
        fft_gamma_complex[k][0] = Psi_complex[k][0] * fft_kappa_complex[k][0]
                - Psi_complex[k][1] * fft_kappa_complex[k][1];
        fft_gamma_complex[k][1] = Psi_complex[k][0] * fft_kappa_complex[k][1]
                + Psi_complex[k][1] * fft_kappa_complex[k][0];
    }

    // Perform inverse Fourier transform to get the shear map
    //logger.info() << "Perform IDFT";
    fftw_execute(plan_g_backward);

    // Fill the shear map
    //logger.info() << "Copy values from complex array to shearMap";
    for (int i = 0; i < m_sizeXaxis; i++)
    {
        for (int j = 0; j < m_sizeYaxis; j++)
        {
            outputShearMap.setBinValue(i, j, 0,
                    gamma_complex[j * m_sizeXaxis + i][0] * fftFactor);
            outputShearMap.setBinValue(i, j, 1,
                    gamma_complex[j * m_sizeXaxis + i][1] * fftFactor);
        }
    }

    fftw_destroy_plan(plan_k_forward);
    fftw_destroy_plan(plan_g_backward);
    fftw_free(gamma_complex);
    fftw_free(Psi_complex);
    fftw_free(kappa_complex);
    fftw_free(fft_gamma_complex);
    fftw_free(fft_kappa_complex);
}

void ConvergenceMap::applyGaussianFilter(float sigma, int k)
{
    applyGaussianFilter(sigma, sigma, k);
}

void ConvergenceMap::applyGaussianFilter(float sigmax, float sigmay, int k)
{
    // Define sizes
    int sizeXaxisFreq = m_sizeXaxis;
    int sizeYaxisFreq = m_sizeYaxis / 2 + 1;

    // Define frequency complex maps
    fftw_complex* gaussianKernalTf = (fftw_complex *) fftw_malloc(
            sizeof(fftw_complex) * sizeXaxisFreq * sizeYaxisFreq);
    fftw_complex* signalTf = (fftw_complex *) fftw_malloc(
            sizeof(fftw_complex) * sizeXaxisFreq * sizeYaxisFreq);

    // Generate a normalized gaussian kernel with sigma
    GenericMap gaussianKernel(m_sizeXaxis, m_sizeYaxis, 1);
    auto gaussianKernalRef = gaussianKernel.getImageAddress(0);
    makeGaussianKernel(gaussianKernalRef, sigmax, sigmay, m_sizeXaxis,
            m_sizeYaxis);

    // Create matrix process
    MatrixProcess mp(m_sizeXaxis, m_sizeYaxis);

    // Get the TF of the kernel (also real gaussian)
    mp.performDFT(gaussianKernalRef, gaussianKernalTf);

    // Will work separately in the k-ith image
    mp.performDFT(getImageAddress(k), signalTf);

    // perform multiplication with the kernel
    int i = 0;
    int sgn = 1;
    for (int x = 0; x < sizeXaxisFreq; x++)
    {
        for (int y = 0; y < sizeYaxisFreq; y++)
        {
            i = x * sizeYaxisFreq + y;
            sgn = pow(-1, x + y); // alternate signs to mimic shift after IDFT
            double real = signalTf[i][0] * gaussianKernalTf[i][0]
                    - signalTf[i][1] * gaussianKernalTf[i][1];
            double imag = signalTf[i][0] * gaussianKernalTf[i][1]
                    + signalTf[i][1] * gaussianKernalTf[i][0];
            signalTf[i][0] = real * sgn;
            signalTf[i][1] = imag * sgn;
        }
    }

    mp.performIDFT(signalTf, getImageAddress(k), m_sizeXaxis * m_sizeYaxis);

    // Free memory
    fftw_free(gaussianKernalTf);
    fftw_free(signalTf);
}

void ConvergenceMap::setKappaBToZero()
{
    clear(1);
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN
