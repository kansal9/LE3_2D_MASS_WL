/**
 * @file src/lib/MatrixProcess.cpp
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

#include "LE3_2D_MASS_WL_CARTESIAN/MatrixProcess.h"
#include "fftw3.h"
#include "math.h"
#include <iostream>

static Elements::Logging logger = Elements::Logging::getLogger("MatrixProcess");

namespace LE3_2D_MASS_WL_CARTESIAN
{

MatrixProcess::MatrixProcess(unsigned int sizeXaxis, unsigned int sizeYaxis) :
        m_sizeXaxis(sizeXaxis), m_sizeYaxis(sizeYaxis)
{
}

void MatrixProcess::performDCT(double* input, double* output)
{
    // Create the plan on which to perform the transform
    fftw_plan DCTplan;

#pragma omp critical
    {
        DCTplan = fftw_plan_r2r_2d(m_sizeXaxis, m_sizeYaxis, input, output,
                FFTW_REDFT10, FFTW_REDFT10,
                FFTW_ESTIMATE);
    }

    logger.debug() << "Will execute DCTplan with sizes: " << m_sizeXaxis << " "
            << m_sizeYaxis;

    // Perform the transformation
    fftw_execute(DCTplan);

    // Free memory
    fftw_destroy_plan(DCTplan);

    // Rescale the output
    double dctFactor = 2 * sqrt(m_sizeXaxis * m_sizeYaxis);
    for (int i = 0; i < m_sizeXaxis * m_sizeYaxis; i++)
    {
        *(output + i) *= 1. / dctFactor;
    }
}

void MatrixProcess::performIDCT(double* input, double* output)
{

    // Create the plan on which to perform the inverse transform
    fftw_plan IDCTplan;

#pragma omp critical
    {
        IDCTplan = fftw_plan_r2r_2d(m_sizeXaxis, m_sizeYaxis, input, output,
                FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
    }

    logger.debug() << "Will execute IDCTplan with sizes: " << m_sizeXaxis << " "
            << m_sizeYaxis;

    // Perform the transformation
    fftw_execute(IDCTplan);

    // Free memory
    fftw_destroy_plan(IDCTplan);

    // Rescale the output
    double dctFactor = 2 * sqrt(m_sizeXaxis * m_sizeYaxis);
    for (int i = 0; i < m_sizeXaxis * m_sizeYaxis; i++)
    {
        *(output + i) *= 1. / dctFactor;
    }
}

void MatrixProcess::transformBspline(Matrix& input,
        std::vector<Matrix>& scalesValues, unsigned int nbScales)
{
    // Init vector of matrices
    scalesValues.resize(nbScales);

    for (size_t i = 0; i < nbScales; i++)
    {
        scalesValues[i].resize(m_sizeXaxis, m_sizeYaxis);
    }

    // Fill first scale
    for (size_t i = 0; i < input.size1(); i++)
    {
        for (size_t j = 0; j < input.size2(); j++)
        {
            scalesValues[0](i, j) = input(i, j);
        }
    }

    // Loop over the number of scales
    for (size_t step = 0; step < nbScales - 1; step++)
    {
        // Apply the b spline transformation with algorithm a trous
        smoothBspline(scalesValues[step], scalesValues[step + 1], step);
        scalesValues[step] -= scalesValues[step + 1];
    }
}

void MatrixProcess::smoothBspline(Matrix& input, Matrix& output,
        unsigned int stepTrou)
{
    int step = 1 << stepTrou;
    float h[5] =
    { 1 / 16., 1 / 4., 3 / 8., 1 / 4., 1 / 16. };
    output.clear();

    for (int i = 0; i < m_sizeXaxis; i++)
    {
        for (int j = 0; j < m_sizeYaxis; j++)
        {
            for (int l = -2; l <= 2; l++)
            {
                for (int c = -2; c <= 2; c++)
                {
                    int x = i + l * step;
                    if (x < 0)
                    {
                        //x += sizeI; // circular boundary condition
                        x = -x - 1; // mirror boundary condition
                    }
                    else if (x >= m_sizeXaxis)
                    {
                        //x -= sizeI;
                        x = 2 * m_sizeXaxis - x - 1;
                    }

                    int y = j + c * step;
                    if (y < 0)
                    {
                        //y += sizeJ;
                        y = -y - 1;
                    }
                    else if (y >= m_sizeYaxis)
                    {
                        //y -= sizeJ;
                        y = 2 * m_sizeYaxis - y - 1;
                    }

                    output(i, j) += h[l + 2] * h[c + 2] * input(x, y);
                }
            }
        }
    }
}

void MatrixProcess::reconsBspline(std::vector<Matrix>& scalesValues,
        Matrix& imageOut)
{
    for (size_t i = 0; i < imageOut.size1(); i++)
    {
        for (size_t j = 0; j < imageOut.size2(); j++)
        {
            double sum = 0;
            for (unsigned int s = 0; s < scalesValues.size(); s++)
            {
                sum += scalesValues[s](i, j);
            }
            imageOut(i, j) = sum;
        }
    }
}

void MatrixProcess::performDFT(double* input, fftw_complex* output,
        double dftFactor)
{
    // Create the plan on which to perform the transform
    fftw_plan DFTplan;

#pragma omp critical
    {
        DFTplan = fftw_plan_dft_r2c_2d(m_sizeXaxis, m_sizeYaxis, input, output,
                FFTW_ESTIMATE);
    }

    // Perform the transformation
    fftw_execute(DFTplan);

    // Free memory
    fftw_destroy_plan(DFTplan);

    // Rescale the output (be careful of the size in he frequency domain)
    for (int i = 0; i < m_sizeXaxis * (m_sizeYaxis / 2 + 1); i++)
    {
        output[i][0] *= 1. / dftFactor;
        output[i][1] *= 1. / dftFactor;
    }
}

void MatrixProcess::performIDFT(fftw_complex* input, double* output,
        double dftFactor)
{

    // Create the plan on which to perform the inverse transform
    fftw_plan IDFTplan;

#pragma omp critical
    {
        IDFTplan = fftw_plan_dft_c2r_2d(m_sizeXaxis, m_sizeYaxis, input, output,
                FFTW_ESTIMATE);
    }

    // Perform the transformation
    fftw_execute(IDFTplan);

    // Free memory
    fftw_destroy_plan(IDFTplan);

    // Rescale the output
    for (int i = 0; i < m_sizeXaxis * m_sizeYaxis; i++)
    {
        output[i] *= 1. / dftFactor;
    }
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN
