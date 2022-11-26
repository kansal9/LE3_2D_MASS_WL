/**
 * @file src/lib/Matrix.cpp
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

#include "LE3_2D_MASS_WL_UTILITIES/Matrix.h"
#include <iostream>
#include "math.h"
static Elements::Logging logger = Elements::Logging::getLogger("Matrix");
namespace LE3_2D_MASS_WL_UTILITIES
{

unsigned int Matrix::getXdim() const
{
    return size1();
}

unsigned int Matrix::getYdim() const
{
    return size2();
}

double Matrix::getMax() const
{
    return *std::max_element(data().begin(), data().end());
}

double Matrix::getMin() const
{
    return *std::min_element(data().begin(), data().end());
}

double Matrix::getFlux() const
{
    return std::accumulate(data().begin(), data().end(), (double) 0);
}

double Matrix::getSigma() const
{
    unsigned int counts(1);
    double x, m, m_old(0), v(0), v_old(0);
    for (unsigned i = 0; i != size1(); ++i)
    {
        for (unsigned j = 0; j != size2(); ++j)
        {
            x = (*this)(i, j);
            if (counts == 1)
            {
                m_old = m = x;
            }
            else
            {
                m = m_old + (x - m_old) / counts;
                v = v_old + (x - m_old) * (x - m);
                m_old = m;
                v_old = v;
            }
            counts++;
        }
    }

    double stdev = sqrt(v / (counts - 1));
    return stdev;
}

void Matrix::applyThreshold(double threshold)
{
    for (unsigned int i = 0; i < size1(); i++)
    {
        for (unsigned int j = 0; j < size2(); j++)
        {
            if (fabs((*this)(i, j)) < threshold)
            {
                (*this)(i, j) = 0;
            }
        }

    }
}

bool Matrix::isLocalMax(unsigned int x, unsigned int y) const
{
    if (x == 0 || y == 0 || x == size1() - 1 || y == size2() - 1)
    {
        return false;
    }
    double localVal = (*this)(x, y);
    // Loop over all neighboors
    for (int i = -1; i < 2; i++)
    {
        for (int j = -1; j < 2; j++)
        {
            // Do not check the pixel itself
            if (i != 0 || j != 0)
            {
                // If a neighboor is greater or equal then this is not a local max
                if ((*this)(x + i, y + j) >= localVal)
                {
                    return false;
                }
            }
        }
    }
    return true;
}

const double* Matrix::getArray() const
{
    return m_values.data().begin();
}

double* Matrix::getArray()
{
    return const_cast<double*>(m_values.data().begin());
}

}  // namespace LE3_2D_MASS_WL_UTILITIES
