/**
 * @file LE3_2D_MASS_WL_CARTESIAN/Matrix.h
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

#ifndef _LE3_2D_MASS_WL_UTILITIES_MATRIX_H
#define _LE3_2D_MASS_WL_UTILITIES_MATRIX_H
#include <boost/numeric/ublas/matrix.hpp>
#include "ElementsKernel/Logging.h"
#include <typeinfo>
#include <cstring>
#include <climits>
#include <cmath>
#include <string>
#include <vector>
#include <numeric>
#include <list>
#include <map>

using namespace boost::numeric::ublas;

namespace LE3_2D_MASS_WL_UTILITIES
{

/**
 * @class Matrix
 * @brief
 *
 */
class Matrix: public matrix<double>
{

public:

    using matrix<double>::matrix;    // use the constructors already defined

    using matrix<double>::operator=; // and the operator= already defined

    /**
     * @brief get the number of pixels on X axis
     * @return the number of pixels on X axis
     */
    unsigned int getXdim() const;

    /**
     * @brief get the number of pixels on Y axis
     * @return the number of pixels on Y axis
     */
    unsigned int getYdim() const;

    /**
     * @brief get the max image value
     */
    double getMax() const;

    /**
     * @brief get the min image value
     */
    double getMin() const;

    /**
     * @brief get the Flux of image
     */
    double getFlux() const;

    /**
     * @brief get the standard deviation in the image
     */
    double getSigma() const;

    /**
     * @brief sets all values below a threshold to zero
     */
    void applyThreshold(double threshold);

    /**
     * @brief returns true if the pixel (x , y) is a local max
     * @param[in] x the position of the pixel on X axis
     * @param[in] y the position of the pixel on Y axis
     * @return true if the pixel (x , y) is a local max, false otherwise
     */
    bool isLocalMax(unsigned int x, unsigned int y) const;

    /**
     * @brief returns the double* const array containing values of the pixels
     * @return a double* array of dimension sizeXaxis*sizeYaxis
     */
    const double* getArray() const;

    /**
     * @brief returns the double* array containing values of the pixels
     * @return a double* array of dimension sizeXaxis*sizeYaxis
     */
    double* getArray();

private:
    /**
     *  @brief <m_values>, matrix values
     */
    matrix<double> m_values;

};
// End of Matrix class

}// namespace LE3_2D_MASS_WL_UTILITIES

#endif
