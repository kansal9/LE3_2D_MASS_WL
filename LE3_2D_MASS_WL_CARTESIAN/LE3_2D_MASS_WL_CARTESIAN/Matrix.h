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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_MATRIX_H
#define _LE3_2D_MASS_WL_CARTESIAN_MATRIX_H
#include "boost/multi_array.hpp"
#include "ElementsKernel/Logging.h"
#include <typeinfo>
#include <cstring>
#include <climits>
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <map>
namespace LE3_2D_MASS_WL_CARTESIAN {

/**
 * @class Matrix
 * @brief
 *
 */
class Matrix {

public:

  /**
   * @brief Destructor
   */
  virtual ~Matrix();

 /**
  * @brief Constructor of an Image Matrix object
  * @param[in] sizeXaxis the number of pixels in the X axis
  * @param[in] sizeYaxis the number of pixels in the Y axis
  * @param[in] values (optional) the input values to fill in the image Matrix (has to be a table
  * of double of expected size sizeXaxis*sizeYaxis). If not provided image filled with zeros.
 */
  Matrix (unsigned int sizeXaxis, unsigned int sizeYaxis, double *values=nullptr);
 /**
  * @brief Copy constructor of an object
  * @param[in] copy the Matrix to copy
 */
  Matrix(const Matrix& copy);

 /**
  * @brief Operator = method
  * @param[in] copy the Image to assign
  * @return a copy of the input image
 */
  Matrix& operator= (const Matrix& copy);
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
   * @brief returns Matrix which is the addition of this and the input image Matrix
   * @param[in] image Matrix the image to add to this
   */
  Matrix add(Matrix image);

  /**
   * @brief returns an matrix which is the substraction of this and the input image matrix
   * @param[in] image matrix the image to substract to this
   */
  Matrix substract(Matrix image);

  /**
   * @brief returns an matrix which is the multiplicator of this and the input factor
   * @param[in] factor the factor to multiply to the image matrix
   */
  Matrix multiply(double factor);

  /**
   * @brief get the image value at position (x, y)
   * @param[in] x the pixel position on X axis
   * @param[in] y the pixel position on Y axis
   * @return the value of the pixel at position (x, y)
   */
  double getValue(int x, int y) const;

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
  double getStandardDeviation() const;

  /**
   * @brief set the image value at position (x, y) to value
   * @param[in] x the pixel position on X axis
   * @param[in] y the pixel position on Y axis
   * @param[in] value the value to assign to the pixel at position (x, y)
   */
  void setValue(unsigned int x, unsigned int y, double value);

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
   * @brief resets all values to zero
   */
  void reset();
  /**
   * @brief returns the double* array containing values of the pixels
   * @return a double* array of dimension sizeXaxis*sizeYaxis
   */
  double* getArray() const;

private:
  /**
   *  @brief <m_values>, matrix values
  */
 double *m_values;
  /**
   *  @brief <m_sizeXaxis>, Xaxis of matrix
  */
 unsigned int m_sizeXaxis;
  /**
   *  @brief <m_sizeYaxis>, Yaxis of matrix
  */
 unsigned int m_sizeYaxis;

};  // End of Matrix class

}  // namespace LE3_2D_MASS_WL_CARTESIAN


#endif
