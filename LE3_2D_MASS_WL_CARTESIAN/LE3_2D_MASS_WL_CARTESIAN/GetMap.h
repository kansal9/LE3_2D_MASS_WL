/**
 * @file LE3_2D_MASS_WL_CARTESIAN/GetMap.h
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_GETMAP_H
#define _LE3_2D_MASS_WL_CARTESIAN_GETMAP_H
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CoordinateBound.h"
#include "boost/multi_array.hpp"
#include <fftw3.h>
#include <utility>
#include "EL_FitsFile/MefFile.h"

namespace LE3_2D_MASS_WL_CARTESIAN {

/**
 * @class GetMap
 * @brief
 *
 */
class GetMap {

public:

 /**
 * @brief Destructor
 */
 virtual ~GetMap();

 /**
 * @brief Constructor of a Map
 * @param[in] array input values to give to the class
 * @param[in] nGalaxies number of galaxies in the map
 *
 * This constructor builds a generic Map of dimensions sizeXaxis, sizeYaxis, sizeZaxis
 * from the input data array
 */
 GetMap(double* array, LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam, int nbGalaxies = 0);

 /**
 * @brief Constructor of a Map
 * @param[in] array input values to give to the class
 * @param[in] sizeXaxis number of pixels in the X axis
 * @param[in] sizeYaxis number of pixels in the Y axis
 * @param[in] sizeZaxis number of pixels in the Z axis
 * @param[in] nGalaxies number of galaxies in the map
 * This constructor builds a generic Map of dimensions sizeXaxis, sizeYaxis, sizeZaxis
 * from the input data array
 */
 GetMap(double* array, int Xaxis, int Yaxis, int Zaxis, int nbGalaxies = 0);

 /**
 * @brief Constructor of a Map
 * @param[in] array input values to give to the class
 * @param[in] sizeXaxis number of pixels in the X axis
 * @param[in] sizeYaxis number of pixels in the Y axis
 * @param[in] sizeZaxis number of pixels in the Z axis
 * @param[in] CoordinateBound contains ra, dec and z information of map
 * @param[in] nGalaxies number of galaxies in the map
 * This constructor builds a generic Map of dimensions sizeXaxis, sizeYaxis, sizeZaxis
 * from the input data array
 */
 GetMap(double* array, int Xaxis, int Yaxis, int Zaxis, LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& coordBound,
        int nbGalaxies = 0);

 /**
 * @brief Constructor of a Map
 * @param[in] filename name of a FITS input file
 * This constructor builds a generic Map from provided FITS input file filename.
 * The provided FITS file should have the image as Primary Header in order to be
 * read properly.
 */
 GetMap(const std::string& filename);

 /**
 * @brief Copy constructor of a Map
 * @param[in] copyMap the Map to copy
 * This copy constructor builds a copy of the input Map
 */
 GetMap(GetMap const& copyMap);

 /**
 * @brief Outputs the mean values for each Z plans of the map
 * @return vector or double of size m_sizeZaxis containing the mean value
 * of the map for each of those Z plans
 * This method returns a vector of the mean values of each of the Z plans as a vector
 */
 std::vector<double> getMeanValues();

 /**
 * @brief Removes an offset value to each Z plans of the map
 * @param[in] a vector of size m_sizeZaxis with an offset value to remove in those plans
 * This method removes an offset value to the entire map in each Z plan
 */
 void removeOffset(std::vector<double> offset);

 /**
 * @brief Applies a gaussian filter of the ConvergenceMap
 * @param sigma sigma of the gaussian kernel
 * This method applies a gaussian filter to the Map of width sigma
 */
 void applyGaussianFilter( float sigma);

 /**
 * @brief Applies a gaussian filter of the ConvergenceMap
 * @param sigmaX sigma of the gaussian kernel on the X direction
 * @param sigmaY sigma of the gaussian kernel on the Y direction
 * This method applies a gaussian filter to the Map of width
 * sigmaX and sigmaY in X and Y directions respectively
 */
 void applyGaussianFilter( float sigmax, float sigmay);

 /**
 * @brief Outputs a gaussian kernel
 * @param[in] sizeX size of the kernel on the X direction
 * @param[in] sizeY size of the kernel on the Y direction
 * @param[in] sigmaX sigma of the gaussian kernel on the X direction
 * @param[in] sigmaY sigma of the gaussian kernel on the Y direction
 * @return a boost multi_array of size {sizeX, sizeY}
 * This method outputs a boost multi_array of size {sizeX, sizeY} which contains a gaussian kernel
 * with sigmaX and sigmaY as sigma in the X and Y directions respectively
 */
 boost::multi_array<double, 2> makeGaussianKernel(long sizeX, long sizeY, double sigmaX, double sigmaY);

 /**
 * @brief Returns the value of the needed bin in the map
 * @param[in] binx bin on x axis
 * @param[in] biny bin on y axis
 * @param[in] binz bin on z axis
 * @return the value of the bin (double)
 * This method returns the value of the bin (binx, biny, binz) in the map
 */
 double getBinValue( int binx, int biny, int binz) const;

  /**
   * @brief get the standard deviation
   */
 double getSigma() const;

  /**
   * @brief sets all values below a threshold to zero
   */
 void thresholding(double value);

 /**
 * @brief Add borders to the map filled with zeros
 */
 void add_borders();

 /**
 * @brief remove borders
 */
 void remove_borders();

 /**
 * @brief Returns the value of the X dimension
 * @return the value of the X dimension (int)
 * This method returns the value of the X dimension
 */
 int getXdim() const;

 /**
 * @brief Returns the value of the Y dimension
 * @return the value of the Y dimension (int)
 * This method returns the value of the Y dimension
 */
 int getYdim() const;

 /**
 * @brief Returns the value of the Z dimension
 * @return the value of the Z dimension (int)
 * This method returns the value of the Z dimension
 */
 int getZdim() const;

 /**
 * @brief returns the coordinate boundaries of the map
 */
 CoordinateBound getCoordinateBound() const;

 /**
 * @brief returns the number of galaxies in the map if known
 * @return the number of galaxies in the map if known, 0 otherwise
 */
 int getNumberOfGalaxies() const;

 /**
 * @brief returns the size of pixel (degrees) in the map if known
 * @return the size of pixel in the map if known, 0 otherwise
 */
 double getPixelSize() const;

 /**
 * @brief Resizes the map to target the galaxiesPerBin number of galaxies per bin
 * @param[in] galaxiesPerBin the target number of galaxies per bin expected
 * @return the binning applied on each dimension, meaning 1 if nothing was done (i.e.
 * the current binning was the closest to the target galaxiesPerBin) or a higher value
 * if binning was applied.
 *
 * This method rebins the map to get a mean number of galaxies per bin the closest
 * to the target input galaxiesPerBin
 *
 */
 unsigned int pixelate(float galaxiesPerBin);

 /**
 * @brief Resizes the map to target the input xBinning and yBinning targets
 * @param[in] xBinning the binning to apply to X dimension as power of two
 * @param[in] yBinning the binning to apply to X dimension as power of two
 * @return true if the binning was applied, false otherwise
 * This method rebins the map with input values xBinning and yBinning
 */
 bool pixelate(int xBinning, int yBinning);

 /**
 * @brief Saves the map as a FITS file
 * @param[in] filename name of the file where to save the map
 * @param[in] cartesianParam input parameters that used to create map
 * This method saves the map as a FITS file into the provided filename
 **/
 void writeMap(const std::string& filename, LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam);

 /**
 * @brief   writes records to the given header
 * @param[in] hdu (e.g. primary) to which records need to be written
 * @param[in] cartesianParam input parameters that used to create map
 * This method writes records to the given header
 **/
 void writeImageHeader(const RecordHdu &hdu, LE3_2D_MASS_WL_CARTESIAN::CartesianParam &cartesianParam);

  /**
   * @brief returns the double* array containing values of the pixels
   * @return a double* array of dimension sizeXaxis*sizeYaxis*sizeZaxis
   */
 void getArray(double *array);

  /**
   * @brief a 3D raster is created using the map values 
   * @return raster
   */
 VecRaster<float, 3> createRaster();

protected:
  /**
   *  @brief <m_mapValues>, map values
  */
  boost::multi_array_ref<double, 3> *m_mapValues;
  //boost::multi_array<double, 3> *m_mapValues;
  /**
   *  @brief <sizeXaxis>, Xaxis of map
  */
  int sizeXaxis;
  /**
   *  @brief <sizeYaxis>, Yaxis of map
  */
  int sizeYaxis;
  /**
   *  @brief <sizeZaxis>, Zaxis of map
  */
  int sizeZaxis;
  /**
   *  @brief <nGalaxies>, number of galaxies in map
  */
  int nGalaxies;

  /**
   *  @brief <m_PixelSize>, size of pixel in map
  */
  double m_PixelSize;
  /**
   *  @brief <m_CB>, CoordinateBound Class object
  */
  CoordinateBound m_CB;

};  // End of GetMap class
}  // namespace LE3_2D_MASS_WL_CARTESIAN

#endif
