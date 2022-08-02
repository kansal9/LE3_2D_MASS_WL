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
#include "LE3_2D_MASS_WL_UTILITIES/Matrix.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MatrixProcess.h"
#include <fftw3.h>
#include <utility>

#include "../../LE3_2D_MASS_WL_UTILITIES/LE3_2D_MASS_WL_UTILITIES/PatchDef.h"

using LE3_2D_MASS_WL_UTILITIES::Matrix;

namespace LE3_2D_MASS_WL_CARTESIAN
{

/**
 * @class GenericMap
 * @brief
 *
 */
class GenericMap
{

public:

    /**
     * @brief Default constructor of a Map
     * This constructor builds a generic Map of dimensions 0, 0, 0
     */
    GenericMap();

    /**
     * @brief Constructor of a Map
     * @param[in] array input values to give to the class
     * @param[in] sizeXaxis number of pixels in the X axis
     * @param[in] sizeYaxis number of pixels in the Y axis
     * @param[in] sizeZaxis number of pixels in the Z axis
     * @param[in] nGalaxies number of galaxies in the map
     * @param[in] pixelSize pixelSize
     * This constructor builds a generic Map of dimensions sizeXaxis, sizeYaxis, sizeZaxis
     * from the input data array
     */
    GenericMap(int Xaxis, int Yaxis, int Zaxis, int nbGalaxies = 0,
               double pixelSize = 0.001);

    /**
     * @brief Constructor of a Map
     * @param[in] filename name of a FITS input file
     * This constructor builds a generic Map from provided FITS input file filename.
     * The provided FITS file should have the image as Primary Header in order to be
     * read properly.
     */
    GenericMap(const std::string& filename);

    /**
     * @brief Constructor of a Map
     * @param[in] filename name of a FITS input file
     * This constructor builds a generic Map from provided FITS input file filename.
     * The provided FITS file should have the image as Primary Header in order to be
     * read properly.
     */
    GenericMap(const fs::path& filename);


    /**
     * @brief Copy constructor of a Map
     * @param[in] copyMap the Map to copy
     * @param[in] copyValues if true does copy data from copyMap, if false copy only metadata
     * This copy constructor builds a copy of the input Map
     */
    GenericMap(GenericMap const& copyMap, bool copyValues = true);

    /**
     * @brief Copy data of a Map (all axis)
     * @param[in] copyMap the Map to copy
     */
    void fullCopy(GenericMap const& copyMap);

    /**
     * @brief Copy data of a Map (single axis)
     * @param[in] copyMap the Map to copy
     * @param[in] k index of the axis to copy
     */
    void singleAxisCopy(GenericMap const& copyMap, int k);

    /**
     * @brief Fill a Map
     * @param[in] filename name of a FITS input file
     * @param[in] indexHdu hdu index number that contains the cube (default 0)
     * This fills the Map from provided FITS input file filename.
     * The provided FITS file should have the image as Primary Header in order to be
     * read properly.
     */
    void fillFromFitsFile(const std::string& filename, int indexHdu = 0);

    /**
     * @brief Outputs the mean values for each Z plans of the map
     * @return vector or double of size m_sizeZaxis containing the mean value
     * of the map for each of those Z plans
     * This method returns a vector of the mean values of each of the Z plans as a vector
     */
    std::vector<double> getMeanValues() const;

    /**
     * @brief Removes an offset value to each Z plans of the map
     * @param[in] a vector of size m_sizeZaxis with an offset value to remove in those plans
     * This method removes an offset value to the entire map in each Z plan
     */
    void removeOffset(std::vector<double> offset);

    /**
     * @brief Returns the value of the needed bin in the map
     * @param[in] binx bin on x axis
     * @param[in] biny bin on y axis
     * @param[in] binz bin on z axis
     * @return the value of the bin (double)
     * This method returns the value of the bin (binx, biny, binz) in the map
     */
    double getBinValue(int binx, int biny, int binz) const;

    /**
     * @brief Return a reference on the Matrix corresponding to the given index
     * @param[in] k index of the matrix to be returned
     */
    Matrix& operator[](int k);

    /**
     * @brief Return a reference on the Matrix corresponding to the given index
     * @param[in] k index of the matrix to be returned
     */
    const Matrix& operator[](int k) const;

    /**
     * @brief Return the address of the Matrix corresponding to the given index
     * @param[in] k index of the matrix to return its address
     */
    double* getImageAddress(int k);

    /**
     * @brief Set the value of the bin in the map
     * @param[in] binx bin on x axis
     * @param[in] biny bin on y axis
     * @param[in] binz bin on z axis
     * @param[in] value
     * This method sets the value of the bin (binx, biny, binz) in the map
     */
    void setBinValue(int binx, int biny, int binz, double value);

    /**
     * @brief reset the map (fill with 0)
     * @param[in] k image to reset
     */
    void clear(int k = 0);

    /**
     * @brief get the standard deviation
     * @param[in] k image on which to compute the standard deviation
     */
    double getSigma(int k = 0) const;

    /**
     * @brief get the minimum
     * @param[in] k image to scan to get the minimum
     */
    double getMin(int k = 0) const;

    /**
     * @brief get the maximum
     * @param[in] k image to scan to get the maximum
     */
    double getMax(int k = 0) const;

    /**
     * @brief get the flux (sum)
     * @param[in] k image to scan to get the flux
     */
    double getFlux(int k = 0) const;

    /**
     * @brief sets all values below a threshold to zero
     * @param[in] threshold
     * @param[in] k image to threshold (default 0)
     */
    void applyThreshold(double threshold, int k = 0);

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
    unsigned int getXdim() const;

    /**
     * @brief Returns the value of the Y dimension
     * @return the value of the Y dimension (int)
     * This method returns the value of the Y dimension
     */
    unsigned int getYdim() const;

    /**
     * @brief Returns the value of the Z dimension
     * @return the value of the Z dimension (int)
     * This method returns the value of the Z dimension
     */
    unsigned int getZdim() const;

    /**
     * @brief returns the number of galaxies in the map if known
     * @return the number of galaxies in the map if known, 0 otherwise
     */
    unsigned int getNumberOfGalaxies() const;

    /**
     * @brief sets the number of galaxies in the map
     */
    void setNumberOfGalaxies(unsigned int N);

    /**
     * @brief returns the size of pixel (degrees) in the map if known
     * @return the size of pixel in the map if known, 0 otherwise
     */
    double getPixelSize() const;

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
    void writeMap(const std::string& filename, CartesianParam &cartesianParam);

    /**
     * @brief Saves the map as a FITS file
     * @param[in] filename name of the file where to save the map
     * @param[in] cartesianParam input parameters that used to create map
     * This method saves the map as a FITS file into the provided filename
     **/
    void writeMap(fs::path filename, CartesianParam &cartesianParam);

    /**
     * @brief   writes records to the given header
     * @param[in] hdu (e.g. primary) to which records need to be written
     * @param[in] cartesianParam input parameters that used to create map
     * This method writes records to the given header
     **/
    void writeImageHeader(const Hdu &hdu, CartesianParam &cartesianParam);

    /**
     * @brief fill with the double* array
     */
    void setArray3D(double *array);

    /**
     * @brief returns the double* array containing values of the pixels
     */
    void getArray(double *array) const;

    /**
     * @brief a 3D raster is created using the map values
     * @return raster
     */
    VecRaster<double, 3> createRaster();

    /**
     * @brief Defines operator = for GenericMap
     */
    void operator=(const GenericMap& map);

    /**
     * @brief Wraps indexes to the actual shape of the map
     * @param[in/out] binx index over the X axis
     * @param[in/out] biny index over the Y axis
     * @param[in/out] binz index over the Z axis
     */
    void wrapIndex(int& binx, int& biny, int& binz) const;

    /**
     * @brief Updates the shape of the map and resize it
     * @param[in] sizeXaxis new size for the X axis
     * @param[in] sizeYaxis new size for the Y axis
     * @param[in] sizeZaxis new size for the Z axis
     */
    void updateSizes(int sizeXaxis, int sizeYaxis, int sizeZaxis);

private:
    /**
     *  @brief <m_mapValues>, map values
     */
    std::vector<Matrix> m_mapValues;

    void initMatrices();

protected:
    /**
     *  @brief <sizeXaxis>, Xaxis of map
     */
    int m_sizeXaxis;
    /**
     *  @brief <sizeYaxis>, Yaxis of map
     */
    int m_sizeYaxis;
    /**
     *  @brief <sizeZaxis>, Zaxis of map
     */
    int m_sizeZaxis;
    /**S
     *  @brief <nGalaxies>, number of galaxies in map
     */
    int m_nGalaxies;

    /**
     *  @brief <m_PixelSize>, size of pixel in map
     */
    double m_PixelSize;
};
// End of GetMap class

}// namespace LE3_2D_MASS_WL_CARTESIAN

#endif
