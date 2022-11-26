/**
 * @file LE3_2D_MASS_WL_UTILITIES/VisibilityMask.h
 * @date 06/16/21
 * @author Vanshika Kansal
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

#ifndef _LE3_2D_MASS_WL_UTILITIES_VISIBILITYMASK_H
#define _LE3_2D_MASS_WL_UTILITIES_VISIBILITYMASK_H
#include "ElementsKernel/Logging.h"
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <boost/filesystem.hpp>

#include "ElementsKernel/ProgramHeaders.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"

namespace fs = boost::filesystem;

namespace LE3_2D_MASS_WL_UTILITIES
{

/**
 * @class VisibilityMask
 * @brief
 *
 */
class VisibilityMask
{

public:

    /**
     * @brief Destructor
     */
    virtual ~VisibilityMask() = default;

    /**
     * @brief Constructor with <newNside> <int> Output resolution for map
     */
    explicit VisibilityMask(int newNside);

    /**
     * @brief         method check catalog format (fits/xml) and returns the visibility Mask data
     * @param[in]     <workdir> It's the working directory
     * @param[in]     <catalogName> It's the filename of the mask
     * @param[in]     <data> <std::vector<std::vector<double>>> data of mask in vectors (ra, dec & weight)
     */
    void readVisibilityMask(fs::path& workdir, fs::path& catalogName,
            std::vector<std::vector<double> >& data);

    /**
     * @brief         check catalog format (xml) and complaint with Data Model Product
     * @param[in]     <workdir> <boost::filesystem::path> work directory
     * @param[in]     <InFile> <boost::filesystem::path> Input Catalog name
     * @param[in]     <data> <std::vector<std::vector<double>>> data of mask in vectors (ra, dec & weight)
     */
    void getMaskData(fs::path& workdir, fs::path& InFile,
            std::vector<std::vector<double> >& data);

    /**
     * @brief         check catalog format (fits) and change pixel index to ra & dec
     * @param[in]     <filename> <const std::string> Input mask filename with path
     * @param[in]     <data> <std::vector<std::vector<double>>> data of mask in vectors (ra, dec, pixel_index & weight)
     * @param[in]     <newNside> <int> Output resolution for map
     */
    void getMaskData(const std::string& filename,
            std::vector<std::vector<double> >& data);

    /**
     * @brief         check catalog format (fits) and change pixel index to ra & dec
     * @param[in]     <filename> <const std::string> Input mask filename with path
     * @param[in]     <newNside> <int> Output resolution for map
     * @param[out]    std::pair<Healpix_Map<double>, Healpix_Map<double>> Output maps
     */
    std::pair<Healpix_Map<double>, Healpix_Map<double> > changeResolution(
            const std::string& filename);

    /**
     * @brief    The parameter which returns the name of fits input Mask
     */
    std::string getMaskFitsFilename();

private:
    /**
     * @brief    The parameter which stores the input filename
     */
    std::string m_inputFile;

    /**
     * @brief    The parameter which stores NSide for visibility that need to be used
     */
    int m_newNside;

};
// End of VisibilityMask class

} // LE3_2D_MASS_WL_UTILITIES
#endif
