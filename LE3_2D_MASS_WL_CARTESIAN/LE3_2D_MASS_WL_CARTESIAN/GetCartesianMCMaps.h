/**
 * @file LE3_2D_MASS_WL_CARTESIAN/GetCartesianMCMaps.h
 * @date 10/13/20
 * @author vkansal
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_GETCARTESIANMCMAPS_H
#define _LE3_2D_MASS_WL_CARTESIAN_GETCARTESIANMCMAPS_H

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GenericMap.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "boost/multi_array.hpp"
#include <utility>
#include <vector>

namespace fs = boost::filesystem;

namespace LE3_2D_MASS_WL_CARTESIAN
{

/**
 * @class GetCartesianMCMaps
 * @brief
 *
 */
class GetCartesianMCMaps
{

public:

    /**
     * @brief Constructor
     */
    GetCartesianMCMaps(CatalogData& inputData, CartesianParam &cartesianParam);

    /**
     * @brief get the Denoised Shear Map
     * @param[out] the denoised shear Map
     */
    void getDeNoisedShearMap(ShearMap& outDenoisedShearMap);

    /**
     * @brief get the Noised Shear Map after randomising the input Shear
     * @param[in] None
     * @return the a vector containing N noised shear Map
     */
    void getNoisedShearMaps(std::vector<ShearMap>& outShearMapList);

    /**
     * @brief perform addition
     * @param[in] DenoisedShearMap
     * @param[in] NoisedShearMap
     * @param[in] filename name of output FITS ShearMap
     * @return true/false
     */
    void performAddition(ShearMap& DenoisedShearMap, ShearMap& NoisedShearMap,
            const std::string& shearMapFilename);

private:
    /**
     *  @brief <m_inData>, Catalog Data
     */
    CatalogData& m_inData;

    /**
     *  @brief <m_cartesianParam>, CartesianParam object with catalog parameters
     */

    CartesianParam m_cartesianParam;

    CartesianAlgoKS m_cartesianAlgo;

}; // End of GetCartesianMCMaps class

}// namespace LE3_2D_MASS_WL_CARTESIAN

#endif
