/**
 * @file LE3_2D_MASS_WL_CARTESIAN/MapMaker.h
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_MAPMAKER_H
#define _LE3_2D_MASS_WL_CARTESIAN_MAPMAKER_H

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include <utility>
#include <vector>

using LE3_2D_MASS_WL_UTILITIES::CatalogData;
using LE3_2D_MASS_WL_UTILITIES::PatchDef;

namespace LE3_2D_MASS_WL_CARTESIAN
{
/**
 * @class MapMaker
 * @brief
 *
 */
class MapMaker
{

public:

    /**
     * @brief Destructor
     */
    virtual ~MapMaker() = default;

public:
    /**
     * @brief constructor
     */
    MapMaker(CatalogData& inData);

    /**
     * @brief  extract Shear Map
     * @param  Catalog Data
     * @param  respective Parameters from Parameter file
     * This method extract the Shear Map from catalog
     * using parameters like min Ra, min Dec, max Ra, max Dec, Pixels values
     * or using parameters like mapSize, mapCenter, pixelSize
     */
    void getShearMap(const PatchDef& patch, ShearMap& outShearMap,
                     double zmin, double zmax);

private:
    /**
     *  @brief <m_inData>, Catalog Data
     */
    CatalogData& inputData;

};
// End of MapMaker class

}// namespace LE3_2D_MASS_WL_CARTESIAN

#endif
