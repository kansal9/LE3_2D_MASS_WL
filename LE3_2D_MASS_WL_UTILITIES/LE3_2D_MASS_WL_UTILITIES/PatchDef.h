/**
 * @file LE3_2D_MASS_WL_CARTESIAN/PatchDef.h
 * @date 02/25/20
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_COORDINATEBOUND_H
#define _LE3_2D_MASS_WL_CARTESIAN_COORDINATEBOUND_H

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

namespace LE3_2D_MASS_WL_UTILITIES
{

/**
 * @class PatchDef
 * @brief
 *
 */
class PatchDef
{

public:
    /**
     * @brief Constructor
     * @param[in] raCtr right ascension of the center of the patch
     * @param[in] decCtr declination of the center of the patch
     * @param[in] patchSize width of the patch
     * @param[in] physical size of the pixel in the patch
     */
    PatchDef(double raCtr, double decCtr, double patchWidth, double pixelSize);

    /**
     * @brief Returns the value of the pixel size
     */
    double getPixelSize() const;

    /**
     * @brief Returns the value of the patch width
     */
    double getPatchWidth() const;

    /**
     * @brief Returns the value of the minimum right ascension
     */
    double getRaMin() const;

    /**
     * @brief Returns the value of the maximum right ascension
     */
    double getRaMax() const;

    /**
     * @brief Returns the value of the minimum declination
     */
    double getDecMin() const;

    /**
     * @brief Returns the value of the maximum declination
     */
    double getDecMax() const;

    /**
     * @brief get Xbin
     * @return Xbin
     */
    unsigned int getXbin() const;

    /**
     * @brief get Ybin
     * @return Ybin
     */
    unsigned int getYbin() const;

private:
    double m_raMin, m_raMax, m_decMin, m_decMax, m_patchWidth, m_pixelSize;

};
// End of PatchDef class

}// namespace LE3_2D_MASS_WL_UTILITIES

#endif
