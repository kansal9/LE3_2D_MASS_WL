/**
 * @file LE3_2D_MASS_WL_PEAK_COUNT/GenericPeakCount.h
 * @date 07/04/22
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

#ifndef _LE3_2D_MASS_WL_PEAK_COUNT_GENERICPEAKCOUNT_H
#define _LE3_2D_MASS_WL_PEAK_COUNT_GENERICPEAKCOUNT_H

#include "LE3_2D_MASS_WL_PEAK_COUNT/PeakCatalog.h"

using LE3_2D_MASS_WL_PEAK_COUNT::PeakCatalog;

namespace LE3_2D_MASS_WL_PEAK_COUNT
{

/**
 * @class GenericPeakCount
 * @brief
 *
 */

class GenericPeakCount
{

public:

    /**
     * @brief Destructor
     */
    virtual ~GenericPeakCount() = default;

    const PeakCatalog& getPeakCatalog() {return m_peakCat;}

protected:
    PeakCatalog m_peakCat;

};
// End of GenericPeakCount class

}// namespace LE3_2D_MASS_WL_PEAK_COUNT

#endif // _LE3_2D_MASS_WL_PEAK_COUNT_GENERICPEAKCOUNT_H
