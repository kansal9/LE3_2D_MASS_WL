/**
 * @file LE3_2D_MASS_WL_PEAK_COUNT/PeakCatalog.h
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

#ifndef _LE3_2D_MASS_WL_PEAK_COUNT_PEAK_CATALOG_H
#define _LE3_2D_MASS_WL_PEAK_COUNT_PEAK_CATALOG_H

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

#include <vector>
#include <map>
#include <string>

namespace LE3_2D_MASS_WL_PEAK_COUNT
{

/**
 * @class PeakCatalog
 * @brief
 *
 */

class PeakCatalog
{

public:

    PeakCatalog();

    unsigned int getNentries() const;

    void getEntry(unsigned int i, double& ra, double& dec,
            double& zmin, double& zmax, double& theta, double& snr) const;

    /**
     * @brief Destructor
     */
    virtual ~PeakCatalog() = default;

    void addPeak(double ra, double dec, double zmin, double zmax, double theta,
            double snr);

    void writePeakCatalog(const std::string filename);

private:
    std::map<std::string, std::vector<double>> m_peakCat;

};
// End of PeakCatalog class

}// namespace LE3_2D_MASS_WL_PEAK_CATALOG

#endif // _LE3_2D_MASS_WL_PEAK_COUNT_PEAK_CATALOG_H
