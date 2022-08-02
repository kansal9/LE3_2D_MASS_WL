/**
 * @file src/lib/NoisyCatalogData.cpp
 * @date 10/13/20
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

#include "LE3_2D_MASS_WL_UTILITIES/NoisyCatalogData.h"
#include <cmath>

static Elements::Logging logger = Elements::Logging::getLogger("RandomizeData");

namespace LE3_2D_MASS_WL_UTILITIES {

NoisyCatalogData::NoisyCatalogData()
{
}

CatalogData NoisyCatalogData::getNoisyCatalog(const CatalogData& inputData)
{
    CatalogData output(inputData);
    long galCount = inputData.getNentries();
    for (unsigned int i = 0; i < galCount; i++)
    {
        double theta, radius, temp;
        temp = pow(inputData["g1"](i), 2) + pow(inputData["g2"](i), 2);
        radius = sqrt(temp);
        theta = M_PI * ((double) rand() / RAND_MAX);
        output["g1"](i) = radius * cos(2 * theta);
        output["g2"](i) = radius * sin(2 * theta);
    }
    return output;
}

}  // namespace LE3_2D_MASS_WL_UTILITIES

