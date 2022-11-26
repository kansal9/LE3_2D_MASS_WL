/**
 * @file src/lib/MapMaker.cpp
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

#include <cmath>

#include "LE3_2D_MASS_WL_CARTESIAN/MapMaker.h"
#include "LE3_2D_MASS_WL_CARTESIAN/Projection.h"

using namespace LE3_2D_MASS_WL_UTILITIES;

static Elements::Logging logger = Elements::Logging::getLogger("MapMaker");

namespace LE3_2D_MASS_WL_CARTESIAN
{

MapMaker::MapMaker(CatalogData& inData) : inputData(inData)
{
}

void MapMaker::getShearMap(const PatchDef& patch, ShearMap& outShearMap,
                           double zmin, double zmax)
{
    logger.info() << "Entering getShearMap function";

    Projection Projection;

    //clear maps
    for (size_t k = 0; k < outShearMap.getZdim(); k++)
    {
        outShearMap.clear(k);
    }

    int xbin = patch.getXbin();
    int ybin = patch.getYbin();
    double halfPatchWidth = patch.getPatchWidth() / 2.;
    double pixelSize = patch.getPixelSize();
    double ra0 = (patch.getRaMax() + patch.getRaMin()) / 2.;
    double dec0 = (patch.getDecMax() + patch.getDecMin()) / 2.;

    unsigned int galCount = 0, selGalCount = 0;
    int tmpx, tmpy;

    for (unsigned int i = 0; i < inputData.getNentries(); i++)
    {
        double g1 = inputData["g1"](i);
        double g2 = inputData["g2"](i);
        double ra = inputData["ra"](i) * deg2rad;
        double dec = inputData["dec"](i) * deg2rad;
        double weight = inputData["w"](i);
        double z = inputData["z"](i) - inputData["z_corr"](i);

        if (z >= zmin && z < zmax)
        {
            // project the selected radec on gnomonic plan
            auto tmpXY = Projection.getGnomonicProjection(ra, dec, ra0, dec0);

            if (abs(tmpXY.first) <= halfPatchWidth
            && abs(tmpXY.second) <= halfPatchWidth)
            {
                //  calculate where it is on the binning
                tmpx = int((tmpXY.first + halfPatchWidth) / pixelSize);
                tmpy = int((tmpXY.second + halfPatchWidth) / pixelSize);

                if (tmpx < 0 || tmpx >= xbin || tmpy < 0 || tmpy >= ybin)
                {
                    continue;
                }

                // Apply correction for the projection
                auto xy1 = Projection.getGnomonicProjection(ra, dec, ra0, dec0);
                auto xy2 = Projection.getGnomonicProjection(ra, (dec + 0.01),
                                                            ra0, dec0);
                double rotationAngle = -atan2((xy2.first - xy1.first),
                        (xy2.second - xy1.second));
                double gamma1cor = -(g1 * cos(2 * rotationAngle)
                                     - g2 * sin(2 * rotationAngle));
                double gamma2cor = -(g1 * sin(2 * rotationAngle)
                                     + g2 * cos(2 * rotationAngle));
                outShearMap[0](tmpx, tmpy) += gamma1cor * weight;
                outShearMap[1](tmpx, tmpy) += gamma2cor * weight;
                outShearMap[2](tmpx, tmpy) += weight;
                selGalCount += weight;
            }
        }
        galCount += weight;
    }
    logger.info() << "number of galaxies selected: " << selGalCount;
    logger.info() << "over the total number of galaxies: " << galCount;

    // Normalize the values to have the mean shear in each bin
    for (int i = 0; i < xbin; i++)
    {
        for (int j = 0; j < ybin; j++)
        {
            if (outShearMap[2](i, j) > 1)
            {
                outShearMap[0](i, j) /= outShearMap[2](i, j);
                outShearMap[1](i, j) /= outShearMap[2](i, j);
            }
        }
    }

    for (size_t k = 0; k < outShearMap.getZdim(); k++)
    {
        logger.info() << "Sum of the map ["<< k <<"]: "<<outShearMap.getFlux(k);
    }

    if (selGalCount == 0)
    {
        logger.info() << "no galaxies are selected for this patch/bin";
    }

    outShearMap.setNumberOfGalaxies(selGalCount);
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN
