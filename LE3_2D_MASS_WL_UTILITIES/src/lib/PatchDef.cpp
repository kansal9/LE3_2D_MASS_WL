/**
 * @file src/lib/PatchDef.cpp
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

#include "LE3_2D_MASS_WL_UTILITIES/PatchDef.h"

namespace LE3_2D_MASS_WL_UTILITIES
{

PatchDef::PatchDef(double raCtr, double decCtr, double patchWidth,
                   double pixelSize) : m_patchWidth(patchWidth),
                   m_pixelSize(pixelSize)
{
    m_raMin = raCtr - patchWidth / 2.;
    m_raMax = m_raMin + patchWidth;
    m_decMin = decCtr - patchWidth / 2.;
    m_decMax = m_decMin + patchWidth;
}

double PatchDef::getPixelSize() const
{
    return m_pixelSize;
}

double PatchDef::getPatchWidth() const
{
    return m_patchWidth;
}

double PatchDef::getRaMin() const
{
    return m_raMin;
}

double PatchDef::getRaMax() const
{
    return m_raMax;
}

double PatchDef::getDecMin() const
{
    return m_decMin;
}

double PatchDef::getDecMax() const
{
    return m_decMax;
}

unsigned int PatchDef::getXbin() const
{
    double xbin = 0;
    if (fabs(getPatchWidth()) > 0 && fabs(getPixelSize()) > 0)
    {
       xbin = int(ceil(getPatchWidth() / getPixelSize()));
    }
    return xbin;
}

unsigned int PatchDef::getYbin() const
{
    double ybin = 0;
    if (fabs(getPatchWidth()) > 0 && fabs(getPixelSize()) > 0)
    {
        ybin = int(ceil(getPatchWidth() / getPixelSize()));
    }
    return ybin;
}


}  // namespace LE3_2D_MASS_WL_UTILITIES
