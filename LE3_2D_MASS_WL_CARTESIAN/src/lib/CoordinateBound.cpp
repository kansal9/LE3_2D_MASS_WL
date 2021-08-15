/**
 * @file src/lib/CoordinateBound.cpp
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

#include "LE3_2D_MASS_WL_CARTESIAN/CoordinateBound.h"

namespace LE3_2D_MASS_WL_CARTESIAN {

CoordinateBound::CoordinateBound(double raMin, double raMax, double decMin, double decMax, double zMin, double zMax):
m_raMin(raMin), m_raMax(raMax), m_decMin(decMin), m_decMax(decMax), m_zMin(zMin), m_zMax(zMax)
{
}

double CoordinateBound::getRaMin() const
{
  return m_raMin;
}

double CoordinateBound::getRaMax() const
{
  return m_raMax;
}

double CoordinateBound::getDecMin() const
{
  return m_decMin;
}

double CoordinateBound::getDecMax() const
{
  return m_decMax;
}

double CoordinateBound::getZMin() const
{
  return m_zMin;
}

double CoordinateBound::getZMax() const
{
  return m_zMax;
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN
