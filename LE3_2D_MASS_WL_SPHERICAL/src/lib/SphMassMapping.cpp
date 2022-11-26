/**
 * @file src/lib/SphMassMapping.cpp
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

#include "LE3_2D_MASS_WL_SPHERICAL/SphMassMapping.h"

#include <cmath>
#include <complex>

#include <iostream>
#include <fstream>

using Euclid::WeakLensing::TwoDMass::Spherical::SphericalIO;
using namespace Euclid::WeakLensing::TwoDMass::Spherical;

static Elements::Logging logger = Elements::Logging::getLogger(
        "SphericalMassMapping");
namespace LE3_2D_MASS_WL_SPHERICAL
{

SphMassMapping::SphMassMapping() :
        m_nside(2048)
{
}

std::pair<Healpix_Map<double>, Healpix_Map<double> > SphMassMapping::create_SheartoConvMap(
        Healpix_Map<double>& mapE, Healpix_Map<double>& mapB)
{
    return sphMassMap(mapType::shearMap, mapE, mapB);
}

std::pair<Healpix_Map<double>, Healpix_Map<double> > SphMassMapping::create_ConvtoShearMap(
        Healpix_Map<double>& mapE, Healpix_Map<double>& mapB)
{
    return sphMassMap(mapType::convMap, mapE, mapB);
}

void SphMassMapping::correctScheme(Healpix_Map<double>* map)
{
    if (map->Scheme() == NEST)
    {
        map->swap_scheme();
    }
}

std::pair<Healpix_Map<double>, Healpix_Map<double> > SphMassMapping::sphMassMap(
        const mapType type,
        Healpix_Map<double>& mapE, Healpix_Map<double>& mapB)
{

    correctScheme(&mapE);
    correctScheme(&mapB);
    m_nside = mapE.Nside();
    int m_nlmax = 3 * m_nside - 1;
    logger.info() << "nside: " << m_nside;

    Healpix_Map<double> T;
    T.SetNside(m_nside, RING);
    T.fill(0.);

    Healpix_Map<double> out_mapE;
    out_mapE.SetNside(m_nside, RING);
    out_mapE.fill(0.);

    Healpix_Map<double> out_mapB;
    out_mapB.SetNside(m_nside, RING);
    out_mapB.fill(0.);

    arr<double> weight;
    weight.alloc(3 * m_nside - 1);
    weight.fill(1.);

    Alm<xcomplex<double> > almT( m_nlmax, m_nlmax);
    Alm<xcomplex<double> > alm_mapE( m_nlmax, m_nlmax);
    Alm<xcomplex<double> > alm_mapB( m_nlmax, m_nlmax);
    almT.SetToZero();
    alm_mapE.SetToZero();
    alm_mapB.SetToZero();

    if (type == mapType::shearMap)
    {
        map2alm_pol_iter(T, mapE, mapB, almT, alm_mapE, alm_mapB, 3, weight);
        // map2alm_pol(T, mapE, mapB, almT, alm_mapE, alm_mapB, weight, false);
        almT.SetToZero();

        for (int l_index = 0; l_index <= m_nlmax; ++l_index)
        {
            for (int m_index = 0; m_index <= l_index; ++m_index)
            {
                if (l_index != 1)
                {
                    double temp = sqrt(
                            ((l_index + 1.0) * l_index)
                                    / ((l_index + 2.) * (l_index - 1.)));
                    if (isnan(temp))
                    {
                        temp = 0;
                    }
                    alm_mapE(l_index, m_index) = alm_mapE(l_index, m_index)
                            * temp;
                    alm_mapB(l_index, m_index) = alm_mapB(l_index, m_index)
                            * temp;
                }
                //if (l_index == 0 || l_index == 1) {
                if (l_index == 1)
                {
                    alm_mapE(l_index, m_index) = 0.0;
                    alm_mapB(l_index, m_index) = 0.0;
                }
            }
        }
        alm2map(alm_mapE, out_mapE);
        alm2map(alm_mapB, out_mapB);
    }
    if (type == mapType::convMap)
    {
        //map2alm(mapE, alm_mapE, weight);
        //map2alm(mapB, alm_mapB, weight);
        map2alm_iter(mapE, alm_mapE, 3, weight);
        map2alm_iter(mapB, alm_mapB, 3, weight);

        for (int l_index = 0; l_index <= m_nlmax; ++l_index)
        {
            for (int m_index = 0; m_index <= l_index; ++m_index)
            {
                if (l_index != 0)
                {
                    double temp = sqrt(
                            ((l_index + 2.) * (l_index - 1.))
                                    / ((l_index + 1.0) * l_index));
                    if (isnan(temp))
                    {
                        temp = 0;
                    }
                    alm_mapE(l_index, m_index) = alm_mapE(l_index, m_index)
                            * temp;
                    alm_mapB(l_index, m_index) = alm_mapB(l_index, m_index)
                            * temp;
                }
                //if (l_index == 0 || l_index == 1) {
                if (l_index == 0)
                {
                    alm_mapE(l_index, m_index) = 0.0;
                    alm_mapB(l_index, m_index) = 0.0;
                }
            }
        }
        alm2map_pol(almT, alm_mapE, alm_mapB, T, out_mapE, out_mapB);
    }
    weight.dealloc();
    return std::pair<Healpix_Map<double>, Healpix_Map<double> >(out_mapE,
            out_mapB);
}

int SphMassMapping::getNside()
{
    return m_nside;
}

void SphMassMapping::computeReducedShear_hp(
        std::pair<Healpix_Map<double>, Healpix_Map<double> >& shearMap,
        std::pair<Healpix_Map<double>, Healpix_Map<double> >& convergenceMap)
{

    int m_npix = shearMap.first.Npix();
    for (int it = 0; it < m_npix; it++)
    {
        shearMap.first[it] = shearMap.first[it]
                * (1 - convergenceMap.first[it]);
        shearMap.second[it] = shearMap.second[it]
                * (1 - convergenceMap.second[it]);
    }

}

}  // namespace LE3_2D_MASS_WL_SPHERICAL
