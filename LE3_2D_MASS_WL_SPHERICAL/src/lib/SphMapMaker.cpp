/**
 * @file src/lib/SphMapMaker.cpp
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

#include "LE3_2D_MASS_WL_SPHERICAL/SphMapMaker.h"

#include <cmath>

using namespace Euclid::WeakLensing::TwoDMass;

static Elements::Logging logger = Elements::Logging::getLogger(
        "SphericalMapMaker");
namespace LE3_2D_MASS_WL_SPHERICAL
{

SphMapMaker::SphMapMaker(LE3_2D_MASS_WL_SPHERICAL::SphericalParam &SphParam) :
        npix(0), m_SphParam(SphParam)
{
}


void SphMapMaker::create_ShearMap(
        const CatalogData& catData,
        Healpix_Map<double>& g1_hmap,
        Healpix_Map<double>& g2_hmap,
        Healpix_Map<double>& ngal_hp)
{
    int order = hb.nside2order(m_SphParam.getNside());
    g1_hmap.Set(order, RING);
    g2_hmap.Set(order, RING);
    ngal_hp.Set(order, RING);

    g1_hmap.fill(0);
    g2_hmap.fill(0);
    ngal_hp.fill(0.);

    npix = g1_hmap.Npix();
    logger.info() << "npix: " << npix;

    int ngal = catData.getNentries();
    logger.info() << "ngal: " << ngal;

    for (int i = 0; i < ngal; i++)
    {
        double theta, phi;
        double dec = catData["dec"](i);
        if (isnan(dec))
        {
            dec = 0.;
        }
        theta = -M_PI / 180.0 * dec + M_PI * double(0.5);

        double ra = catData["ra"](i);
        if (isnan(ra))
        {
            ra = 0.;
        }
        phi = M_PI / 180.0 * (ra);

        pointing ptg = pointing(theta, phi);
        ptg.normalize();
        auto id_pix = g2_hmap.ang2pix(ptg);
        auto id = g1_hmap.ang2pix(ptg);
        g1_hmap[id] = g1_hmap[id] + catData["g1"](i);
        g2_hmap[id_pix] = g2_hmap[id_pix] + catData["g2"](i);
        ngal_hp[id_pix] = ngal_hp[id_pix] + 1.;
    }
    for (unsigned int id_pix = 0; id_pix < npix; id_pix++)
    {
        if (ngal_hp[id_pix] != 0)
        {
            g1_hmap[id_pix] = g1_hmap[id_pix] / ngal_hp[id_pix];
            g2_hmap[id_pix] = g2_hmap[id_pix] / ngal_hp[id_pix];
        }
    }
}

} // namespace LE3_2D_MASS_WL_SPHERICAL
