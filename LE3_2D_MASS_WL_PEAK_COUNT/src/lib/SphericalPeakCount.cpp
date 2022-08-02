/**
 * @file src/lib/SphericalPeakCount.cpp
 * @date 06/18/21
 * @author vanshika kansal
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

#include "LE3_2D_MASS_WL_PEAK_COUNT/SphericalPeakCount.h"

static Elements::Logging logger = Elements::Logging::getLogger(
        "Spherical PeakCount");

using namespace Euclid::WeakLensing::TwoDMass;
using namespace Euclid::WeakLensing::TwoDMass::Spherical;
using namespace LE3_2D_MASS_WL_UTILITIES;

namespace LE3_2D_MASS_WL_PEAK_COUNT
{

SphericalPeakCount::SphericalPeakCount(Healpix_Map<double>& convergenceE,
        int nPeakScale) : m_kappaE(convergenceE), m_nPeakScale(nPeakScale)
{
    m_nside = m_kappaE.Nside();
    m_nlmax = 3 * m_nside - 1;
    m_npix = m_kappaE.Npix();
    m_order = m_kappaE.Order();
    logger.info() << "lmax: " << m_nlmax;
    logger.info() << "npix: " << m_npix;
}

Healpix_Map<double> SphericalPeakCount::getSNRimage_hp(
        Healpix_Map<double>& inputKappa, double globalNoise)
{
    Healpix_Map<double> mySNRimage(m_nside, RING, nside_dummy());
    mySNRimage.fill(0.);
    // divide the signal map by the global noise
    for (int i = 0; i < m_npix; i++)
    {
        mySNRimage[i] = inputKappa[i] * (1. / globalNoise);
    }
    return mySNRimage;
}

void SphericalPeakCount::findPeaks_hp()
{
    double stdev = getStd(m_kappaE);
    auto kappaBands = transformBspline_hp(m_kappaE, m_nPeakScale);

    // loop over scales
    for (unsigned int scale = 0; scale < kappaBands.size() - 1; scale++)
    {
        // Divide the image by the noise, to get an SNR image
        double globalNoise = stdev * SphericalNorm[scale];
        auto snrImage = getSNRimage_hp(kappaBands[scale], globalNoise);
        findPeaksAtTheta(snrImage, (double)scale);
    }
}

void SphericalPeakCount::findPeaksAtTheta(Healpix_Map<double> &map, double theta)
{
    fix_arr<int, 8> neighbours;
    std::pair<double, double> radec;

    for (int ipix = 0; ipix < m_npix; ipix++)
    {
        // fill neighbours array
        map.neighbors(ipix, neighbours);

        // get neighbours maximum value
        double maxNeighboursVal = map[neighbours[0]];
        for (size_t i = 0; i < neighbours.size(); i++)
        {
            if(neighbours[i] < 0)
            {
                continue;
            }
            if(map[neighbours[i]] > maxNeighboursVal)
            {
                maxNeighboursVal = map[neighbours[i]];
            }
        }

        if(map[ipix] > maxNeighboursVal)
        {
            radec = getIndex2RaDec(map, ipix);
            m_peakCat.addPeak(radec.first, radec.second, 0, 0, theta, map[ipix]);
        }
    }
}

}  // namespace LE3_2D_MASS_WL_PEAK_COUNT
