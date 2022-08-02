/**
 * @file src/lib/PeakCatalog.cpp
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

#include "LE3_2D_MASS_WL_PEAK_COUNT/PeakCatalog.h"

using namespace LE3_2D_MASS_WL_UTILITIES;

namespace LE3_2D_MASS_WL_PEAK_COUNT {

PeakCatalog::PeakCatalog()
{
    m_peakCat["RA_OBJ"] = std::vector<double>();
    m_peakCat["DEC_OBJ"] = std::vector<double>();
    m_peakCat["ZMIN"] = std::vector<double>();
    m_peakCat["ZMAX"] = std::vector<double>();
    m_peakCat["THETA"] = std::vector<double>();
    m_peakCat["SNR"] = std::vector<double>();
}

unsigned int PeakCatalog::getNentries() const
{
    return m_peakCat.at("RA_OBJ").size();
}

void PeakCatalog::getEntry(unsigned int i, double& ra, double& dec,
        double& zmin, double& zmax, double& theta, double& snr) const
{
    ra = m_peakCat.at("RA_OBJ")[i];
    dec = m_peakCat.at("DEC_OBJ")[i];
    zmin = m_peakCat.at("ZMIN")[i];
    zmax = m_peakCat.at("ZMAX")[i];
    theta = m_peakCat.at("THETA")[i];
    snr = m_peakCat.at("SNR")[i];
}

void PeakCatalog::addPeak(double ra, double dec, double zmin, double zmax,
        double theta, double snr)
{
    m_peakCat["RA_OBJ"].push_back(ra);
    m_peakCat["DEC_OBJ"].push_back(dec);
    m_peakCat["ZMIN"].push_back(zmin);
    m_peakCat["ZMAX"].push_back(zmax);
    m_peakCat["THETA"].push_back(theta);
    m_peakCat["SNR"].push_back(snr);
}

void PeakCatalog::writePeakCatalog(const std::string filename)
{
    MefFile peakfile(filename, FileMode::Overwrite);
    const auto col1 = generateColumn("RA_OBJ", m_peakCat["RA_OBJ"]);
    const auto col2 = generateColumn("DEC_OBJ", m_peakCat["DEC_OBJ"]);
    const auto col3 = generateColumn("ZMIN", m_peakCat["ZMIN"]);
    const auto col4 = generateColumn("ZMAX", m_peakCat["ZMAX"]);
    const auto col5 = generateColumn("THETA", m_peakCat["THETA"]);
    const auto col6 = generateColumn("SNR", m_peakCat["SNR"]);
    peakfile.assignBintableExt("PEAK_CATALOG", col1, col2, col3, col4, col5,
            col6);
}


}  // namespace LE3_2D_MASS_WL_PEAK_CATALOG



