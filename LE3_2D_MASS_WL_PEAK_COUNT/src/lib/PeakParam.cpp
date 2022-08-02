/*
 * Copyright (C) 2012-2020 Euclid Science Ground Segment
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
 */

/**
 * @file src/lib/PeakParam.cpp
 * @date 09/17/19
 * @author user
 */

#include "LE3_2D_MASS_WL_PEAK_COUNT/PeakParam.h"
using namespace LE3_2D_MASS_WL_UTILITIES;

using namespace dpd::le3::wl::twodmass::inp::paramspeakcountconvergence;
using namespace dpd::le3::wl::twodmass::inp::paramspeakcountmassap;

static Elements::Logging logger = Elements::Logging::getLogger("PeakParam");

namespace LE3_2D_MASS_WL_PEAK_COUNT
{

PeakParam::PeakParam() : m_NPeakScale(0), m_ApPeakRadius()
{
}

void PeakParam::readWaveletPeakXML(const fs::path& paramPeakFile)
{
    try
    {
        auto PeakConvParam_xml = DpdTwoDMassParamsPeakCatalogConvergence(
               paramPeakFile.native(), xml_schema::flags::dont_validate);
        m_NPeakScale = PeakConvParam_xml->Data().NPeakScale();
    } catch (const xml_schema::exception& e)
    {
        std::cerr << e << std::endl;
    }
}

void PeakParam::readMassApertureXML(const fs::path& paramPeakFile)
{
    try
    {
        auto PeakMassApParam_xml = DpdTwoDMassParamsPeakCatalogMassAperture(
                  paramPeakFile.native(), xml_schema::flags::dont_validate);

        auto& ApPeakRadius_list = PeakMassApParam_xml->Data().ApPeakRadius();
        logger.info() << "Number of ApPeakRadius: " << ApPeakRadius_list.size();
        for (auto &item : ApPeakRadius_list)
        {
            m_ApPeakRadius.push_back(item);
        }
    } catch (const xml_schema::exception& e)
    {
        std::cerr << e << std::endl;
    }
}

int PeakParam::getNPeakScale() const
{
    return m_NPeakScale;
}

std::vector<double> PeakParam::getApPeakRadius() const
{
    return m_ApPeakRadius;
}

void readPeakParamFile(const fs::path& PeakParamFile, PeakParam &params)
{
    logger.debug() << "Parsing file " << PeakParamFile << " ...";
    if (checkFileType(PeakParamFile.native(), signXML))
    {
        if (getXmlProductType(PeakParamFile) ==
                                     "DpdTwoDMassParamsPeakCatalogMassAperture")
        {
            logger.info() << "Parameter file for aperture mass peak count";
            params.readMassApertureXML(PeakParamFile.native());
        }
        else if (getXmlProductType(PeakParamFile) ==
                                      "DpdTwoDMassParamsPeakCatalogConvergence")
        {
            logger.info() << "Parameter file for peak catalog convergence";
            params.readWaveletPeakXML(PeakParamFile.native());
        }
    }
    else
    {
        throw Elements::Exception() << "Parameter file is not in XML format..";
    }
}

} // LE3_2D_MASS_WL_PEAK_COUNT namespace
