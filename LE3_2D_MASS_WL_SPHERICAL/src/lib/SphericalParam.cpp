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
 * @file src/lib/SphericalParam.cpp
 * @date 07/14/19
 * @author user
 */

#include "LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h"

// Input namespace and classes
using namespace dpd::le3::wl::twodmass::inp::paramsconvergencesphere;
using namespace dpd::le3::wl::twodmass::inp::paramsconvergencempatchestosphere;
using namespace pro::le3::wl::twodmass;

static Elements::Logging logger = Elements::Logging::getLogger(
        "SphericalParam");

namespace LE3_2D_MASS_WL_SPHERICAL
{

SphericalParam::SphericalParam()
{
}


void SphericalParam::readConvSphereXMLFile(const fs::path& filepath,
                                           CatalogData& cat)
{
    logger.info() << "Getting information from input parameter XML file "
                  << filepath << " ...";
    try
    {
        auto param_xml = DpdTwoDMassParamsConvergenceSphere(
                        filepath.native(), xml_schema::flags::dont_validate);

        auto& data_xml = param_xml->Data();

        readBaseParams<twoDMassParamsConvergenceSphere>(data_xml);
        readRedshiftBaseParams<twoDMassParamsConvergenceSphere>(data_xml, cat);

        setExtName(std::string("KAPPA_SPHERE"));
        setNside(param_xml->Data().Nside());


    } catch (const xml_schema::exception& e)
    {
        std::cerr << e << std::endl;
    }
}

} // LE3_2D_MASS_WL_SPHERICAL namespace
