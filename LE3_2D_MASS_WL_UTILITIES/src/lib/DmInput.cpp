/**
 * @file src/lib/DmInput.cpp
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

#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"

namespace fs = boost::filesystem;

// using namespace dpd::le3::wl::inp::clustercatalogs;

static Elements::Logging logger = Elements::Logging::getLogger("DmInput");

namespace LE3_2D_MASS_WL_UTILITIES {

DmInput::DmInput()
{

}

DmInput::~DmInput()
{

}

void DmInput::readXmlFile(const fs::path& filename)
{
    read_xml(filename.native(), m_tree);
    std::string header_tag = m_tree.begin()->first;
    m_catalogClassName = splitString(header_tag, ":").back();

    // Choose proper dataTagName corresponding to catalog type
    std::string dataTagName;
    if(m_catalogClassName == "DpdWLLE2Catalog")
    {
        dataTagName = "LE2Catalog";
    }
    else if(m_catalogClassName == "DpdTwoDMassKSBCatalog"
         or m_catalogClassName == "DpdTwoDMassRegaussCatalog"
         or m_catalogClassName == "DpdTwoDMassMomentsMLCatalog"
         or m_catalogClassName == "DpdTwoDMassLensMCCatalog")
    {
        dataTagName = "ShearCatalog";
    }
    else if(m_catalogClassName == "DpdTwoDMassClusterCatalog")
    {
        dataTagName = "ClusterCatalog";
    }
    else if(m_catalogClassName == "DpdTwoDMassVisibilityMask")
    {
        dataTagName = "VisibilityMaskHealpix";
    }

    // Parse corresponding section
    std::string nodeName = header_tag + ".Data." + dataTagName
            + ".DataContainer.FileName";
    m_catalogFilename = m_tree.get<std::string>(nodeName);
}

std::string DmInput::getCatalogClassName() const
{
    return m_catalogClassName;
}

std::string DmInput::getCatalogFilename() const
{
    return m_catalogFilename;
}

} // namespace LE3_2D_MASS_WL_UTILITIES

