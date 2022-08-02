/**
 * @file src/lib/CartesianParam.cpp
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

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include <cmath>

using namespace dpd::le3::wl::twodmass::inp::paramsconvergencepatch;
using namespace dpd::le3::wl::twodmass::inp::paramsconvergenceclusters;
using namespace dpd::le3::wl::twodmass::inp::paramsconvergencempatchestosphere;
using namespace pro::le3::wl::twodmass;

static Elements::Logging logger = Elements::Logging::getLogger(
        "CartesianParam");

using namespace LE3_2D_MASS_WL_UTILITIES;

namespace LE3_2D_MASS_WL_CARTESIAN
{

/**
 * @brief   main parser function for cartesian parameters
 */
void readParameterFile(const fs::path& ParamFile, CartesianParam &params,
                       CatalogData& shearCat,
                       double catRamin, double catRamax,
                       double catDecmin, double catDecmax,
                       const CatalogData& clusterCat)
{
    logger.info() << "Using file: " << ParamFile << " as DM input product";
    if (!checkFileType(ParamFile, signXML))
    {
        logger.error() << "File is not in xml format";
        return;
    }
    if (getXmlProductType(ParamFile) == "DpdTwoDMassParamsConvergenceClusters")
    {
        logger.info() << "Parameter file: Cluster Convergence Patches";
        params.readConvClustersXMLFile(ParamFile, clusterCat);
        logger.info() << "Done reading Convergence Clusters parameter";
    }
    else if (getXmlProductType(ParamFile) == "DpdTwoDMassParamsConvergencePatch")
    {
        logger.info() << "Parameter file: Convergence Patch";
        params.readConvPatchXMLFile(ParamFile, shearCat);
        logger.info() << "Done reading Convergence Patch parameter";
    }
    else if (getXmlProductType(ParamFile) ==
                                  "DpdTwoDMassParamsConvergencePatchesToSphere")
    {
        logger.info() << "Parameter file: Convergence Patches To Sphere";
        params.readConvPatchesToSphereXMLFile(ParamFile, catRamin, catRamax,
                                              catDecmin, catDecmax, shearCat);
        logger.info() << "Done reading Convergence Patches To Sphere parameter";
    }
}

/**
 * @brief   constructor (to create an object with default values)
 */
CartesianParam::CartesianParam()
{

}

/**
 * @brief   function to read Convergence Patches parameter XML file
 */
void CartesianParam::readConvPatchXMLFile(const fs::path& filepath,
                                          CatalogData& cat)
{
    try
    {
        auto param_xml = DpdTwoDMassParamsConvergencePatch(
                         filepath.native(), xml_schema::flags::dont_validate);

        auto& data_xml = param_xml->Data();

        readBaseParams<twoDMassParamsConvergencePatch>(data_xml);

        setExtName(std::string("KAPPA_PATCHES"));
        setParaFileType(std::string("Conv_Patch"));
        setAddBorder(data_xml.GapsParams().AddBorder());

        // redshift parameters: fill zmin and zmax vectors
        readRedshiftBaseParams<twoDMassParamsConvergencePatch>(data_xml, cat);

        // patch parameters
        auto& pp_xml = data_xml.PatchParams();
        m_Patches.clear();

        // iterate over patch list
        for (auto i = pp_xml.begin(); i != pp_xml.end(); ++i)
        {
            auto& pl(i->PatchList());
            for (auto ii = pl.begin(); ii != pl.end(); ++ii)
            {
                patchDefinition& ptr(*ii);
                m_Patches.emplace_back(ptr.ProjCtr().Longitude() * deg2rad,
                                       ptr.ProjCtr().Latitude() * deg2rad,
                                       ptr.PatchWidth() * deg2rad,
                                       ptr.PixelSize() / 60 * deg2rad);
            }
        }
        setNPatches(m_Patches.size());
    } catch (const xml_schema::exception& e)
    {
        std::cerr << e << std::endl;
    }
}

/**
 * @brief   function to read the convergence Clusters parameter file
 */
void CartesianParam::readConvClustersXMLFile(const fs::path& filepath,
                                             const CatalogData& clusterCat)
{
    try
    {
        auto param_xml = DpdTwoDMassParamsConvergenceClusters(
                         filepath.native(), xml_schema::flags::dont_validate);

        auto& data_xml = param_xml->Data();

        readBaseParams<twoDMassParamsConvergenceClusters>(data_xml);

        setExtName(std::string("KAPPA_PATCHES"));
        setParaFileType(std::string("Conv_Cluster"));
        setAddBorder(data_xml.GapsParams().AddBorder());

        setMassThreshold(data_xml.MassThreshold());
        setZMargin(data_xml.ZMargin());
        setZMaxHalo(data_xml.ZMaxHalo());

        // fill patches and redshift according to clusters (1 bin, but different
        // for each patch: [z+zMargin, zMax])
        setNbins(1);
        double width = data_xml.PatchWidth() * deg2rad;
        double size = data_xml.PixelSize() / 60. * deg2rad;
        double zMaxParam = data_xml.ZMax();
        m_Patches.clear();
        for (long i = 0; i < clusterCat.getNentries(); ++i)
        {
            if((clusterCat["richness"](i) > getMassThreshold())
            && (clusterCat["z"](i) < getZMaxHalo()))
            {
                double raCtr = clusterCat["ra"](i);
                double decCtr = clusterCat["dec"](i);
                m_Patches.emplace_back(raCtr, decCtr, width, size);
                m_ZMin.push_back(clusterCat["z"](i) + getZMargin());
                m_ZMax.push_back(zMaxParam);
            }
        }
        logger.info() << "Number of selected clusters: " << m_Patches.size();
        setNPatches(m_Patches.size());

    } catch (const xml_schema::exception& e)
    {
        std::cerr << e << std::endl;
    }
}

void CartesianParam::readConvPatchesToSphereXMLFile(
        const fs::path& filepath, double& catRamin,
        double& catRamax, double& catDecmin, double& catDecmax,
        CatalogData& cat)
{
    try
    {
        auto param_xml = DpdTwoDMassParamsConvergencePatchesToSphere(
                filepath.native(), xml_schema::flags::dont_validate);
        auto& data_xml = param_xml->Data();

        readBaseParams<twoDMassParamsConvergencePatchesToSphere>(data_xml);

        setExtName(std::string("KAPPA_PATCHES"));
        setParaFileType(std::string("Conv_PatchesToSphere"));
        setAddBorder(data_xml.GapsParams().AddBorder());

        // redshift parameters: fill zmin and zmax vectors
        readRedshiftBaseParams<twoDMassParamsConvergencePatchesToSphere>(data_xml, cat);

        // multipatch parameters
        setNside(data_xml.MultiPatchParams().Nside());
        setOverlap(data_xml.MultiPatchParams().Overlap());

        double width = data_xml.MultiPatchParams().PatchWidth() * deg2rad;
        double size = data_xml.MultiPatchParams().PixelSize() / 60. * deg2rad;
        double d_overlap = 1 + getOverlap();
        double ddec = catDecmax - catDecmin;
        double dra = catRamax - catRamin;

        int nlat = int((ddec * d_overlap) / width) + 1;
        double dtheta_lat = ddec / nlat;

        m_Patches.clear();
        for (int i = 0; i < nlat; i++)
        {
            double lat = i * dtheta_lat + dtheta_lat / 2. + catDecmin;

            // number of patches in longitude wrt overlap and size of the survey
            int nlon = int((dra * cos((i+0.5)*width) * d_overlap) / width) + 1;
            double dtheta_lon = (dra * cos((i+0.5)*width)) / nlon;
            for (int j = 0; j < nlon; j++)
            {
                double lon = j * dtheta_lon + dtheta_lon / 2. + catRamin;
                m_Patches.emplace_back(lon, lat, width, size);
            }
        }
        setNPatches(m_Patches.size());
        logger.info() << "Number of Patches: " << getNPatches();

    } catch (const xml_schema::exception& e)
    {
        std::cerr << "XML ERROR: " << e << std::endl;
    }
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN
