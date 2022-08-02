/**
 * @file src/lib/CartesianProcessor.cpp
 * @date 02/15/22
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

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianProcessor.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

using namespace LE3_2D_MASS_WL_UTILITIES;
using Logging = Elements::Logging;

namespace LE3_2D_MASS_WL_CARTESIAN {

static Logging logger = Logging::getLogger("CartesianProcessor");

CartesianProcessor::CartesianProcessor(
        std::map<std::string, std::string>& args) : m_args(args),
                m_shearCatalogPath(), m_shearMapPath(), m_clusterCatalog(),
                m_produceMcMaps(false)
{
    parseOptions();
}

CartesianProcessor::CartesianProcessor(
        std::map<std::string, variable_value>& args) : m_args(),
                m_shearCatalogPath(), m_shearMapPath(), m_clusterCatalog(),
                m_produceMcMaps(false)
{
    // iterate over map and turn variable_value to string
    for(auto iter = args.begin(); iter != args.end(); ++iter)
    {
        m_args[iter->first] = iter->second.as<std::string>();
    }
    parseOptions();
}

void CartesianProcessor::parseOptions()
{
    m_workdir = m_args["workdir"];
    m_paramfile = m_args["paramfile"];
    m_paramtype = getXmlProductType(m_paramfile);

    // if parameter product type not in available ones, throw an error
    if(ParameterTypesAvail.count(m_paramtype) == 0)
    {
        logger.fatal("Bad parameter type!");
        throw ExitCode::USAGE;
    }

    logger.info("Parameter type: " + m_paramtype);

    if(m_paramtype == "DpdTwoDMassParamsConvergencePatch")
    {
        // Case: small patches defined in parameters, need a shear catalog.
        // May process multiple patches in multiple redshift bins.
        m_shearCatalogPath = m_args["shear"];
    }
    else if(m_paramtype == "DpdTwoDMassParamsConvergenceClusters")
    {
        // Case: small patches defined in cluster catalog, need a shear catalog
        // and a cluster catalog. Will process multiple patches in a single
        // redshift bins.
        m_shearCatalogPath = m_args["shear"];
        m_clusterCatalog = m_args["clusterCatalog"];
    }
    else if(m_paramtype == "DpdTwoDMassParamsConvergencePatchesToSphere")
    {
        // Case: multiple patches covering the sphere, need a shear catalog.
        // Will process multiple patches in multiple redshift bins (TBC, maybe
        // only 1 due to computation time).
        // TODO
    }
    else if(m_paramtype == "DpdTwoDMassParamsPeakCatalogConvergence")
    {
        // Case: TODO
    }
    else if(m_paramtype == "DpdTwoDMassParamsPeakCatalogMassAperture")
    {
        // Case: TODO
    }
    else if(m_paramtype == "DpdTwoDMassParamsConvergenceSphere")
    {
        // Case: TODO
    }
}

void CartesianProcessor::process()
{
    CartesianParam params;
    CatalogData cat, clusterCat;

    // Load shear catalog
    if(not m_shearCatalogPath.empty())
    {
        cat.getCatalogData(m_workdir, m_shearCatalogPath);
        logger.info() << "Number of galaxies in catalog: " << cat.getNentries();
    }

    // Load cluster catalog if exists
    if(not m_clusterCatalog.empty())
    {
        clusterCat.getCatalogData(m_workdir, m_clusterCatalog);
        logger.info() << "Number of clusters in catalog: " << cat.getNentries();
    }

    // Read parameters (and patches definitions)
    readParameterFile(m_paramfile, params, cat, 0, 0, 0, 0, clusterCat);

    // Build cartesian algo
    CartesianAlgoKS cartesianAlgoKS(params);

    // Iterate over patches
    int Npatches = params.getNPatches();
    logger.info() << "Start iteration over NPatches = " << Npatches;
    for(int ip=0; ip<Npatches; ip++)
    {
        auto& patch = params.getPatch(ip);
        logger.info() << "Patch bounds:";
        logger.info() << "Dec: " << patch.getDecMin() << " "
                                 << patch.getDecMax();
        logger.info() << "Ra: " << patch.getRaMin() << " "
                                << patch.getRaMax();
        logger.info() << "Width: " << patch.getPatchWidth();
        logger.info() << "Pixel size: " << patch.getPixelSize();

        // Iterate over redshift bins
        int Nbins = params.getNbins();
        logger.info() << "Start iteration over Nbins of redshift = " << Nbins;
        logger.info() << "----------------------------------------------------";
        for (int iz=0; iz<Nbins; iz++)
        {
            double zmin = params.getZMin(iz);
            double zmax = params.getZMax(iz);
            logger.info() << "Redshift bounds:";
            logger.info() << "zMin: " << zmin;
            logger.info() << "zMax: " << zmax;

            logger.info() << "Create shearMap before filling with data";
            logger.info() << "SizeXaxis: " << patch.getXbin();
            logger.info() << "SizeYaxis: " << patch.getYbin();

            // extract patch into shear map
            ShearMap shearMap(patch.getXbin(), patch.getYbin(), 3);
            cartesianAlgoKS.extractShearMap(shearMap, cat, patch, zmin, zmax);

            if(shearMap.getNumberOfGalaxies() == 0)
            {
                logger.info() << "Nothing to do: continue...";
                logger.info() << "--------------------------------------------";
                continue;
            }

            // saved shear map [test only, not a PF product]
            fs::path outShearMapFilename("ShearMap_extracted_Patch"
                    + std::to_string(ip) + "_Zbin" + std::to_string(iz) + "_"
                    + getDateTimeString() + ".fits");
            std::string name = "SHEAR_PATCH";
            params.setExtName(name);
            shearMap.writeMap(outShearMapFilename, params);

            // perform mass mapping

            logger.info() << "------------------------------------------------";
        }
    }
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN


