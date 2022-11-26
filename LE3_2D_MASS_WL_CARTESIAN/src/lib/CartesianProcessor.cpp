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

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianProcessor.h"
#include "LE3_2D_MASS_WL_UTILITIES/NoisyCatalogData.h"
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

#include "ST_DM_FilenameProvider/FilenameProvider.h"

using namespace LE3_2D_MASS_WL_UTILITIES;
using Logging = Elements::Logging;
using FileNameProvider = Euclid::DataModelTools::FileNameProvider;

namespace LE3_2D_MASS_WL_CARTESIAN {

static Logging logger = Logging::getLogger("CartesianProcessor");

CartesianProcessor::CartesianProcessor() : m_args(),
                m_shearCatalogPath(), m_shearMapPath(), m_clusterCatalog(),
                m_produceMcMaps(false)
{
}

CartesianProcessor::CartesianProcessor(
        std::map<std::string, variable_value>& args) : m_args(args),
                m_shearCatalogPath(), m_shearMapPath(), m_clusterCatalog(),
                m_produceMcMaps(false)
{
    parseOptions();
}

void CartesianProcessor::checkOption(const std::string & key)
{
    if(m_args.count(key) == 0)
    {
        logger.fatal() << "Option not found: " << key;
        throw ExitCode::USAGE;
    }
}

void CartesianProcessor::parseOptions()
{
    m_workdir = m_args["workdir"].as<std::string>();
    m_datadir = m_workdir / "data";
    m_paramfile = m_args["paramfile"].as<std::string>();
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
        // aka. Patch-Based E/B Convergence Maps Parameters
        checkOption("shear");
        m_shearCatalogPath = m_args["shear"].as<std::string>();
    }
    else if(m_paramtype == "DpdTwoDMassParamsConvergenceClusters")
    {
        // Case: small patches defined in cluster catalog, need a shear catalog
        // and a cluster catalog. Will process multiple patches in a single
        // redshift bins.
        // aka. Patch-Based E/B Convergence Maps Parameters for Clusters
        checkOption("shear");
        checkOption("clusterCatalog");
        m_shearCatalogPath = m_args["shear"].as<std::string>();
        m_clusterCatalog = m_args["clusterCatalog"].as<std::string>();
    }
    else if(m_paramtype == "DpdTwoDMassParamsConvergencePatchesToSphere")
    {
        // Case: multiple patches covering the sphere, need a shear catalog.
        // Will process multiple patches in multiple redshift bins (TBC, maybe
        // only 1 due to computation time).
        // aka. MultiPatch-Based E/B Convergence Maps Parameters
        checkOption("shear");
        m_shearCatalogPath = m_args["shear"].as<std::string>();
    }
    else if(m_paramtype == "DpdTwoDMassParamsPeakCatalogConvergence")
    {
        // Case: TODO
        // aka. Wavelet-based Peak Catalog Parameters
    }
    else if(m_paramtype == "DpdTwoDMassParamsPeakCatalogMassAperture")
    {
        // Case: TODO
        // aka. Peak Catalog Based on Aperture Mass Parameters
    }
    else if(m_paramtype == "DpdTwoDMassParamsConvergenceSphere")
    {
        // Case: TODO
        // aka. Spherical E/B Convergence Maps Parameters
    }
}

void CartesianProcessor::processPatch(CartesianAlgoKS& cartesianAlgoKS, CatalogData& cat)
{
    CartesianParam& params = cartesianAlgoKS.getCartesianParam();
    int ip = params.m_currIpatch;
    auto& patch = params.getPatch(ip);

    logger.info() << "Patch bounds:";
    logger.info() << "Dec: " << patch.getDecMin() << " "
                             << patch.getDecMax();
    logger.info() << "Ra: " << patch.getRaMin() << " "
                            << patch.getRaMax();
    logger.info() << "Width: " << patch.getPatchWidth();
    logger.info() << "Pixel size: " << patch.getPixelSize();

    // Set filenames for this patch and all redshift bins
    // TODO: fix this and use SWVersion
    FileNameProvider filename_provider("LE3", "2D-MASS-WL", "2.6");
    params.m_currOutFilesPath["ConvergenceNoisy"] =
            m_datadir / filename_provider.getFilename("ConvergenceNoisy_Patch" + std::to_string(ip), ".fits", true);
    params.m_currOutFilesPath["ConvergenceDenoised"] =
            m_datadir / filename_provider.getFilename("ConvergenceDenoised_Patch" + std::to_string(ip), ".fits", true);
    params.m_currOutFilesPath["SNR"] =
            m_datadir / filename_provider.getFilename("SNR_Patch" + std::to_string(ip), ".fits", true);

    // Iterate over redshift bins
    int Nbins = params.getNbins();
    logger.info() << "Start iteration over Nbins of redshift = " << Nbins;
    logger.info() << "----------------------------------------------------";
    for (int iz=0; iz<Nbins; iz++)
    {
        params.m_currIzbin = iz;
        processZbin(cartesianAlgoKS, cat);
    }

    // End of this patch, now, we have a fits saved to create the
    // product xml file

    // Create a DMOutput object
    DmOutput dm;
    const std::string product_type = "DpdTwoDMassConvergencePatch";
    auto product = initProduct<dpdTwoDMassConvergencePatch,
                               twoDMassCollectConvergencePatch,
                               int>(product_type, params.getNResamples());

    dm.createPatchXml(product, outputType::NoisedPatch, params.m_currOutFilesPath["ConvergenceNoisy"]);
    dm.createPatchXml(product, outputType::DenoisedPatch, params.m_currOutFilesPath["ConvergenceDenoised"]);
    dm.createPatchXml(product, outputType::SNRPatch, params.m_currOutFilesPath["SNR"]);

    auto outfile = m_workdir / fs::path(filename_provider.getFilename(
            "DpdTwoDMassConvergencePatch_" + std::to_string(ip), ".xml", true));
    writeProduct<dpdTwoDMassConvergencePatch>(product, outfile);
}

void CartesianProcessor::processZbin(CartesianAlgoKS& cartesianAlgoKS, CatalogData& cat)
{
    CartesianParam& params = cartesianAlgoKS.getCartesianParam();
    int ip = params.m_currIpatch;
    int iz = params.m_currIzbin;
    auto& patch = params.getPatch(ip);

    double zmin = params.getZMin(iz);
    double zmax = params.getZMax(iz);
    logger.info() << "Redshift bounds:";
    logger.info() << "zMin: " << zmin;
    logger.info() << "zMax: " << zmax;

    logger.info() << "Create shearMap before filling with data";
    logger.info() << "SizeXaxis: " << patch.getXbin();
    logger.info() << "SizeYaxis: " << patch.getYbin();

    // Extract patch into shear map
    ShearMap shearMap(patch.getXbin(), patch.getYbin(), 3);
    cartesianAlgoKS.extractShearMap(shearMap, cat, patch, zmin, zmax);

    if(shearMap.getNumberOfGalaxies() == 0)
    {
        logger.info() << "Nothing to do: continue...";
        logger.info() << "--------------------------------------------";
        return;
    }

    // Save shear map [test only, not a PF product]
    fs::path outShearMapFilename = getFitsFilenameFromBase("ShearMap_extracted_Patch"
            + std::to_string(ip) + "_Zbin" + std::to_string(iz));
    params.setExtName("SHEAR");
    shearMap.writeMap(m_workdir / outShearMapFilename, params);

    // Perform direct mass mapping (either simple or reduced) to get the
    // convergence map. First define files names, then create working
    // convergence map (kappa) without copying values
    ConvergenceMap convergenceMap(shearMap, false);
    cartesianAlgoKS.performReducedShear(shearMap, convergenceMap);

    // Noisy convergence and denoised convergence are output products for
    // twoDMassConvergencePatch and need to be saved in separate fits

    // Perform inverse mass mapping to get the denoised shear map
    cartesianAlgoKS.performInverseKSMassMapping(convergenceMap, shearMap);

    // Create SNR map to be filled inside loop over Nresamples
    // copy only density axis
    ConvergenceMap snrMap(convergenceMap, false);
    snrMap.singleAxisCopy(convergenceMap, 2);

    // Computes N monte-carlo maps by adding to the last convergence map
    // N shear map with randomized ellipticities
    ShearMap noisyShearMap(patch.getXbin(), patch.getYbin(), 3);
    int Nresamples = params.getNResamples();
    for(int is=0; is<Nresamples; is++)
    {
        // Create noisy shear map
        CatalogData randomizeCat = NoisyCatalogData::getNoisyCatalog(cat);
        cartesianAlgoKS.extractShearMap(noisyShearMap, randomizeCat,
                patch, zmin, zmax);

        // Add denoised shear
        ShearMap shearMapNoise(shearMap, false);
        cartesianAlgoKS.sumShear(shearMap, noisyShearMap, shearMapNoise);

        // Perform inversion with no denoising
        cartesianAlgoKS.performMassMapping(shearMapNoise, convergenceMap);

        // Accumulate in SNR map
        int denom = (is == Nresamples - 1) ? 1 : Nresamples;
        cartesianAlgoKS.accumulateSquareAndComputeSnr(convergenceMap,
                snrMap, denom, (is == Nresamples - 1));
    } // end loop over Nresamples

    // Save SNR map
    params.setExtName("SNR_ZBIN_" + std::to_string(iz));
    snrMap.writeMap(params.m_currOutFilesPath["SNR"], params);

    logger.info() << "------------------------------------------------";
}

void CartesianProcessor::process()
{
    CartesianParam params;
    CatalogData cat, clusterCat;

    // Load shear catalog (should always exist)
    cat.getCatalogData(m_workdir, m_shearCatalogPath);
    logger.info() << "Number of galaxies in catalog: " << cat.getNentries();

    // Load cluster catalog if exists
    if(not m_clusterCatalog.empty())
    {
        clusterCat.getCatalogData(m_workdir, m_clusterCatalog);
        logger.info() << "Number of clusters in catalog: " << clusterCat.getNentries();
    }

    // Read parameters (and patches definitions)
    readParameterFile(m_paramfile, params, cat, clusterCat);

    // Build cartesian algo
    CartesianAlgoKS cartesianAlgoKS(params);

    // Iterate over patches
    int Npatches = params.getNPatches();
    logger.info() << "Start iteration over NPatches = " << Npatches;
    for(int ip=0; ip<Npatches; ip++)
    {
        params.m_currIpatch = ip;
        processPatch(cartesianAlgoKS, cat);
    }
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN


