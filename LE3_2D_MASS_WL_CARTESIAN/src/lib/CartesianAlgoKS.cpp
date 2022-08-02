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
 * @file src/lib/CartesianAlgoKS.cpp
 * @date 10/21/19
 * @author user
 */

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "ElementsKernel/Temporary.h"
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

static Elements::Logging logger =
                                Elements::Logging::getLogger("CartesianAlgoKS");

using LE3_2D_MASS_WL_UTILITIES::DmInput;

namespace LE3_2D_MASS_WL_CARTESIAN
{

CartesianAlgoKS::CartesianAlgoKS(CartesianParam &params) :
        m_cartesianParam(params)
{
}

////////////////////////////////////////////////////////////////////////
// Extract ShearMap function
////////////////////////////////////////////////////////////////////////
void CartesianAlgoKS::extractShearMap(const fs::path& shearMapFilename,
        CatalogData& cat, const PatchDef& patch)
{
    extractShearMap(shearMapFilename.native(), cat, patch);
}

void CartesianAlgoKS::extractShearMap(ShearMap& shearMap,
        CatalogData& cat, const PatchDef& patch, double zmin, double zmax)
{
    MapMaker map(cat);
    map.getShearMap(patch, shearMap, zmin, zmax);

    // Pixelate X and Y axis
    if ((shearMap.getXdim()) == 2048 && (shearMap.getYdim()) == 2048)
    {
        shearMap.pixelate(1, 1);
    }
    if ((shearMap.getXdim()) == 4096 && (shearMap.getYdim()) == 4096)
    {
        shearMap.pixelate(2, 2);
    }
}

void CartesianAlgoKS::extractShearMap(const std::string& shearMapFilename,
        CatalogData& cat, const PatchDef& patch)
{
    logger.info() << "Create shearMap before filling with data";
    logger.info() << "SizeXaxis: " << patch.getXbin();
    logger.info() << "SizeYaxis: " << patch.getYbin();

    ShearMap shearMap(patch.getXbin(), patch.getYbin(), 3);

    extractShearMap(shearMap, cat, patch, 0, 10);

    // Writing Shear Map
    std::string name = "SHEAR_PATCH";
    m_cartesianParam.setExtName(name);
    shearMap.writeMap(shearMapFilename, m_cartesianParam);
}

////////////////////////////////////////////////////////////////////////
// Perform Inverse KS Mass Mapping
////////////////////////////////////////////////////////////////////////
void CartesianAlgoKS::performInverseKSMassMapping(const ConvergenceMap& convMap,
        ShearMap& shearMap)
{
    // Fill density image from conv map
    shearMap.singleAxisCopy(convMap, 2);

    convMap.getShearMap(shearMap);
}

////////////////////////////////////////////////////////////////////////
// Perform Reduced Shear Computation
////////////////////////////////////////////////////////////////////////
void CartesianAlgoKS::performReducedShear(fs::path& inShearMap,
        const fs::path& workdir, const fs::path& outputConvMapsJson)
{
    // get vector of input shear map filepath
    std::vector<fs::path> filenames = readFilenames(
            workdir / "data" / inShearMap);

    // iterate over input shear map(s)
    for (size_t i = 0; i < filenames.size(); i++)
    {
        logger.info() << "performReducedShear on shear map " << filenames[i];

        // load input shear map (gamma)
        ShearMap reducedShear(filenames[i]);

        if (m_cartesianParam.isAddBorder())
        {
            reducedShear.add_borders();
        }

        // create working convergence map (kappa) without copying values
        ConvergenceMap convergenceMap(reducedShear, false);
        ConvergenceMap tildeConvergenceMap;

        // fill density image from shear map
        convergenceMap.singleAxisCopy(reducedShear, 2);

        // create vector of output convergence maps filenames (same for noisy and denoised)
        std::vector<fs::path> convMapsFilenames;

        size_t pos = (filenames[i].string()).find("ShearMap");

        // loop over reduce iteration(s)
        for (size_t it = 0; it < m_cartesianParam.getRsNItReducedShear(); it++)
        {
            // write current reduced shear
            fs::path reducedShearMap(
                    "EUC_LE3_WL_ShearMap_" + std::to_string(it) + "_Reduced"
                            + (filenames[i].string()).substr(pos));
            reducedShear.writeMap(workdir / "data" / reducedShearMap,
                    m_cartesianParam);

            // perform either KS or KS+ mass mapping (no filtering performed)
            performMassMapping(reducedShear, convergenceMap, workdir);

            logger.info() << "come back after performMassMapping";

            // copy convergence to tilde convergence
            tildeConvergenceMap = ConvergenceMap(convergenceMap);

            logger.info() << "create copy into tildeConvergenceMap";

            // filter kappaE using RSSigmaGauss and apply threshold
            getTildeConvergence(tildeConvergenceMap);

            logger.info() << "come back after getTildeConvergence";

            // perform correction (after this: tildeConvergenceMap is not used anymore)
            reducedShear.correctReducedShear(tildeConvergenceMap);

            logger.info() << "come back after correctReducedShear";

        } // end of reduced shear iteration numbers

        logger.info() << "end of reduction";

        // in case borders were added, remove them
        if (m_cartesianParam.isAddBorder())
        {
            convergenceMap.remove_borders();
        }

        // set convergence filenames
        auto outNoisyConvMap = fs::path(
                "EUC_LE3_WL_NoisyConvergenceMapKS_"
                        + (filenames[i].string()).substr(pos));
        auto outDenoisedConvMap = fs::path(
                "EUC_LE3_WL_DenoisedConvergenceMapKS_"
                        + (filenames[i].string()).substr(pos));
        std::string name = "KAPPA_PATCH";
        m_cartesianParam.setExtName(name);

        // Update convergence map filename with KSPlus instead of KS
        if (m_cartesianParam.getNInpaint() != 0)
        {
            replaceSubStringInPath(outNoisyConvMap, "KS", "KSPlus");
            replaceSubStringInPath(outDenoisedConvMap, "KS", "KSPlus");
        }

        // write noisy convergence
        logger.info() << "Will write noisy convergence map to "
                << outNoisyConvMap;
        convergenceMap.writeMap(workdir / "data" / outNoisyConvMap,
                m_cartesianParam);
        convMapsFilenames.push_back(outNoisyConvMap);

        // if filter is activated, filter and write denoised convergence map
        if (fabs(m_cartesianParam.getGaussStd()) > 0.001)
        {
            convergenceMap.applyGaussianFilter(
                    m_cartesianParam.getGaussStd());
            convergenceMap.setKappaBToZero();

            logger.info() << "Will write denoised convergence map to "
                    << outDenoisedConvMap;
            convergenceMap.writeMap(workdir / "data" / outDenoisedConvMap,
                    m_cartesianParam);
            convMapsFilenames.push_back(outDenoisedConvMap);
        }

        // write json file with convergence maps
        fillJsonOutput(workdir, outputConvMapsJson, convMapsFilenames);

    } //end of number of maps iterations

    logger.info() << "end of performReducedShear";
}

////////////////////////////////////////////////////////////////////////
// Perform Mass Mapping Function
////////////////////////////////////////////////////////////////////////
void CartesianAlgoKS::performMassMapping(const ShearMap& shearMap,
        ConvergenceMap& convMap, const fs::path& workdir)
{
    logger.info("KS/KSPlus conversion from Shear to Convergence Map");

    // First mass mapping
    shearMap.getConvMap(convMap);

    // Perform inpainting
    if (m_cartesianParam.getNInpaint() != 0)
    {
        performInPainting(shearMap, convMap, workdir);
    }

    logger.info() << "Exit performMassMapping";
}

////////////////////////////////////////////////////////////////////////
// Get tilde convergence
////////////////////////////////////////////////////////////////////////
void CartesianAlgoKS::getTildeConvergence(ConvergenceMap& ConvMap)
{
    // apply gaussian filter on the map
    if (fabs(m_cartesianParam.getRsGaussStd()) > 0.001)
    {
        ConvMap.applyGaussianFilter(m_cartesianParam.getRsGaussStd());
    }

    // get kappaB std
    double std_dev = ConvMap.getSigma(1);
    logger.info() << "std dev: " << std_dev;

    // apply threshold
    ConvMap.applyThreshold(m_cartesianParam.getRsThreshold() * std_dev, 0);
}

////////////////////////////////////////////////////////////////////////
// Perform Inpainting
////////////////////////////////////////////////////////////////////////
void CartesianAlgoKS::performInPainting(const ShearMap& shearMap,
        ConvergenceMap& outConvMap, const fs::path& workdir)
{
    GenericMap mask(shearMap.getXdim(), shearMap.getYdim(), 1);
    InpaintingAlgo myIPalgo(shearMap, mask, m_cartesianParam);

    myIPalgo.performInPaintingAlgo(outConvMap);

    auto maskFilepath = fs::path(
            "InPaintingMask_" + getDateTimeString() + ".fits");
    mask.writeMap(workdir / "data" / maskFilepath, m_cartesianParam);

    logger.info() << "Exit performInPainting";
}

////////////////////////////////////////////////////////////////////////
void CartesianAlgoKS::getSNRMap(fs::path& workdir, fs::path& NReMaps,
        fs::path& outConvergenceMap)
{
    std::vector<fs::path> filenames = readFilenamesInJson(NReMaps);

    // extract NReSample maps
    std::vector<fs::path> filenamesNReMaps;
    for (size_t i = 0; i < filenames.size(); i++)
    {
        if ((filenames[i].string()).find("NReSample") != std::string::npos)
        {
            filenamesNReMaps.push_back(filenames[i]);
        }
    }

    size_t N = filenamesNReMaps.size();
    logger.info() << "number of sampled maps: " << N;

    // create ouput SNR convergenceMap
    ConvergenceMap snrMap;

    for (size_t i = 0; i < N; i++)
    {
        logger.info() << "filenames: " << filenamesNReMaps[i];

        // load convergence map
        ConvergenceMap SampledMap(filenamesNReMaps[i]);

        // resize SNR map to size of first one
        if (i == 0)
        {
            size_t xbin = SampledMap.getXdim();
            size_t ybin = SampledMap.getYdim();
            size_t zbin = SampledMap.getZdim();
            snrMap = ConvergenceMap(xbin, ybin, zbin);
            for (size_t z = 0; z < snrMap.getZdim(); z++)
            {
                snrMap.clear(z);
            }
        }

        // accumulate square of the values
        for (size_t z = 0; z < snrMap.getZdim(); z++)
        {
            for (size_t y = 0; y < snrMap.getYdim(); y++)
            {
                for (size_t x = 0; x < snrMap.getXdim(); x++)
                {
                    double val = snrMap.getBinValue(x, y, z)
                            + pow(SampledMap.getBinValue(x, y, z), 2);
                    snrMap.setBinValue(x, y, z, val);

                    // at last iteration, compute square root of mean
                    if (i == N - 1)
                    {
                        double mean = snrMap.getBinValue(x, y, z) / double(N);
                        snrMap.setBinValue(x, y, z, sqrt(mean));
                    }
                }
            }
        }
    }  // end loop on input convergence map

    // set default name if empty
    if (outConvergenceMap.string().empty())
    {
        outConvergenceMap = fs::path(
                "EUC_LE3_WL_SNRConvergenceMap_" + getDateTimeString()
                        + ".fits");
    }

    std::string name = "KAPPA_PATCH";
    m_cartesianParam.setExtName(name);
    snrMap.writeMap((workdir / "data" / outConvergenceMap).native(),
            m_cartesianParam);

}
////////////////////////////////////////////////////////////////////////
// write xml files
////////////////////////////////////////////////////////////////////////
void CartesianAlgoKS::writeXMLfile(fs::path& inputConvMapsJson,
        fs::path& outputConvMapsXml)
{
    // Create a DMOutput object
    DmOutput dm;

    std::string file = outputConvMapsXml.string();
    int NResamples = m_cartesianParam.getNResamples();

    auto filenameJson = inputConvMapsJson.filename();
    // auto workdir = outConvergenceMapJson.parent_path();

    std::vector<fs::path> filenames = readFilenamesInJson(filenameJson);
    const fs::path outfile { file };

    if ((m_cartesianParam.getParaFileType()).compare("Conv_Patch") == 0)
    {
        const std::string product_type = "DpdTwoDMassConvergencePatch";
        auto product = initProduct<dpdTwoDMassConvergencePatch,
                twoDMassCollectConvergencePatch, int>(product_type, NResamples);
        for (size_t i = 0; i < filenames.size(); i++)
        {
            fs::path mapOut(filenames[i].native());
            if ((mapOut.string()).find("NoisyConvergence") != std::string::npos)
            {
                logger.info() << "writing file: " << mapOut;
                dm.createNoisedPatchXml(product, mapOut);
            }
            if ((mapOut.string()).find("DenoisedConvergence")
                    != std::string::npos)
            {
                dm.createDenoisedPatchXml(product, mapOut);
            }
            if ((mapOut.string()).find("SNRConvergence") != std::string::npos)
            {
                dm.createSNRPatchOutputXml(product, mapOut);
            }
        }
        writeProduct<dpdTwoDMassConvergencePatch>(product, outfile);
    }

    logger.info() << "DM output products created in: " << file;
}

}// namespace LE3_2D_MASS_WL_CARTESIAN

