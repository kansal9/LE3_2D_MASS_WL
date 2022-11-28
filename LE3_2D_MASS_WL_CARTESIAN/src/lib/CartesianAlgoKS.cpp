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
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

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

CartesianParam& CartesianAlgoKS::getCartesianParam() const
{
    return m_cartesianParam;
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
// Perform Mass Mapping Function
////////////////////////////////////////////////////////////////////////
void CartesianAlgoKS::performMassMapping(const ShearMap& shearMap,
                                         ConvergenceMap& convMap )
{
    logger.info("........ KS/KSPlus conversion from Shear to Convergence Map");

    // First mass mapping
    shearMap.getConvMap(convMap);

    // Perform inpainting
    if (m_cartesianParam.getNInpaint() != 0)
    {
        performInPainting(shearMap, convMap);
        logger.info() << ".......... exit performInPainting";
    }

    logger.info() << "........ exit performMassMapping";
}

////////////////////////////////////////////////////////////////////////
// Perform Reduced Shear Computation
////////////////////////////////////////////////////////////////////////
void CartesianAlgoKS::performReducedShear(const ShearMap& inShearMap,
                                          ConvergenceMap& outConvergenceMap)
{
    // copy the input of shear map to working one, so it is not modified
    ShearMap shearMap = ShearMap(inShearMap);

    if (m_cartesianParam.isAddBorder())
    {
        shearMap.add_borders();
    }

    // fill density image from shear map, other axis are filled with 0
    outConvergenceMap.singleAxisCopy(inShearMap, 2);

    for (size_t it = 0; it <= REDUCESHEARNITER; it++)
    {
        // do not perform reduction at first iteration because convergence is
        // initialised at 0
        if(it > 0)
        {
            // filter kappaE using RSSigmaGauss and apply threshold to get an estimate
            // of the E-mode to be used to correct the reduced shear
            getTildeConvergence(outConvergenceMap);
            logger.info() << "...... come back after getTildeConvergence";

            // perform correction (after this: shearMap is a better estimate of the true shear
            shearMap.correctReducedShear(outConvergenceMap);
            logger.info() << "...... come back after correctReducedShear";
        }

        // perform either KS or KS+ mass mapping (no filtering performed) to
        // get new convergence
        performMassMapping(shearMap, outConvergenceMap);
        logger.info() << ".... come back after performMassMapping";

        // if correction is not activated, break the loop here, else do 3 iterations
        if(!m_cartesianParam.getRsCorrection())
        {
            break;
        }
    } // end of reduced shear iteration
    logger.info() << ".. end of reduction";

    // in case borders were added, remove them
    if (m_cartesianParam.isAddBorder())
    {
        outConvergenceMap.remove_borders();
    }

    // save the noisy convergence map before denoising (backup original gaussStd)
    if(m_cartesianParam.getParaFileType() == parameterType::DpdTwoDMassParamsConvergencePatch)
    {
        double gaussStd = m_cartesianParam.getGaussStd();
        m_cartesianParam.setGaussStd(0);
        logger.info() << "Saving noisy convergence map";
        m_cartesianParam.setExtName("KAPPA_ZBIN_" + std::to_string(m_cartesianParam.m_currIzbin));
        outConvergenceMap.writeMap(m_cartesianParam.m_currOutFilesPath["ConvergenceNoisy"], m_cartesianParam);
        // restore gaussStd value
        m_cartesianParam.setGaussStd(gaussStd);
    }

    // if filter is activated, filter and write denoised convergence map
    if (fabs(m_cartesianParam.getGaussStd()) > 0.001)
    {
        outConvergenceMap.applyGaussianFilter(m_cartesianParam.getGaussStd());
        outConvergenceMap.setKappaBToZero();

        // save the denoised convergence map
        if(m_cartesianParam.getParaFileType() == parameterType::DpdTwoDMassParamsConvergencePatch)
        {
            logger.info() << "saving denoised convergence map";
            m_cartesianParam.setExtName("KAPPA_ZBIN_" + std::to_string(m_cartesianParam.m_currIzbin));
            outConvergenceMap.writeMap(m_cartesianParam.m_currOutFilesPath["ConvergenceDenoised"], m_cartesianParam);
        }
    }

    logger.info() << "end of performReducedShear";
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
    logger.debug() << "std dev: " << std_dev;

    // apply threshold
    ConvMap.applyThreshold(m_cartesianParam.getRsThreshold() * std_dev, 0);
}

////////////////////////////////////////////////////////////////////////
// Perform Inpainting
////////////////////////////////////////////////////////////////////////
void CartesianAlgoKS::performInPainting(const ShearMap& shearMap,
                                        ConvergenceMap& outConvMap)
{
    GenericMap mask(shearMap.getXdim(), shearMap.getYdim(), 1);

    // In all Inpainting methods, shearMap is const ref to the reduced shear
    InpaintingAlgo myIPalgo(shearMap, mask, m_cartesianParam);

    myIPalgo.performInPaintingAlgo(outConvMap);
}

////////////////////////////////////////////////////////////////////////
// Sum two shear maps
////////////////////////////////////////////////////////////////////////
void CartesianAlgoKS::sumShear(const ShearMap& shearMapA, const ShearMap& shearMapB,
        ShearMap& shearMapSum)
{
    assert (shearMapA.getXdim() == shearMapB.getXdim());
    assert (shearMapA.getXdim() == shearMapSum.getXdim());
    assert (shearMapA.getYdim() == shearMapB.getYdim());
    assert (shearMapA.getYdim() == shearMapSum.getYdim());

    for (size_t j = 0; j < shearMapA.getXdim(); j++)
    {
        for (size_t i = 0; i < shearMapA.getYdim(); i++)
        {
            std::complex<double> gA(shearMapA.getBinValue(i, j, 0),
                                    shearMapA.getBinValue(i, j, 1));

            std::complex<double> gB(shearMapB.getBinValue(i, j, 0),
                                    shearMapB.getBinValue(i, j, 1));

            std::complex<double> gSum = (gA + gB) / (std::complex<double>(1, 0) + std::conj(gA) * gB);
            double real = std::real(gSum);
            double imag = (m_cartesianParam.isForceBMode() == true) ? 0 : std::imag(gSum);

            shearMapSum.setBinValue(i, j, 0, real);
            shearMapSum.setBinValue(i, j, 1, imag);
        }
    }
}

////////////////////////////////////////////////////////////////////////
// Accumulate squares
////////////////////////////////////////////////////////////////////////
void CartesianAlgoKS::accumulateSquareAndComputeSnr(const ConvergenceMap& convMap,
                                                    ConvergenceMap& snrMap,
                                                    int denom = 1,
                                                    bool computeSnr = false)
{
    // accumulate square of the values
    for (size_t z = 0; z < snrMap.getZdim(); z++)
    {
        for (size_t y = 0; y < snrMap.getYdim(); y++)
        {
            for (size_t x = 0; x < snrMap.getXdim(); x++)
            {
                double val = snrMap.getBinValue(x, y, z)
                        + pow(convMap.getBinValue(x, y, z), 2);
                snrMap.setBinValue(x, y, z, val);

                if (computeSnr)
                {
                    snrMap.setBinValue(x, y, z, sqrt(val / double(denom)));
                }
            }
        }
    }
}

}// namespace LE3_2D_MASS_WL_CARTESIAN

