/**
 * @file src/lib/GenericParam.cpp
 * @date 02/04/22
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

#include "LE3_2D_MASS_WL_UTILITIES/GenericParam.h"

static Elements::Logging logger =  Elements::Logging::getLogger("GenericParam");

namespace LE3_2D_MASS_WL_UTILITIES {

/**
 * @brief   constructor (to create an object with default values)
 */
GenericParam::GenericParam()
{
    // parameters used to handle reduce shear
    m_RSCorrection = false;
    m_RSThreshold = 5.0;
    m_RSGaussSTD = 0.0;

    // denoising parameters for final map
    m_DenoisingAlgo = std::string("");
    m_GaussSTD = 0.0;
    m_ThresholdFDR = 0.0;

    // number of resampling that needs to be performed for Noise/SNR maps
    m_NResamples = 0;

    // how to handle redshift
    m_ZMin = {0.0};  // vector
    m_ZMax = {10.0}; // vector
    m_Nbins = 1;
    m_BalancedBins = false;
    m_ZMargin = 0.0;  // only for clusters
    m_ZMaxHalo = 1; // only for clusters

    // how to handle gaps in data
    m_NInpaint = 100;
    m_NInpScale = 0;
    m_EqualVarPerScale = false;
    m_ForceBMode = true;
    m_AddBorder = false; // only for cartesian

    // how to handle the projection
    m_Project = std::string("TAN");
    m_NPatches = 1;
    m_Patches = {};

    // how to handle the pixelization
    m_Nside = 0; // only for sphere
    m_overlap = 0; // only for patch2sphere

    // how to select from the cluster catalogue
    m_MassThreshold = 0.0;  // only for cluster

    // code defined parameters
    m_ExtName = std::string("");
    m_ParaFileType = std::string("");
}

void GenericParam::print() const
{
    std::cout << "RSCorrection: " << m_RSCorrection << std::endl;
    std::cout << "RSThreshold: " << m_RSThreshold << std::endl;
    std::cout << "RSGaussSTD: " << m_RSGaussSTD << std::endl;
    std::cout << "DenoisingAlgo: " << m_DenoisingAlgo << std::endl;
    std::cout << "GaussSTD: " << m_GaussSTD << std::endl;
    std::cout << "ThresholdFDR: " << m_ThresholdFDR << std::endl;
    std::cout << "NResamples: " << m_NResamples << std::endl;
    for(size_t i = 0; i < m_ZMin.size(); i++)
    {
        std::cout << "ZMin: " << m_ZMin[i] << std::endl;
        std::cout << "ZMax: " << m_ZMax[i] << std::endl;
    }
    std::cout << "Nbins: " << m_Nbins << std::endl;
    std::cout << "BalancedBins: " << m_BalancedBins << std::endl;
    std::cout << "ZMargin: " << m_ZMargin << std::endl;
    std::cout << "ZMaxHalo: " << m_ZMaxHalo << std::endl;
    std::cout << "NInpaint: " << m_NInpaint << std::endl;
    std::cout << "NInpScale: " << m_NInpScale << std::endl;
    std::cout << "EqualVarPerScale: " << m_EqualVarPerScale << std::endl;
    std::cout << "ForceBMode: " << m_ForceBMode << std::endl;
    std::cout << "AddBorder: " << m_AddBorder << std::endl;
    std::cout << "Project: " << m_Project << std::endl;
    std::cout << "NPatches: " << m_NPatches << std::endl;
    std::cout << "Nside: " << m_Nside << std::endl;
    std::cout << "overlap: " << m_overlap << std::endl;
    std::cout << "MassThreshold: " << m_MassThreshold << std::endl;
    std::cout << "ExtName: " << m_ExtName << std::endl;
    std::cout << "ParaFileType: " << m_ParaFileType << std::endl;
}

const std::vector<PatchDef>& GenericParam::getPatches() const
{
    return m_Patches;
}

const PatchDef& GenericParam::getPatch(int i) const
{
    return m_Patches[i];
}

bool GenericParam::isAddBorder() const
{
    return m_AddBorder;
}

void GenericParam::setAddBorder(bool addBorder)
{
    m_AddBorder = addBorder;
}

bool GenericParam::isBalancedBins() const
{
    return m_BalancedBins;
}

void GenericParam::setBalancedBins(bool balancedBins)
{
    m_BalancedBins = balancedBins;
}

const std::string& GenericParam::getDenoisingAlgo() const
{
    return m_DenoisingAlgo;
}

void GenericParam::setDenoisingAlgo(const std::string& denoisingAlgo)
{
    m_DenoisingAlgo = denoisingAlgo;
}

bool GenericParam::isEqualVarPerScale() const
{
    return m_EqualVarPerScale;
}

void GenericParam::setEqualVarPerScale(bool equalVarPerScale)
{
    m_EqualVarPerScale = equalVarPerScale;
}

const std::string& GenericParam::getExtName() const
{
    return m_ExtName;
}

void GenericParam::setExtName(const std::string& extName)
{
    m_ExtName = extName;
}

bool GenericParam::isForceBMode() const
{
    return m_ForceBMode;
}

void GenericParam::setForceBMode(bool forceBMode)
{
    m_ForceBMode = forceBMode;
}

double GenericParam::getGaussStd() const
{
    return m_GaussSTD;
}

void GenericParam::setGaussStd(double gaussStd)
{
    m_GaussSTD = gaussStd;
}

double GenericParam::getMassThreshold() const
{
    return m_MassThreshold;
}

void GenericParam::setMassThreshold(double massThreshold)
{
    m_MassThreshold = massThreshold;
}

int GenericParam::getNbins() const
{
    return m_Nbins;
}

void GenericParam::setNbins(int nbins)
{
    m_Nbins = nbins;
}

int GenericParam::getNInpaint() const
{
    return m_NInpaint;
}

void GenericParam::setNInpaint(int nInpaint)
{
    m_NInpaint = nInpaint;
}

int GenericParam::getNInpScale() const
{
    return m_NInpScale;
}

void GenericParam::setNInpScale(int nInpScale)
{
    m_NInpScale = nInpScale;
}

unsigned int GenericParam::getNPatches() const
{
    return m_NPatches;
}

void GenericParam::setNPatches(int nPatches)
{
    m_NPatches = nPatches;
}

int GenericParam::getNResamples() const
{
    return m_NResamples;
}

void GenericParam::setNResamples(int nResamples)
{
    m_NResamples = nResamples;
}

int GenericParam::getNside() const
{
    return m_Nside;
}

void GenericParam::setNside(int nside)
{
    m_Nside = nside;
}

int GenericParam::getOverlap() const
{
    return m_overlap;
}

void GenericParam::setOverlap(int overlap)
{
    m_overlap = overlap;
}

const std::string& GenericParam::getParaFileType() const
{
    return m_ParaFileType;
}

void GenericParam::setParaFileType(const std::string& paraFileType)
{
    m_ParaFileType = paraFileType;
}

const std::string& GenericParam::getProject() const
{
    return m_Project;
}

void GenericParam::setProject(const std::string& project)
{
    m_Project = project;
}

double GenericParam::getRsGaussStd() const
{
    return m_RSGaussSTD;
}

void GenericParam::setRsGaussStd(double rsGaussStd)
{
    m_RSGaussSTD = rsGaussStd;
}

bool GenericParam::getRsCorrection() const
{
    return m_RSCorrection;
}

void GenericParam::setRsCorrection(bool rsCorrection)
{
    m_RSCorrection = rsCorrection;
}

double GenericParam::getRsThreshold() const
{
    return m_RSThreshold;
}

void GenericParam::setRsThreshold(double rsThreshold)
{
    m_RSThreshold = rsThreshold;
}

double GenericParam::getThresholdFdr() const
{
    return m_ThresholdFDR;
}

void GenericParam::setThresholdFdr(double thresholdFdr)
{
    m_ThresholdFDR = thresholdFdr;
}

double GenericParam::getZMargin() const
{
    return m_ZMargin;
}

void GenericParam::setZMargin(double zMargin)
{
    m_ZMargin = zMargin;
}

double GenericParam::getZMax(int i) const
{
    return m_ZMax[i];
}


double GenericParam::getZMaxHalo() const
{
    return m_ZMaxHalo;
}

void GenericParam::setZMaxHalo(double zMaxHalo)
{
    m_ZMaxHalo = zMaxHalo;
}

double GenericParam::getZMin(int i) const
{
    return m_ZMin[i];
}

}  // namespace LE3_2D_MASS_WL_UTILITIES



