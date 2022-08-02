/**
 * @file src/lib/GetCartesianMCMaps.cpp
 * @date 10/13/20
 * @author vkansal
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

#include "LE3_2D_MASS_WL_CARTESIAN/GetCartesianMCMaps.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MapMaker.h"
#include "LE3_2D_MASS_WL_UTILITIES/NoisyCatalogData.h"
#include "LE3_2D_MASS_WL_UTILITIES/Matrix.h"
#include "ElementsKernel/Temporary.h"
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

static Elements::Logging logger = Elements::Logging::getLogger(
        "CartesianMCMaps");

using namespace LE3_2D_MASS_WL_UTILITIES;

namespace LE3_2D_MASS_WL_CARTESIAN
{

GetCartesianMCMaps::GetCartesianMCMaps(CatalogData& inputData,
        CartesianParam &cartesianParam) :
        m_inData(inputData), m_cartesianParam(cartesianParam),
        m_cartesianAlgo(cartesianParam)
{
}

void GetCartesianMCMaps::getDeNoisedShearMap(ShearMap& outDenoisedShearMap)
{
    Elements::TempDir one;

    fs::path temp_path = one.path();
    // fs::path temp_path = "/home/user/Work/workdir/out_MC";

    auto datadir = temp_path / "data";
    fs::create_directories(datadir);

    // Save this map back in a new file
    fs::path shearMapFilename = "ShearMapFits.fits";
    fs::path shearMapPath = datadir / shearMapFilename;

    fs::path convMapsFilePath = "convMapFits.json";

    for (size_t it = 0; it < m_cartesianParam.getNPatches(); it++)
    {
        auto& patch = m_cartesianParam.getPatches()[it];
        ShearMap deNoisedShearMap;

        // perform extraction and get (reduced) noised shear map
        m_cartesianAlgo.extractShearMap(shearMapPath.native(), m_inData, patch);

        //perform mass mapping (either KS or KS+) on one shear map in fits
        m_cartesianAlgo.performReducedShear(shearMapFilename, temp_path,
            convMapsFilePath);

        // load final denoised convergence (second entry of vector)
        // filter must be activated
        std::vector<fs::path> convMapsJson = readFilenamesInJson(temp_path / convMapsFilePath);
        ConvergenceMap convMap(temp_path / "data" / convMapsJson.back());

        // perform inversion to get denoised shear map (gamma1, gamma2)
        outDenoisedShearMap = ShearMap(convMap, false);
        outDenoisedShearMap.singleAxisCopy(convMap, 2);
        convMap.getShearMap(outDenoisedShearMap);
    }
}

void GetCartesianMCMaps::getNoisedShearMaps(
        std::vector<ShearMap>& outShearMapList)
{
    int Niter = m_cartesianParam.getNResamples();
    if (Niter <= 0)
    {
        Niter = 1;
    }

    outShearMapList.resize(Niter);
    NoisyCatalogData randomise;
    auto& patch = m_cartesianParam.getPatches()[0];
    for (int i = 0; i != Niter; ++i)
    {
        outShearMapList[i].updateSizes(patch.getXbin(),patch.getYbin(),3);
        CatalogData randomizeCat = randomise.getNoisyCatalog(m_inData);
        MapMaker mapMaker(randomizeCat);
        // TODO: fix redshift here
        mapMaker.getShearMap(patch, outShearMapList[i], 0, 10);
    }
}

void GetCartesianMCMaps::performAddition(ShearMap& DenoisedShearMap,
        ShearMap& NoisedShearMap, const std::string& shearMapFilename)
{
    ShearMap shearMap(DenoisedShearMap, false);
    auto& patch = m_cartesianParam.getPatches()[0];
    for (size_t j = 0; j < patch.getXbin(); j++)
    {
        for (size_t i = 0; i < patch.getYbin(); i++)
        {
            std::complex<double> gamma(DenoisedShearMap.getBinValue(i, j, 0),
                    DenoisedShearMap.getBinValue(i, j, 1));
            std::complex<double> e_obs(NoisedShearMap.getBinValue(i, j, 0),
                    NoisedShearMap.getBinValue(i, j, 1));

            std::complex<double> eps_gamma = (e_obs + gamma)
                    / (std::complex<double>(1, 0) + std::conj(gamma) * e_obs);
            double gamma1_noise = std::real(eps_gamma);
            double gamma2_noise = std::imag(eps_gamma);

            shearMap.setBinValue(i, j, 0, gamma1_noise);
            if (m_cartesianParam.isForceBMode() == true)
            {
                shearMap.setBinValue(i, j, 1, 0);
            }
            else
            {
                shearMap.setBinValue(i, j, 1, gamma2_noise);
            }
        }
    }

    shearMap.writeMap(shearMapFilename, m_cartesianParam);
}

}  // namespace LE3_2D_MASS_WL_CARTESIAN
