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
 * @file LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h
 * @date 10/21/19
 * @author user
 */

#ifndef _CARTESIANALGOKS_H
#define _CARTESIANALGOKS_H

#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/InpaintingAlgo.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MapMaker.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/GenericMap.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmOutput.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"
#include "LE3_2D_MASS_WL_UTILITIES/PatchDef.h"
#include <string>
#include "boost/multi_array.hpp"
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <experimental/filesystem>
#include <algorithm>
#include <utility>

#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-ConvergencePatch.h"

using namespace dpd::le3::wl::twodmass::out::convergencepatch;
using namespace pro::le3::wl::twodmass;
using namespace LE3_2D_MASS_WL_UTILITIES;

namespace LE3_2D_MASS_WL_CARTESIAN
{
/**
 * @class CartesianAlgoKS
 * @brief
 *
 */
class CartesianAlgoKS
{

public:
    /**
     * @brief Constructor
     */
    CartesianAlgoKS(CartesianParam &params);
    /**
     * @brief Destructor
     */
    virtual ~CartesianAlgoKS() = default;

    /**
     * @brief     It extracts the shear Map from catalog
     * @param     <shearMap>, <path> name of the output shearMap
     * @param     <Data>, <std::vector<std::vector<double> >> Catalog data
     */
    void extractShearMap(const fs::path& shearMap, CatalogData& cat,
            const PatchDef& patch);

    /**
     * @brief     It extracts the shear Map from catalog
     * @param     <shearMap>, shear map to be filled
     * @param     <Data>, <std::vector<std::vector<double> >> Catalog data
     */
    void extractShearMap(ShearMap& shearMap, CatalogData& cat,
            const PatchDef& patch, double zmin, double zmax);

    /**
     * @brief     It extracts the shear Map from catalog
     * @param     <shearMap>, <string> name of the output shearMap
     * @param     <Data>, <std::vector<std::vector<double> >> Catalog data
     */
    void extractShearMap(const std::string& shearMap, CatalogData& cat,
            const PatchDef& patch);

    /**
     * @brief     This function performs inverse KS mass mapping
     * @param     <convMap>, input convergence Map
     * @param     <shearMap>, output shearMap
     */
    void performInverseKSMassMapping(const ConvergenceMap& convMap,
            ShearMap& outShearMap);

    /**
     * @brief     This function performs Inpainting
     * @param     <shearMap>, input shear map
     * @param     <outConvMap>, output convergence map
     * @param   <workdir>, <boost::filesystem::path> work directory
     */
    void performInPainting(const ShearMap& shearMap, ConvergenceMap& outConvMap);

    /**
     * @brief     This function sums to shear
     * @param     <shearMapA>, input shear map A
     * @param     <shearMapB>, input shear map B
     * @param     <shearMapSum>, output shear map
     */
    void sumShear(const ShearMap& shearMapA, const ShearMap& shearMapB,
                  ShearMap& shearMapSum);

    /**
     * @brief   This function performs mass mapping main function
     * @param   <shearMap>, shear map
     * @param   <convMap>, output convergence Map
     * @param   <workdir>, <boost::filesystem::path> work directory
     */
    void performMassMapping(const ShearMap& shearMap, ConvergenceMap& convMap);

    /**
     * @brief     This function performs reduced shear function
     * @param     <inShearMap>, <boost::filesystem::path> name of the shearMap
     * @param     <workdir>, <boost::filesystem::path> work directory
     * @param     <outputConvMaps>, <boost::filesystem::path> name of the json/
     *            xml file which contains convergence maps name
     */
    void performReducedShear(const ShearMap& inShearMap,
                             ConvergenceMap& outConvergenceMap);

    /**
     * @brief     This function turns input convergence into tilde convergence
     *            using filtering and thresholding according sigma kappaB
     * @param     <ConvMap>, input/output convergence map
     */
    void getTildeConvergence(ConvergenceMap& ConvMap);

    /**
     * @brief     This function accumulates square of convMap into snrMap
     *            and compute snr when computeSnr is true
     * @param     <ConvMap>, input convergence map
     * @param     <snrMap>, output snr map
     * @param     <denom>, denominator used is snr computation
     * @param     <computeSnr>, do compute snr
     */
    void accumulateSquareAndComputeSnr(const ConvergenceMap& convMap,
                                       ConvergenceMap& snrMap,
                                       int denom,
                                       bool computeSnr);


    /**
     * @brief     This function provides snr map
     * @param     <workdir>, <boost::filesystem::path> work directory
     * @param     <NReMaps>, <boost::filesystem::path> names of the sampled convergence Maps
     * @param     <outConvergenceMap>, <boost::filesystem::path> name of the output SNR convergence Map
     */
    void getSNRMap(fs::path& workdir, fs::path& NReMaps,
            fs::path& outConvergenceMap);

    /**
     * @brief     This function writess xml files
     * @param     <outConvergenceMapJson>, <boost::filesystem::path> name of the output convergence Maps in json file
     * @param     <outXMLConvergenceMap>, <boost::filesystem::path> name of the output XML file
     */
    void writeXMLfile(fs::path& inputConvMapsJson, fs::path& outputConvMapsXml);

    /**
     * @brief     This function returns a reference to the internal CartesianParam
     */
    CartesianParam& getCartesianParam() const;

private:
    /**
     *  @brief <m_cartesianParam>, CartesianParam object with catalog parameters
     */
    CartesianParam& m_cartesianParam;

}; /* End of CartesianAlgoKS class */

}// namespace LE3_2D_MASS_WL_CARTESIAN

#endif
