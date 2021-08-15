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

#include "LE3_2D_MASS_WL_CARTESIAN/GetMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ConvergenceMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/ShearMap.h"
#include "LE3_2D_MASS_WL_CARTESIAN/InpaintingAlgo.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MapMaker.h"
#include "LE3_2D_MASS_WL_CARTESIAN/MassMapping.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h"
#include "LE3_2D_MASS_WL_CARTESIAN/CoordinateBound.h"

#include <string>
#include "boost/multi_array.hpp"
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <experimental/filesystem>
#include <algorithm>
#include <utility>
// Datamodel for OUTPUT products
#include "LE3_2D_MASS_WL_UTILITIES/DmOutput.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-ConvergencePatch.h"

namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {
   namespace CartesianKS {
/**
 * @class CartesianAlgoKS
 * @brief
 *
 */
class CartesianAlgoKS {

public:
    /**
     * @brief Constructor
     */
    CartesianAlgoKS(LE3_2D_MASS_WL_CARTESIAN::CartesianParam &params);
    /**
     * @brief Destructor
     */
    virtual ~CartesianAlgoKS() = default;
   /**
    * @brief     It extracts the shear Map from catalog
    * @param     <shearMap>, <string> name of the output shearMap
    * @param     <Data>, <std::vector<std::vector<double> >> Catalog data
    * @return    <bool> true if shear map well extracted/created
   */
    bool extractShearMap(const std::string& shearMap, std::vector<std::vector<double> >& Data,
                         LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& CB);
   /**
    * @brief     It extracts the convergence Map from catalog
    * @param     <convergenceMap>, <string> name of the output convergenceMap
    * @param     <Data>, <std::vector<std::vector<double> >> Catalog data
    * @return    <bool> true if map is well extracted/created
   */
    bool extractConvergnceMap (const std::string& convergenceMap, std::vector<std::vector<double> >& Data,
                               LE3_2D_MASS_WL_CARTESIAN::CoordinateBound& CB);
   /**
    * @brief     This function performs KS mass mapping
    * @param     <shearMap>, <string> name of the shearMap
    * @param     <outConvMap>, <string> name of the output convergence Map
    * @return    <bool> true if convergence map well created
   */
    bool performKSMassMapping(const std::string& shearMap, const std::string& outConvMap);
   /**
    * @brief     This function performs inverse KS mass mapping
    * @param     <outConvMap>, <string> name of the input convergence Map
    * @param     <shearMap>, <string> name of the output shearMap
    * @return    <bool> true if shear map well created
   */
    bool performInverseKSMassMapping(const std::string& convMap, const std::string& outShearMap);
   /**
    * @brief     This function performs reduced shear computation
    * @param     <convMap>, <string> name of the convergence Map file with respected path
    * @param     <shearMap>, <string> name of the shearMap file with respected path
    * @param     <reducedShearMap>, <string> name of the output reduced shearMap file with respected path
    * @return    <bool> true if Reduced shear map well created
   */
 //   bool performReducedShear(const std::string& convMap, const std::string& shearMap,
  //        const std::string& reducedShearMap);
   /**
    * @brief     This function performs Inpainting
    * @param     <convMap>, <string> name of the convergence Map file with respected path
    * @param     <shearMap>, <string> name of the shearMap file with respected path
    * @param     <outConvMap>, <string> name of the output Convergence Map file with respected path
    * @return    <bool> true if convergence map well created
   */
    bool performInPainting(const std::string& convMap, const std::string& shearMap, const std::string& outConvMap);

   /**
    * @brief   This function performs mass mapping main function
    * @param   <InShearMap>, <boost::filesystem::path> name of the shearMap
    * @param   <outConvergenceMap>, <boost::filesystem::path> name of the output convergence Map
    * @param   <workdir>, <boost::filesystem::path> work directory
    * @return  <bool> true if convergsence map well created
   */
    bool perform_MassMapping_Function(boost::filesystem::path& InShearMap, boost::filesystem::path& outConvergenceMap,
          boost::filesystem::path& workdir);
   /**
    * @brief     This function performs reduced shear function
    * @param     <InShearMap>, <boost::filesystem::path> name of the shearMap
    * @param     <outConvergenceMap>, <boost::filesystem::path> name of the output convergence Map
    * @param     <workdir>, <boost::filesystem::path> work directory
    * @param     <outputConvMaps>, <boost::filesystem::path> name of the json file which contains convergence maps name
    * @return    <bool> true if convergence map well created
   */
    bool performReducedShear(boost::filesystem::path& InShearMap, boost::filesystem::path& outConvergenceMap,
                boost::filesystem::path& workdir, boost::filesystem::path outputConvMaps);
   //bool perform_ReducedShear_Computation(boost::filesystem::path& InShearMap,
   //                  boost::filesystem::path& outConvergenceMap, boost::filesystem::path& workdir,
   //                                                   boost::filesystem::path outXMLConvergenceMap);

   void getReducedShearMap(LE3_2D_MASS_WL_CARTESIAN::ShearMap& reducedShear,
                LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap& ConvMap);
   void getTildeConvergence(LE3_2D_MASS_WL_CARTESIAN::ConvergenceMap& ConvMap);

   /**
    * @brief     This function sets name of output file for noisy convergence map and perform mass mapping
    * @param     <workdir>, <boost::filesystem::path> work directory
    * @param     <InShearMap>, <boost::filesystem::path> name of the shearMap
    * @param     <outConvergenceMap>, <boost::filesystem::path> name of the output convergence Map
    * @return    <bool> true if convergence map well created
   */
    bool getNoisyConvergenceMap(boost::filesystem::path& workdir, boost::filesystem::path& InShearMap,
                                 boost::filesystem::path& outConvergenceMap);

   /**
    * @brief     This function sets name of output file for Denoised convergence map and perform mass mapping
    * @param     <workdir>, <boost::filesystem::path> work directory
    * @param     <InShearMap>, <boost::filesystem::path> name of the shearMap
    * @param     <outConvergenceMap>, <boost::filesystem::path> name of the output convergence Map
    * @return    <bool> true if convergence map well created
   */
    bool getDenoisedConvergenceMap(boost::filesystem::path& workdir, boost::filesystem::path& InShearMap,
                                   boost::filesystem::path& outConvergenceMap);

   /**
    * @brief     This function returns sampled convergence maps in json file
    * @param     <workdir>, <boost::filesystem::path> work directory
    * @param     <filenames>, vector<boost::filesystem::path> vector of names of the shearMap
    * @return    <boost::filesystem::path> returns json file containing names of the convergence maps
   */
    boost::filesystem::path getResampledMaps(boost::filesystem::path& workdir, std::vector<fs::path>& filenames);

   /**
    * @brief     This function provides snr map
    * @param     <workdir>, <boost::filesystem::path> work directory
    * @param     <NReMaps>, <boost::filesystem::path> names of the sampled convergence Maps
    * @param     <outConvergenceMap>, <boost::filesystem::path> name of the output SNR convergence Map
    * @return    <bool> returns true if snr convergence map well created
   */
    bool getSNRMap(boost::filesystem::path& workdir, boost::filesystem::path& NReMaps,
                   boost::filesystem::path& outConvergenceMap);

   /**
    * @brief     This function writess xml files
    * @param     <outConvergenceMapJson>, <boost::filesystem::path> name of the output convergence Maps in json file
    * @param     <outXMLConvergenceMap>, <boost::filesystem::path> name of the output XML file
    * @return    <bool> true if convergence map well created
   */
    bool writeXMLfile (boost::filesystem::path& outConvergenceMapJson, boost::filesystem::path& outXMLConvergenceMap);

private:
    /**
     *  @brief <m_cartesianParam>, CartesianParam object with catalog parameters
    */
    LE3_2D_MASS_WL_CARTESIAN::CartesianParam m_cartesianParam;

}; /* End of CartesianAlgoKS class */
   } /* namespace CartesianKS */
  } /* namespace TwoDMass */
 } /* namespace WeakLensing */
} /* namespace Euclid */

#endif
