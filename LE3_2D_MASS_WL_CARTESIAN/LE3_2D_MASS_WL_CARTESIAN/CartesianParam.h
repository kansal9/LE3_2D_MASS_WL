/**
 * @file LE3_2D_MASS_WL_CARTESIAN/CartesianParam.h
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

#ifndef _LE3_2D_MASS_WL_CARTESIAN_CARTESIANPARAM_H
#define _LE3_2D_MASS_WL_CARTESIAN_CARTESIANPARAM_H

#include <string>
#include <boost/variant.hpp>
#include <boost/filesystem.hpp>
#include "ElementsKernel/Logging.h"

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_UTILITIES/GenericParam.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"

#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergencePatch.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergenceClusters.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergencePatchesToSphere.h"
#include "ST_DataModelBindings/dictionary/pro/le3/wl/twodmass/euc-test-le3-wl-twodmass.h"

using LE3_2D_MASS_WL_UTILITIES::GenericParam;
using LE3_2D_MASS_WL_UTILITIES::CatalogData;

namespace LE3_2D_MASS_WL_CARTESIAN
{

/**
 * @class CartesianParam
 * @brief
 */
class CartesianParam : public GenericParam
{

public:
    /**
     * @brief   Default constructor
     */
    CartesianParam();

    /**
     * @brief   function to read the convergence patches parameter file in XML wrt dpd
     * @param   <paramConvergencePatch> parameter filename in <string> format
     * @param   <cat> shear catalog to define balances redshift
     */
    void readConvPatchXMLFile(const fs::path& paramConvergencePatch,
                              const CatalogData& cat);

    /**
     * @brief   function to read the convergence Clusters parameter file in XML wrt dpd
     * @param   <paramConvClusters> parameter filename in <string> format
     * @param   <clusterCat> cluster catalog to define patches
     */
    void readConvClustersXMLFile(const fs::path& paramConvClusters,
                                 const CatalogData& clusterCat = CatalogData());

    /**
     * @brief   function to read the parameter file to create convergence Patches to sphere in XML wrt dpd
     * @param   <paramConvClusters> parameter filename in <string> format
     */
    void readConvPatchesToSphereXMLFile(
            const fs::path& paramConvPatchesToSphere, const CatalogData& cat);

    /**
     * current redshift bin index
     */
    int m_currIzbin = 0;

    /**
     * current patch index
     */
    int m_currIpatch = 0;

    /**
     * map of product filenames
     */
    std::map<std::string, fs::path> m_currOutFilesPath;


private:

}; // End of CartesianParam class

/**
 * @brief     reads the parameter file and fill the param object
 * @param     <ParamFile>, <boost::filesystem::path> Parameter filename with path
 * @param     <params>, <LE3_2D_MASS_WL_CARTESIAN::CartesianParam> object to return parameters
 */
void readParameterFile(const fs::path& ParamFile,
        CartesianParam &params, const CatalogData& shearCat,
        const CatalogData& clusterCat = CatalogData());

}  // namespace LE3_2D_MASS_WL_CARTESIAN

#endif
