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
 * @file LE3_2D_MASS_WL_SPHERICAL/SphericalParam.h
 * @date 07/14/19
 * @author user
 */

#ifndef _LE3_2D_MASS_WL_SPHERICAL_SPHERICALPARAM_H
#define _LE3_2D_MASS_WL_SPHERICAL_SPHERICALPARAM_H

#include <string>
#include <boost/filesystem.hpp>
#include "ElementsKernel/Logging.h"

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"
#include "LE3_2D_MASS_WL_UTILITIES/GenericParam.h"
#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"

#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergenceSphere.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsConvergencePatchesToSphere.h"

using LE3_2D_MASS_WL_UTILITIES::GenericParam;
using LE3_2D_MASS_WL_UTILITIES::CatalogData;

namespace LE3_2D_MASS_WL_SPHERICAL
{

/**
 * @class SphericalParam
 * @brief
 */
class SphericalParam : public GenericParam
{

public:

    /**
     * @brief   Default empty constructor
     */
    SphericalParam();

    /**
     * @brief   function to read the parameter file in XML wrt dpd for convergence sphere
     * @param   <paramFile> parameter filename in
     */
    void readConvSphereXMLFile(const fs::path& paramConvFile,
                               CatalogData& cat);

private:

}; /* End of SphericalParam class */

} /* namespace LE3_2D_MASS_WL_SPHERICAL */

#endif
