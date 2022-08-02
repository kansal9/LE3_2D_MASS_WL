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
 * @file LE3_2D_MASS_WL_PEAK_COUNT/PeakParam.h
 * @date 09/17/19
 * @author user
 */

#ifndef _LE3_2D_MASS_WL_PEAK_COUNT_PEAKPARAM_H
#define _LE3_2D_MASS_WL_PEAK_COUNT_PEAKPARAM_H

#include <string>
#include <boost/filesystem.hpp>
#include "ElementsKernel/Logging.h"
#include "ElementsKernel/ProgramHeaders.h"
#include <vector>
#include <sstream>

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

// Datamodel for INPUT products
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsPeakCountConvergence.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ParamsPeakCountMassAperture.h"
#include "ST_DataModelBindings/dictionary/pro/le3/wl/twodmass/euc-test-le3-wl-twodmass.h"

namespace LE3_2D_MASS_WL_PEAK_COUNT
{

/**
 * @class PeakParam
 * @brief
 *
 */
class PeakParam
{

public:

    /**
     * @brief   Default empty constructor
     */
    PeakParam();

    /**
     * @brief Destructor
     */
    virtual ~PeakParam() = default;

    /**
     * @brief   function to read the Peak count using convergence parameter file
     *          in XML wrt dpd
     * @param   <paramPeakConvergence> parameter filename in <string> format
     * @return  parameters from the file
     */
    void readWaveletPeakXML(const fs::path& paramPeakConvergence);

    /**
     * @brief   function to read the Peak count using mass aperture parameter
     *          file in XML wrt dpd
     * @param   <paramPeakMassAper> parameter filename in <string> format
     * @return  parameters from the file
     */
    void readMassApertureXML(const fs::path& paramPeakMassAper);

    /**
     * @brief   function to return number of Scales from Parameter file
     * @return  number of scales
     */
    int getNPeakScale() const;

    /**
     * @brief   function to return Peak Scale from Parameter file (Mass Aperture)
     * @return  Peak Scale
     */
    std::vector<double> getApPeakRadius() const;

private:
    int m_NPeakScale;

    std::vector<double> m_ApPeakRadius;

};
/* End of PeakParam class */

void readPeakParamFile(const fs::path& PeakParamFile, PeakParam &params);

} /* namespace LE3_2D_MASS_WL_PEAK_COUNT */

#endif
