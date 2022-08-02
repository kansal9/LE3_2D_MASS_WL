/**
 * @file LE3_2D_MASS_WL_CARTESIAN/PatchesToSphereAlgo.h
 * @date 08/17/21
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

#ifndef _PATCHESTOSPHEREALGO_H
#define _PATCHESTOSPHEREALGO_H

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "ElementsKernel/Temporary.h"

#include "LE3_2D_MASS_WL_CARTESIAN/CartesianAlgoKS.h"
#include "LE3_2D_MASS_WL_CARTESIAN/Projection.h"
#include "LE3_2D_MASS_WL_UTILITIES/DmOutput.h"

#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-ConvergencePatchesToSphere.h"

using namespace dpd::le3::wl::twodmass::out::convergencepatchestosphere;

namespace fs = boost::filesystem;

namespace LE3_2D_MASS_WL_CARTESIAN
{

/**
 * @class PatchesToSphereAlgo
 * @brief
 *
 */
class PatchesToSphereAlgo
{

public:

    /**
     * @brief Constructor
     */
    PatchesToSphereAlgo(CartesianParam &params);

    /**
     * @brief Destructor
     */
    virtual ~PatchesToSphereAlgo() = default;

    /**
     * @brief     This function sets name of output file for noisy convergence maps and perform mass mapping
     * @param     <workdir>, <boost::filesystem::path> work directory
     * @param     <InShearMap>, <boost::filesystem::path> name of the json file which contains name of the shearMaps
     * @param     <outputMaps>, <boost::filesystem::path> json file containing names of the noisy convergence maps
     * @return    <bool> true if convergence maps well created
     */
    bool getNoisyConvergenceMaps(fs::path& workdir, fs::path& InShearMap,
            fs::path& outputMaps);

    /**
     * @brief     This function sets name of output file for denoised convergence maps and perform mass mapping
     * @param     <workdir>, <boost::filesystem::path> work directory
     * @param     <InShearMap>, <boost::filesystem::path> name of the json file which contains name of the shearMaps
     * @param     <outputMaps>, <boost::filesystem::path> json file containing names of the denoised convergence maps
     * @return    <bool> true if convergence maps well created
     */
    bool getDeNoisyConvergenceMaps(fs::path& workdir, fs::path& InShearMap,
            fs::path& outputMaps);

    void getPatchData(fs::path &convMaps, fs::path &workdir,
            const std::string& outCenterPatchPositionFileFits,
            std::vector<std::vector<double>>& data);

    bool getCenterPatchPositionsFile(const std::string& outFileFits,
            std::vector<double>& ra_ctr, std::vector<double>& dec_ctr,
            PatchDef& m_CB);

    void getHealpixFormatMap(std::vector<std::vector<double>>& data,
            Healpix_Map<double>& mapE, Healpix_Map<double>& mapB,
            Healpix_Map<double>& mapGalCount);

    /**
     @brief  This method write Fits BinTable (according to DataModel) using input healpix map and filename
     @param  <filename> name of the map in fits format
     @param  <map> healpix map
     @param  <colname> column name in table
     @retun  true if map is wriiten correctly
     */
    bool write_Map(const std::string& filename, Healpix_Map<double>& map,
            const std::string &colname);

    /**
     * @brief   writes records to the Primary header
     * @param[in] Primary hdu to which records need to be written
     * This method writes records to the Primary header
     **/
    void writePrimaryHeader(const Hdu &hdu);

    /**
     * @brief   writes records to the given header
     * @param[in] hdu (other than primary) to which records need to be written
     * This method writes records to the given header
     **/
    void writeHdu(const Hdu &hdu, int Nside);

private:
    /**
     *  @brief <m_cartesianParam>, CartesianParam object with catalog parameters
     */
    CartesianParam m_cartesianParam;

    /**
     *  @brief <hb>, Healpix_Base object
     */
    Healpix_Base hb;

};
// End of PatchesToSphereAlgo class

}// namespace LE3_2D_MASS_WL_CARTESIAN

#endif
