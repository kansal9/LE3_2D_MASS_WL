/**
 * @file LE3_2D_MASS_WL_UTILITIES/DmOutput.h
 * @date 10/13/20
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

#ifndef _DMOUTPUT_H
#define _DMOUTPUT_H

#include <string>
#include <boost/filesystem.hpp>
#include "ElementsKernel/Logging.h"

#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-ConvergenceClusters.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-ConvergencePatch.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-ConvergencePatchesToSphere.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-ConvergenceSphere.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-PeakCatalog.h"
#include "ST_DataModelBindings/dictionary/pro/le3/wl/twodmass/euc-test-le3-wl-twodmass.h"
#include "ST_DataModelBindings/dictionary/sys/euc-test-sys.h"
#include "ST_DataModelBindings/dictionary/bas/fit/euc-test-fit.h"
#include "ST_DM_HeaderProvider/GenericHeaderProvider.h"

using namespace dpd::le3::wl::twodmass::out::convergenceclusters;
using namespace dpd::le3::wl::twodmass::out::convergencepatch;
using namespace dpd::le3::wl::twodmass::out::convergencepatchestosphere;
using namespace dpd::le3::wl::twodmass::out::convergencesphere;
using namespace dpd::le3::wl::twodmass::out::peakcatalog;
using namespace pro::le3::wl::twodmass;
using namespace sys::dss;

using Euclid::DataModel::GenericHeaderGenerator;

namespace LE3_2D_MASS_WL_UTILITIES
{

/**
 * @class DmOutput
 * @brief This class writes the output products
 *
 */
class DmOutput
{

public:

    /**
     * @brief Destructor
     */
    virtual ~DmOutput() = default;

    /**
     * @brief     brief Constructor
     * @details   create output XML product
     */
    DmOutput();

    /**
     * @brief     Creates DataContainer
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     */
    sys::dss::dataContainer createDataContainer(const fs::path& fits_out_filename);

    /**
     * @brief     helper function to convergencepatch output xml
     */
    static void helpWriteProduct(std::ofstream& out,
            std::unique_ptr<dpdTwoDMassConvergencePatch> &product)
    {
        DpdTwoDMassConvergencePatch(out, *(product.get()));
    }

    /**
     * @brief     helper function to convergencesphere output xml
     */
    static void helpWriteProduct(std::ofstream& out,
            std::unique_ptr<dpdTwoDMassConvergenceSphere> &product)
    {
        DpdTwoDMassConvergenceSphere(out, *(product.get()));
    }

    /**
     * @brief     helper function to convergencepatchestosphere output xml
     */
    static void helpWriteProduct(std::ofstream& out,
            std::unique_ptr<dpdTwoDMassConvergencePatchesToSphere> &product)
    {
        DpdTwoDMassConvergencePatchesToSphere(out, *(product.get()));
    }

    /**
     * @brief     helper function to convergence single cluster output xml
     */
    /*
    static void helpWriteProduct(std::ofstream& out,
            std::unique_ptr<dpdTwoDMassConvergenceSingleCluster> &product)
    {
        DpdTwoDMassConvergenceSingleCluster(out, *(product.get()));
    }
    */

    /**
     * @brief     Add fits file (cartesian patch) to XML product
     * @param[in] product, <dpdTwoDMassConvergencePatch> convergence patch product
     * @param[in] datatype, type of the product to add (entry in outputType enum)
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     * @details   This method creates output product in xml format where it stores input
     *            fits output file in a data container of the XML file
     */
    void createPatchXml(
            std::unique_ptr<dpdTwoDMassConvergencePatch> &product,
            outputType datatype,
            const fs::path& fits_out_filename);

    /**
     * @brief     Add fits file (noised sphere convergence map) to XML product
     * @param[in] product, <dpdTwoDMassConvergenceSphere> Convergence Sphere product
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     * @details   This method creates output product in xml format where it stores input
     *            fits output file in a data container of the XML file
     */
    void createNoisedSphereXml(
            std::unique_ptr<dpdTwoDMassConvergenceSphere> &product,
            const fs::path& fits_out_filename);

    /**
     * @brief     Add fits file (Denoised sphere convergence map) to XML product
     * @param[in] product, <dpdTwoDMassConvergenceSphere> Convergence Sphere product
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     * @details   This method creates output product in xml format where it stores input
     *            fits output file in a data container of the XML file
     */
    void createDenoisedSphereXml(
            std::unique_ptr<dpdTwoDMassConvergenceSphere> &product,
            const fs::path& fits_out_filename);

    /**
     * @brief     Add fits file (SNR sphere convergence map) to XML product
     * @param[in] product, <dpdTwoDMassConvergenceSphere> Convergence Sphere product
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     * @details   This method creates output product in xml format where it stores input
     *            fits output file in a data container of the XML file
     */
    void createSNRSphereOutputXml(
            std::unique_ptr<dpdTwoDMassConvergenceSphere> &product,
            const fs::path& fits_out_filename);

    /**
     * @brief     Add fits file (Galxies count on sphere convergence map) to XML product
     * @param[in] product, <dpdTwoDMassConvergenceSphere> Galxies count on sphere convergence map product
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     * @details   This method creates output product in xml format where it stores input
     *            fits output file in a data container of the XML file
     */
    void createSphereGalCountXml(
            std::unique_ptr<dpdTwoDMassConvergenceSphere> &product,
            const fs::path& fits_out_filename);

    /**
     * @brief     Add fits file (Cluster Catalog) to XML product
     * @param[in] product, <dpdTwoDMassConvergenceSingleClusters> Cluster Catalog product
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     * @details   This method creates output product in xml format where it stores input
     *            fits output file in a data container of the XML file
     */
    /*
    void createSingleClusterXml(
            std::unique_ptr<dpdTwoDMassConvergenceSingleCluster> &product,
            const fs::path& fits_out_filename);
    */

    /**
     * @brief     Add fits file (Sphere Monte Carlo convergence maps) to XML product
     * @param[in] product, <dpdTwoDMassConvergenceSphere> Sphere convergence map
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     * @details   This method creates output product in xml format where it stores input
     *            fits output file in a data container of the XML file
     */
    void createSphereMCXml(
            std::unique_ptr<dpdTwoDMassConvergenceSphere> &product,
            const fs::path& fits_out_filename);

    /**
     * @brief     Add fits file (Noised cartesian Patch convergence maps to sphere) to XML product
     * @param[in] product, <dpdTwoDMassConvergencePatchesToSphere> output spherical convergence from patches
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     * @details   This method creates output product in xml format where it stores input
     *            fits output file in a data container of the XML file
     */
    void createNoisedPatchtoSphereXml(
            std::unique_ptr<dpdTwoDMassConvergencePatchesToSphere> &product,
            const fs::path& fits_out_filename);

    /**
     * @brief     Add fits file (DeNoised cartesian Patch convergence maps to sphere) to XML product
     * @param[in] product, <dpdTwoDMassConvergencePatchesToSphere> output spherical convergence from patches
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     * @details   This method creates output product in xml format where it stores input
     *            fits output file in a data container of the XML file
     */
    void createDenoisedPatchtoSphereXml(
            std::unique_ptr<dpdTwoDMassConvergencePatchesToSphere> &product,
            const fs::path& fits_out_filename);

    /**
     * @brief     Add fits file (SNR cartesian Patch convergence maps to sphere) to XML product
     * @param[in] product, <dpdTwoDMassConvergencePatchesToSphere> output spherical convergence from patches
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     * @details   This method creates output product in xml format where it stores input
     *            fits output file in a data container of the XML file
     */
    void createSNRPatchtoSphereOutputXml(
            std::unique_ptr<dpdTwoDMassConvergencePatchesToSphere> &product,
            const fs::path& fits_out_filename);

    /**
     * @brief     Add fits file (Galaxies count on cartesian Patch convergence maps to sphere) to XML product
     * @param[in] product, <dpdTwoDMassConvergencePatchesToSphere> output spherical Galaxies count from patches
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     * @details   This method creates output product in xml format where it stores input
     *            fits output file in a data container of the XML file
     */
    void createPatchtoSphereGalCountXml(
            std::unique_ptr<dpdTwoDMassConvergencePatchesToSphere> &product,
            const fs::path& fits_out_filename);
    /**
     * @brief     Add fits file (cartesian to Sphere Monte Carlo convergence map) to XML product
     * @param[in] product, <dpdTwoDMassConvergencePatchesToSphere> output spherical MC convergence map from patches
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     * @details   This method creates output product in xml format where it stores input
     *            fits output file in a data container of the XML file
     */
    void createPatchtoSphereMCXml(
            std::unique_ptr<dpdTwoDMassConvergencePatchesToSphere> &product,
            const fs::path& fits_out_filename);

    /**
     * @brief     Add fits file (center position of cartesian convergence patches) to XML product
     * @param[in] product, <dpdTwoDMassConvergencePatchesToSphere> center position of cartesian convergence patches
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     * @details   This method creates output product in xml format where it stores input
     *            fits output file in a data container of the XML file
     */
    void createPatchtoSphereProjCPosXml(
            std::unique_ptr<dpdTwoDMassConvergencePatchesToSphere> &product,
            const fs::path& fits_out_filename);

    /**
     * @brief     create output XML product for output peak catalog
     * @param[in] out_dm_xml_filename, <filesystem::path> output xml filename including path
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     * @details   This method creates output product in xml format where it stores input
     *            fits output file in a data container of the XML file
     */
    void createPeakOutputXml(const boost::filesystem::path& out_dm_xml_filename,
            const fs::path& fits_out_filename);

    /**
     * @brief     create output XML product for Cluster Catalog
     * @param[in] out_dm_xml_filename, <filesystem::path> output xml filename including path
     * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
     * @return    None
     * @details   This method creates output product in xml format where it stores input
     *            tar output file in a data container of the XML file
     */
    void createClusterCatalogXml(
            const fs::path& out_dm_xml_filename,
            const fs::path& fits_out_filename);
private:

};
/* End of DmOutput class */

/**
 * @brief  function to generate a generic header
 * @return generic header
 */
sys::genericHeader* getGenericHeader(const std::string& productType);

/**
 * @brief     init output XML product for product_type
 * @param[in] out_dm_xml_filename, <filesystem::path> output xml filename including path
 * @param[in] product_type, <std::string> product type
 * @return    None
 */
template<typename Tdpd, typename Tpro, typename Tdata>
std::unique_ptr<Tdpd> initProduct(const std::string product_type,
                                  Tdata data_in)
{
    auto header = getGenericHeader(product_type);
    Tpro data(data_in);
    Tdpd product(*header, data);
    delete header;
    return std::make_unique<Tdpd>(product);
}

/**
 * @brief     write xml product
 * @param[in] out_dm_xml_filename, <filesystem::path> output xml filename including path
 * @return    None
 */
template<typename Tdpd>
void writeProduct(std::unique_ptr<Tdpd> &product,
                  const fs::path& out_dm_xml_filename)
{
    std::ofstream out;
    if (fs::exists(out_dm_xml_filename))
    {
        out.open(out_dm_xml_filename.string(), std::ios_base::app);
    }
    else
    {
        out.open(out_dm_xml_filename.string());
    }
    DmOutput::helpWriteProduct(out, product);
}

} // namespace LE3_2D_MASS_WL_UTILITIES

#endif
