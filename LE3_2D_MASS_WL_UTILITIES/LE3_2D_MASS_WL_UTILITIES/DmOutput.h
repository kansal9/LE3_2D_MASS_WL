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

// Datamodel OUTPUT products

#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-ConvergenceClusters.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-ConvergenceSingleCluster.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-ConvergencePatch.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-ConvergencePatchesToSphere.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-ConvergenceSphere.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/out/euc-test-le3-wl-twodmass-PeakCatalog.h"

// Generic header provider based on PF
//#include "HeaderProvider/GenericHeaderProvider.h"
#include "ST_DM_HeaderProvider/GenericHeaderProvider.h"

#include "ST_DataModelBindings/dictionary/pro/le3/wl/twodmass/euc-test-le3-wl-twodmass.h"
#include "ST_DataModelBindings/dictionary/sys/euc-test-sys.h"

namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {

/**
 * @class DmOutput
 * @brief This class writes the output products
 *
 */
class DmOutput {

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
  sys::dss::dataContainer createDataContainer(const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     helper function to convergencepatch output xml
   */
  static void helpWriteProduct(std::ofstream& out,
          std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch> &product){
	  dpd::le3::wl::twodmass::out::convergencepatch::DpdTwoDMassConvergencePatch(out, *(product.get()));
  }

  /**
   * @brief     helper function to convergencesphere output xml
   */
  static void helpWriteProduct(std::ofstream& out,
        std::unique_ptr<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere> &product){
	  dpd::le3::wl::twodmass::out::convergencesphere::DpdTwoDMassConvergenceSphere(out, *(product.get()));
  }

  /**
   * @brief     helper function to convergencepatchestosphere output xml
   */
  static void helpWriteProduct(std::ofstream& out, std::unique_ptr
      <dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere> &product){
    dpd::le3::wl::twodmass::out::convergencepatchestosphere::DpdTwoDMassConvergencePatchesToSphere
                                                                                        (out, *(product.get()));
  }

  /**
   * @brief     helper function to convergence single cluster output xml
   */
  static void helpWriteProduct(std::ofstream& out, std::unique_ptr
          <dpd::le3::wl::twodmass::out::convergencesinglecluster::dpdTwoDMassConvergenceSingleCluster> &product){
   dpd::le3::wl::twodmass::out::convergencesinglecluster::DpdTwoDMassConvergenceSingleCluster (out, *(product.get()));
  }

  /**
   * @brief     Add fits file (Cartesian Patch Noised convergence map) to XML product
   * @param[in] product, <dpdTwoDMassConvergencePatch> Convergence patch product
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createNoisedPatchXml
       (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch> &product,
                                                            const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     Add fits file (Cartesian Patch DeNoised convergence map) to XML product
   * @param[in] product, <dpdTwoDMassConvergencePatch> Convergence patch product
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createDenoisedPatchXml
       (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch> &product,
                                                            const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     Add fits file (Cartesian Patch SNR convergence map) to XML product
   * @param[in] product, <dpdTwoDMassConvergencePatch> Convergence patch product
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createSNRPatchOutputXml
       (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch> &product,
                                                            const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     Add fits file (noised sphere convergence map) to XML product
   * @param[in] product, <dpdTwoDMassConvergenceSphere> Convergence Sphere product
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createNoisedSphereXml
       (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere> &product,
                                                             const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     Add fits file (Denoised sphere convergence map) to XML product
   * @param[in] product, <dpdTwoDMassConvergenceSphere> Convergence Sphere product
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createDenoisedSphereXml
       (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere> &product,
                                                             const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     Add fits file (SNR sphere convergence map) to XML product
   * @param[in] product, <dpdTwoDMassConvergenceSphere> Convergence Sphere product
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createSNRSphereOutputXml
       (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere> &product,
                                                             const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     Add fits file (Galxies count on sphere convergence map) to XML product
   * @param[in] product, <dpdTwoDMassConvergenceSphere> Galxies count on sphere convergence map product
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createSphereGalCountXml
       (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere> &product,
                                                             const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     Add fits file (Cluster Catalogs zip) to XML product
   * @param[in] product, <dpdTwoDMassConvergenceClusters> Cluster Catalogs product
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  //void createClusterCatalogXml
  //     (std::unique_ptr<dpd::le3::wl::twodmass::out::convergenceclusters::dpdTwoDMassConvergenceClusters> &product,
  //                                                           const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     Add fits file (Cluster Catalog) to XML product
   * @param[in] product, <dpdTwoDMassConvergenceSingleClusters> Cluster Catalog product
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createSingleClusterXml
       (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencesinglecluster::dpdTwoDMassConvergenceSingleCluster>
                        &product, const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     Add fits file (Sphere Monte Carlo convergence maps) to XML product
   * @param[in] product, <dpdTwoDMassConvergenceSphere> Sphere convergence map
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createSphereMCXml(std::unique_ptr<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere>
            &product, const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     Add fits file (output peak catalog) to XML product
   * @param[in] product, <dpdTwoDMassPeakCatalog> output peak catalog
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  //void createPeakOutputXml(std::unique_ptr<dpd::le3::wl::twodmass::out::peakcatalog::dpdTwoDMassPeakCatalog>
  //                      &product, const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     Add fits file (Noised cartesian Patch convergence maps to sphere) to XML product
   * @param[in] product, <dpdTwoDMassConvergencePatchesToSphere> output spherical convergence from patches
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createNoisedPatchtoSphereXml
      (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>
                        &product, const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     Add fits file (DeNoised cartesian Patch convergence maps to sphere) to XML product
   * @param[in] product, <dpdTwoDMassConvergencePatchesToSphere> output spherical convergence from patches
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createDenoisedPatchtoSphereXml
      (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>
                        &product, const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     Add fits file (SNR cartesian Patch convergence maps to sphere) to XML product
   * @param[in] product, <dpdTwoDMassConvergencePatchesToSphere> output spherical convergence from patches
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createSNRPatchtoSphereOutputXml
      (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>
                        &product, const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     Add fits file (Galaxies count on cartesian Patch convergence maps to sphere) to XML product
   * @param[in] product, <dpdTwoDMassConvergencePatchesToSphere> output spherical Galaxies count from patches
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createPatchtoSphereGalCountXml
      (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>
                        &product, const boost::filesystem::path& fits_out_filename);
  /**
   * @brief     Add fits file (cartesian to Sphere Monte Carlo convergence map) to XML product
   * @param[in] product, <dpdTwoDMassConvergencePatchesToSphere> output spherical MC convergence map from patches
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createPatchtoSphereMCXml
      (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>
                        &product, const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     Add fits file (center position of cartesian convergence patches) to XML product
   * @param[in] product, <dpdTwoDMassConvergencePatchesToSphere> center position of cartesian convergence patches
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createPatchtoSphereProjCenterPosXml
      (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>
                        &product, const boost::filesystem::path& fits_out_filename);



//=================================================================================================
  /**
   * @brief     create output XML product for sphere Monte Carlo convergence maps
   * @param[in] out_dm_xml_filename, <filesystem::path> output xml filename including path
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  static void createSphereMCXml(const boost::filesystem::path& out_dm_xml_filename,
      const boost::filesystem::path& fits_out_filename, int NResamples=0);

  /**
   * @brief     create output XML product for Monte Carlo cartesian Patch convergence maps to sphere
   * @param[in] out_dm_xml_filename, <filesystem::path> output xml filename including path
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  static void createPatchtoSphereMCXml(const boost::filesystem::path& out_dm_xml_filename,
      const boost::filesystem::path& fits_out_filename, int NResamples=0);

//=================================================================================================

  /**
   * @brief     create output XML product for output peak catalog
   * @param[in] out_dm_xml_filename, <filesystem::path> output xml filename including path
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            fits output file in a data container of the XML file
   */
  void createPeakOutputXml(const boost::filesystem::path& out_dm_xml_filename,
      const boost::filesystem::path& fits_out_filename);

  /**
   * @brief     create output XML product for Cluster Catalog
   * @param[in] out_dm_xml_filename, <filesystem::path> output xml filename including path
   * @param[in] fits_out_filename, <filesystem::path> output fits filename including path
   * @return    None
   * @details   This method creates output product in xml format where it stores input
   *            tar output file in a data container of the XML file
   */
  void createClusterCatalogXml(const boost::filesystem::path& out_dm_xml_filename,
      const boost::filesystem::path& fits_out_filename);

private:

}; /* End of DmOutput class */

  /**
   * @brief  function to generate a generic header
   * @return generic header
   */
  sys::genericHeader* getGenericHeader(const std::string& productType, const boost::filesystem::path& xmlFile);

  /**
   * @brief     init output XML product for product_type
   * @param[in] out_dm_xml_filename, <filesystem::path> output xml filename including path
   * @param[in] product_type, <std::string> product type
   * @return    None
   */
template <typename Tdpd, typename Tpro, typename Tdata>
std::unique_ptr<Tdpd> initProduct(const std::string product_type,
											  	   const boost::filesystem::path& out_dm_xml_filename,
												   Tdata data_in){

  //logger.info() << "Creating PF output " << product_type << " XML product file " << out_dm_xml_filename << "...";

  sys::genericHeader *header = getGenericHeader(product_type, out_dm_xml_filename);
  Tpro data(data_in);
  Tdpd product(*header, data);
  return std::make_unique<Tdpd>(product);
}

  /**
   * @brief     write xml product
   * @param[in] out_dm_xml_filename, <filesystem::path> output xml filename including path
   * @return    None
   */
template <typename Tdpd>
void writeProduct(std::unique_ptr<Tdpd> &product, const boost::filesystem::path& out_dm_xml_filename){
	  std::ofstream out;
	  if(fs::exists(out_dm_xml_filename)){
		  out.open(out_dm_xml_filename.string(), std::ios_base::app);
	  } else {
		  out.open(out_dm_xml_filename.string());
	  }
	  DmOutput::helpWriteProduct(out, product);
	  //logger.info() << "Finished writing XML in " << out_dm_xml_filename;
}


} /* namespace TwoDMass */
} /* namespace WeakLensing */
} /* namespace Euclid */
#endif
