/**
 * @file src/lib/DmOutput.cpp
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

#include "LE3_2D_MASS_WL_UTILITIES/DmOutput.h"
#include "ST_DM_HeaderProvider/GenericHeaderProvider.h"
#include <fstream>
#include <boost/filesystem.hpp>


namespace fs = boost::filesystem;

// DM Output namespace and classes
using namespace dpd::le3::wl::twodmass::out::convergenceclusters;
using namespace dpd::le3::wl::twodmass::out::convergencepatch;
using namespace dpd::le3::wl::twodmass::out::convergencepatchestosphere;
using namespace dpd::le3::wl::twodmass::out::convergencesphere;
using namespace dpd::le3::wl::twodmass::out::peakcatalog;
using namespace Euclid::DataModel;
using Euclid::DataModel::GenericHeaderGenerator;

static Elements::Logging logger = Elements::Logging::getLogger("DmOutput");
namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {

 sys::genericHeader* getGenericHeader(const std::string& productType, const boost::filesystem::path& xmlFile){
	GenericHeaderGenerator* generator = new GenericHeaderGenerator("LE3Product");
    generator->setTagValue("ProdSDC", "SDC-FR");
    generator->setTagValue("SoftwareName", "2D-MASS-WL");
    generator->setTagValue("Curator", "VKansal");
    generator->changeProductType(productType);
    ELEMENTS_UNUSED sys::genericHeader *header = generator->generate();
  return header;
 }

 DmOutput::DmOutput(){
   logger.info() << "Create a DmOuput object for product type ";
 }

sys::dss::dataContainer DmOutput::createDataContainer(const boost::filesystem::path& fits_out_filename){
   logger.info() << "Create a data container pointing to the " << fits_out_filename << "......";
   boost::filesystem::path fits_file {fits_out_filename};
   sys::dss::dataContainer output(fits_file.filename().string(), "PROPOSED");
   return output;
}

// output 1
  void DmOutput::createNoisedPatchXml
          (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch> &product,
                                                            const boost::filesystem::path& fits_out_filename) {
     auto output = createDataContainer(fits_out_filename);
     pro::le3::wl::twodmass::twoDMassConvergencePatch NoisyMap(output, "le3.wl.2dmass.output.patchconvergence",
                                                              SWVersion);
     product->Data().NoisyConvergence
                      (pro::le3::wl::twodmass::twoDMassCollectConvergencePatch::NoisyConvergence_type{NoisyMap});
     auto& noisemaplist = product->Data().NoisyConvergence().get().DataContainer();
     noisemaplist.FileName(fits_out_filename.string());
     logger.info() << "Finished adding NoisyConvergence data container and file ";
  }

// output 2
  void DmOutput::createDenoisedPatchXml
       (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch> &product,
                                                            const boost::filesystem::path& fits_out_filename) {
     auto output = createDataContainer(fits_out_filename);
     pro::le3::wl::twodmass::twoDMassConvergencePatch DenoisedMap(output, "le3.wl.2dmass.output.patchconvergence",
                                                              SWVersion);
     product->Data().DenoisedConvergence
                 (pro::le3::wl::twodmass::twoDMassCollectConvergencePatch::DenoisedConvergence_type{DenoisedMap});
     auto& denoisedmaplist = product->Data().DenoisedConvergence().get().DataContainer();
     denoisedmaplist.FileName(fits_out_filename.string());
     logger.info() << "Finished adding DenoisedConvergence data container and file ";
  }

// output 3
  void DmOutput::createSNRPatchOutputXml
       (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatch::dpdTwoDMassConvergencePatch> &product,
                                                            const boost::filesystem::path& fits_out_filename) {
     auto output = createDataContainer(fits_out_filename);
     pro::le3::wl::twodmass::twoDMassSNRPatch snrMap(output, "le3.wl.2dmass.output.patchconvergence",
                                                              SWVersion);
     product->Data().SNRPatch
                 (pro::le3::wl::twodmass::twoDMassCollectConvergencePatch::SNRPatch_type{snrMap});
     auto& snrmaplist = product->Data().SNRPatch().get().DataContainer();
     snrmaplist.FileName(fits_out_filename.string());
     logger.info() << "Finished adding SNR Convergence data container and file ";
  }

// output 4
  void DmOutput::createNoisedSphereXml
       (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere> &product,
                                                             const boost::filesystem::path& fits_out_filename) {
     auto output = createDataContainer(fits_out_filename);
     pro::le3::wl::twodmass::twoDMassConvergenceSphere SphNoisedMap(output, "le3.wl.2dmass.output.sphereconvergence",
                                                                           SWVersion);
     product->Data().NoisyConvergence
                    (pro::le3::wl::twodmass::twoDMassCollectConvergenceSphere::NoisyConvergence_type{SphNoisedMap});
     auto& sphNoisedmaplist = product->Data().NoisyConvergence().get().DataContainer();
     sphNoisedmaplist.FileName(fits_out_filename.string());
     logger.info() << "Finished adding Noised Spherical Convergence data container and file ";
  }

// output 5
  void DmOutput::createDenoisedSphereXml
       (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere> &product,
                                                             const boost::filesystem::path& fits_out_filename) {
     auto output = createDataContainer(fits_out_filename);
     pro::le3::wl::twodmass::twoDMassConvergenceSphere SphDenoisedMap(output, "le3.wl.2dmass.output.sphereconvergence",
                                                                           SWVersion);
     product->Data().DenoisedConvergence
               (pro::le3::wl::twodmass::twoDMassCollectConvergenceSphere::DenoisedConvergence_type{SphDenoisedMap});
     auto& sphDenoisedmaplist = product->Data().DenoisedConvergence().get().DataContainer();
     sphDenoisedmaplist.FileName(fits_out_filename.string());
     logger.info() << "Finished adding De-Noised Spherical Convergence data container and file ";
  }

// output 6
  void DmOutput::createSNRSphereOutputXml
       (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere> &product,
                                                            const boost::filesystem::path& fits_out_filename) {
     auto output = createDataContainer(fits_out_filename);
     pro::le3::wl::twodmass::twoDMassSNRSphere SphsnrMap(output, "le3.wl.2dmass.output.sphereconvergence",
                                                                           SWVersion);
     product->Data().SNRMaps
                 (pro::le3::wl::twodmass::twoDMassCollectConvergenceSphere::SNRMaps_type{SphsnrMap});
     auto& sphsnrmaplist = product->Data().SNRMaps().get().DataContainer();
     sphsnrmaplist.FileName(fits_out_filename.string());
     logger.info() << "Finished adding SNR Spherical Convergence data container and file ";
  }

// output 7
  void DmOutput::createSphereGalCountXml
       (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere> &product,
                                                            const boost::filesystem::path& fits_out_filename) {
     auto output = createDataContainer(fits_out_filename);
     pro::le3::wl::twodmass::twoDMassGalCountSphere galCountMap(output, "le3.wl.2dmass.output.sphereconvergence",
                                                                           SWVersion);
     product->Data().GalCount
                 (pro::le3::wl::twodmass::twoDMassCollectConvergenceSphere::GalCount_type{galCountMap});
     auto& galcountmaplist = product->Data().GalCount().get().DataContainer();
     galcountmaplist.FileName(fits_out_filename.string());
     logger.info() << "Finished adding Spherical Galaxies count data container and file ";
  }

// output 8
  void DmOutput::createPeakOutputXml(const boost::filesystem::path& out_dm_xml_filename,
      const boost::filesystem::path& fits_out_filename) {
    logger.info() << "Creating PF output XML Peak Catalog product file " << out_dm_xml_filename << "...";

    logger.warn() << "Create a data container pointing to the" << fits_out_filename << "......";
    auto output = createDataContainer(fits_out_filename);
    std::ofstream out;
     const std::string produc = "DpdTwoDMassPeakCatalog";
     sys::genericHeader *header = getGenericHeader(produc, out_dm_xml_filename);
      pro::le3::wl::twodmass::twoDMassPeakCatalog out_cat(output, "le3.wl.2dmass.output.peakcatalog", SWVersion);
      pro::le3::wl::twodmass::twoDMassPeakCatalog data (out_cat);

      // Create the final object
      dpd::le3::wl::twodmass::out::peakcatalog::dpdTwoDMassPeakCatalog product(*header, data);
      //std::ofstream out;
      out.open (out_dm_xml_filename.string());
      dpd::le3::wl::twodmass::out::peakcatalog::DpdTwoDMassPeakCatalog(out, product);
     logger.info() << "Finished creating file " << out_dm_xml_filename;
  }

// output 9
  void DmOutput::createNoisedPatchtoSphereXml
      (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>
                        &product, const boost::filesystem::path& fits_out_filename) {
    auto output = createDataContainer(fits_out_filename);
    pro::le3::wl::twodmass::twoDMassConvergencePatchesToSphere NoisedMapP2S
                                               (output, "le3.wl.2dmass.output.patchessphereconvergence", SWVersion);
    product->Data().NoisyConvergence
           (pro::le3::wl::twodmass::twoDMassCollectConvergencePatchesToSphere::NoisyConvergence_type{NoisedMapP2S});
    auto& noisedmapP2Slist = product->Data().NoisyConvergence().get().DataContainer();
    noisedmapP2Slist.FileName(fits_out_filename.string());
    logger.info() << "Finished adding Noised patchestoSphere convergence data container and file";
  }

// output 10
  void DmOutput::createDenoisedPatchtoSphereXml
      (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>
                        &product, const boost::filesystem::path& fits_out_filename) {
    auto output = createDataContainer(fits_out_filename);
    pro::le3::wl::twodmass::twoDMassConvergencePatchesToSphere DenoisedMapP2S
                                  (output, "le3.wl.2dmass.output.patchessphereconvergence", SWVersion);
    product->Data().DenoisedConvergence
        (pro::le3::wl::twodmass::twoDMassCollectConvergencePatchesToSphere::DenoisedConvergence_type{DenoisedMapP2S});
    auto& denoisedmapP2Slist = product->Data().DenoisedConvergence().get().DataContainer();
    denoisedmapP2Slist.FileName(fits_out_filename.string());
    logger.info() << "Finished adding Denoised patchestoSphere convergence data container and file";
  }

// output 11
  void DmOutput::createSNRPatchtoSphereOutputXml
      (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>
                        &product, const boost::filesystem::path& fits_out_filename) {
	auto output = createDataContainer(fits_out_filename);
    pro::le3::wl::twodmass::twoDMassSNRPatchesToSphere snrMapP2S
                               (output, "le3.wl.2dmass.output.patchesspheresnr", SWVersion);
    product->Data().SNRMaps
           (pro::le3::wl::twodmass::twoDMassCollectConvergencePatchesToSphere::SNRMaps_type{snrMapP2S});
    auto& snrmapP2Slist = product->Data().SNRMaps().get().DataContainer();
    snrmapP2Slist.FileName(fits_out_filename.string());
    logger.info() << "Finished adding SNR patchestoSphere convergence data container and file";
  }

// output 12
  void DmOutput::createPatchtoSphereGalCountXml
      (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>
                        &product, const boost::filesystem::path& fits_out_filename) {
    auto output = createDataContainer(fits_out_filename);
    pro::le3::wl::twodmass::twoDMassGalCountPatchesToSphere galCountMap
                                    (output, "le3.wl.2dmass.output.patchesspheregalcount", SWVersion);
    product->Data().GalCount
               (pro::le3::wl::twodmass::twoDMassCollectConvergencePatchesToSphere::GalCount_type{galCountMap});
    auto& galcountmaplist = product->Data().GalCount().get().DataContainer();
    galcountmaplist.FileName(fits_out_filename.string());
    logger.info() << "Finished adding patchestoSphere GalCount data container and file";
  }

// output 13
  void DmOutput::createPatchtoSphereProjCenterPosXml
      (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>
                        &product, const boost::filesystem::path& fits_out_filename) {
	auto output = createDataContainer(fits_out_filename);
    pro::le3::wl::twodmass::twoDMassProjPosPatchesToSphere posi
                                  (output, "le3.wl.2dmass.output.patchesspherepos", SWVersion);
    product->Data().ProjCenterPos
             (pro::le3::wl::twodmass::twoDMassCollectConvergencePatchesToSphere::ProjCenterPos_type{posi});
    auto& poislist = product->Data().ProjCenterPos().get().DataContainer();
    poislist.FileName(fits_out_filename.string());;
    logger.info() << "Finished adding patchestoSphere Center Positions of patches data container and file";
  }

// output 14
  void DmOutput::createSingleClusterXml
       (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencesinglecluster::dpdTwoDMassConvergenceSingleCluster>
                        &product, const boost::filesystem::path& fits_out_filename) {
    auto output = createDataContainer(fits_out_filename);
    pro::le3::wl::twodmass::twoDMassConvergenceSingleClusterFile clusterconvergence
                                                    (output, "le3.wl.2dmass.output.clusterconvergence", SWVersion);
    product->Data().ClusterFile
                  (pro::le3::wl::twodmass::twoDMassConvergenceSingleCluster::ClusterFile_type{clusterconvergence});
    auto& clusterconvergencelist = product->Data().ClusterFile().get().DataContainer();
    clusterconvergencelist.FileName(fits_out_filename.string());
  }

// output 15
  void DmOutput::createClusterCatalogXml(const boost::filesystem::path& out_dm_xml_filename,
      const boost::filesystem::path& fits_out_filename) {
      logger.info() << "Creating PF output XML Peak Catalog product file " << out_dm_xml_filename << "...";

      logger.warn() << "Create a data container pointing to the" << fits_out_filename << "......";
      auto output = createDataContainer(fits_out_filename);
      std::ofstream out;
      const std::string produc = "DpdTwoDMassConvergenceClusters";
      sys::genericHeader *header = getGenericHeader(produc, out_dm_xml_filename);

      pro::le3::wl::twodmass::twoDMassConvergenceClusters data(output, fits_out_filename.string());

      dpd::le3::wl::twodmass::out::convergenceclusters::dpdTwoDMassConvergenceClusters
                                                                             product(*header, data);
      //std::ofstream out;
      out.open (out_dm_xml_filename.string());
      dpd::le3::wl::twodmass::out::convergenceclusters::DpdTwoDMassConvergenceClusters(out, product);
  }

// output 16
  void createSphereMCXml(std::unique_ptr<dpd::le3::wl::twodmass::out::convergencesphere::dpdTwoDMassConvergenceSphere>
            &product, const boost::filesystem::path& fits_out_filename) {
   //TODO
  }

// output 17
  void createPatchtoSphereMCXml
      (std::unique_ptr<dpd::le3::wl::twodmass::out::convergencepatchestosphere::dpdTwoDMassConvergencePatchesToSphere>
                        &product, const boost::filesystem::path& fits_out_filename) {
   //TODO
  }


} /* namespace TwoDMass */
} /* namespace WeakLensing */
} /* namespace Euclid */
