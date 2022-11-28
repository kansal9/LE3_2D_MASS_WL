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

#include <fstream>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

static Elements::Logging logger = Elements::Logging::getLogger("DmOutput");

namespace LE3_2D_MASS_WL_UTILITIES
{

using DataContainer_type = sys::dss::dataContainer;
using NoisyConvergence_type = twoDMassCollectConvergencePatch::NoisyConvergence_type;
using DenoisedConvergence_type = twoDMassCollectConvergencePatch::DenoisedConvergence_type;
using SNRPatch_type = twoDMassCollectConvergencePatch::SNRPatch_type;

sys::genericHeader* getGenericHeader(const std::string& productType)
{
    GenericHeaderGenerator generator("LE3Product");
    generator.setTagValue("ProdSDC", "SDC-FR");
    generator.setTagValue("SoftwareName", "2D-MASS-WL");
    generator.setTagValue("Curator", "VKansal");
    generator.changeProductType(productType);
    return generator.generate();
}

DmOutput::DmOutput()
{
    logger.info() << "Create a DmOuput object for product type ";
}

dataContainer DmOutput::createDataContainer(const fs::path& fits_out_filename)
{
    logger.info() << "Create a data container pointing to the "
                  << fits_out_filename << "......";
    dataContainer output(fits_out_filename.filename().string(), "PROPOSED");
    return output;
}

void DmOutput::createPatchXml(
        std::unique_ptr<dpdTwoDMassConvergencePatch> &product,
        outputType datatype,
        const fs::path& fits_out_filename)
{
    if(!fs::exists(fits_out_filename))
    {
        logger.debug() << "Not creating product for non-existant file " << fits_out_filename;
        return;
    }

    auto output = createDataContainer(fits_out_filename);
    std::unique_ptr<DataContainer_type> maplist_ptr;
    if(datatype == outputType::NoisedPatch)
    {
        twoDMassConvergencePatch map(output, "le3.wl.2dmass.output.patchconvergence");
        product->Data().NoisyConvergence(NoisyConvergence_type{ map });
        maplist_ptr = std::make_unique<DataContainer_type>(
                product->Data().NoisyConvergence().get().DataContainer());

    }
    else if(datatype == outputType::DenoisedPatch)
    {
        twoDMassConvergencePatch map(output, "le3.wl.2dmass.output.patchconvergence");
        product->Data().DenoisedConvergence(DenoisedConvergence_type{ map });
        maplist_ptr = std::make_unique<DataContainer_type>(
                product->Data().DenoisedConvergence().get().DataContainer());
    }
    else if(datatype == outputType::SNRPatch)
    {
        twoDMassSNRPatch map(output, "le3.wl.2dmass.output.patchsnr");
        product->Data().SNRPatch(SNRPatch_type{ map });
        maplist_ptr = std::make_unique<DataContainer_type>(
                product->Data().SNRPatch().get().DataContainer());
    }
    maplist_ptr->FileName(fits_out_filename.filename().string());
    logger.info() << "Added data container and file";
}

// output 4
void DmOutput::createNoisedSphereXml(
        std::unique_ptr<dpdTwoDMassConvergenceSphere> &product,
        const fs::path& fits_out_filename)
{
    auto output = createDataContainer(fits_out_filename);
    twoDMassConvergenceSphere SphNoisedMap(output,
            "le3.wl.2dmass.output.sphereconvergence");
    product->Data().NoisyConvergence(
            twoDMassCollectConvergenceSphere::NoisyConvergence_type
            { SphNoisedMap });
    auto& sphNoisedmaplist =
            product->Data().NoisyConvergence().get().DataContainer();
    sphNoisedmaplist.FileName(fits_out_filename.string());
    logger.info()
            << "Added Noised Spherical Convergence data container and file ";
}

// output 5
void DmOutput::createDenoisedSphereXml(
        std::unique_ptr<dpdTwoDMassConvergenceSphere> &product,
        const fs::path& fits_out_filename)
{
    auto output = createDataContainer(fits_out_filename);
    twoDMassConvergenceSphere SphDenoisedMap(output,
            "le3.wl.2dmass.output.sphereconvergence");
    product->Data().DenoisedConvergence(
            twoDMassCollectConvergenceSphere::DenoisedConvergence_type
            { SphDenoisedMap });
    auto& sphDenoisedmaplist =
            product->Data().DenoisedConvergence().get().DataContainer();
    sphDenoisedmaplist.FileName(fits_out_filename.string());
    logger.info()
            << "Added De-Noised Spherical Convergence data container and file ";
}

// output 6
void DmOutput::createSNRSphereOutputXml(
        std::unique_ptr<dpdTwoDMassConvergenceSphere> &product,
        const fs::path& fits_out_filename)
{
    auto output = createDataContainer(fits_out_filename);
    twoDMassSNRSphere SphsnrMap(output,
            "le3.wl.2dmass.output.sphereconvergence");
    product->Data().SNRMaps(twoDMassCollectConvergenceSphere::SNRMaps_type
            { SphsnrMap });
    auto& sphsnrmaplist = product->Data().SNRMaps().get().DataContainer();
    sphsnrmaplist.FileName(fits_out_filename.string());
    logger.info()
            << "Added SNR Spherical Convergence data container and file ";
}

// output 7
void DmOutput::createSphereGalCountXml(
        std::unique_ptr<dpdTwoDMassConvergenceSphere> &product,
        const fs::path& fits_out_filename)
{
    auto output = createDataContainer(fits_out_filename);
    twoDMassGalCountSphere galCountMap(output,
            "le3.wl.2dmass.output.sphereconvergence");
    product->Data().GalCount(twoDMassCollectConvergenceSphere::GalCount_type
            { galCountMap });
    auto& galcountmaplist = product->Data().GalCount().get().DataContainer();
    galcountmaplist.FileName(fits_out_filename.string());
    logger.info()
            << "Added Spherical Galaxies count data container and file ";
}

// output 8
void DmOutput::createPeakOutputXml(
        const fs::path& out_dm_xml_filename,
        const fs::path& fits_out_filename)
{
    logger.info() << "Creating PF output XML Peak Catalog product file "
                  << out_dm_xml_filename << "...";

    logger.warn() << "Create a data container pointing to the"
                  << fits_out_filename << "......";
    dataContainer output = createDataContainer(fits_out_filename);
    twoDMassPeakCatalog data(output, "le3.wl.2dmass.output.peakcatalog");
    const std::string produc = "DpdTwoDMassPeakCatalog";
    auto header = getGenericHeader(produc);
    dpdTwoDMassPeakCatalog product(*header, data);
    std::ofstream out;
    out.open(out_dm_xml_filename.string());
    DpdTwoDMassPeakCatalog(out, product);
    logger.info() << "Finished creating file " << out_dm_xml_filename;
    delete header;
}

// output 9
void DmOutput::createNoisedPatchtoSphereXml(
        std::unique_ptr<dpdTwoDMassConvergencePatchesToSphere> &product,
        const fs::path& fits_out_filename)
{
    auto output = createDataContainer(fits_out_filename);
    twoDMassConvergencePatchesToSphere NoisedMapP2S(
            output, "le3.wl.2dmass.output.patchessphereconvergence");
    product->Data().NoisyConvergence(
            twoDMassCollectConvergencePatchesToSphere::NoisyConvergence_type
            { NoisedMapP2S });
    auto& noisedmapP2Slist =
            product->Data().NoisyConvergence().get().DataContainer();
    noisedmapP2Slist.FileName(fits_out_filename.string());
    logger.info()
            << "Added Noised patchestoSphere convergence data container and file";
}

// output 10
void DmOutput::createDenoisedPatchtoSphereXml(
        std::unique_ptr<dpdTwoDMassConvergencePatchesToSphere> &product,
        const fs::path& fits_out_filename)
{
    auto output = createDataContainer(fits_out_filename);
    twoDMassConvergencePatchesToSphere DenoisedMapP2S(
            output, "le3.wl.2dmass.output.patchessphereconvergence");
    product->Data().DenoisedConvergence(
            twoDMassCollectConvergencePatchesToSphere::DenoisedConvergence_type
            { DenoisedMapP2S });
    auto& denoisedmapP2Slist =
            product->Data().DenoisedConvergence().get().DataContainer();
    denoisedmapP2Slist.FileName(fits_out_filename.string());
    logger.info()
            << "Added Denoised patchestoSphere convergence data container and file";
}

// output 11
void DmOutput::createSNRPatchtoSphereOutputXml(
        std::unique_ptr<dpdTwoDMassConvergencePatchesToSphere> &product,
        const fs::path& fits_out_filename)
{
    auto output = createDataContainer(fits_out_filename);
    twoDMassSNRPatchesToSphere snrMapP2S(output,
            "le3.wl.2dmass.output.patchesspheresnr");
    product->Data().SNRMaps(
            twoDMassCollectConvergencePatchesToSphere::SNRMaps_type
            { snrMapP2S });
    auto& snrmapP2Slist = product->Data().SNRMaps().get().DataContainer();
    snrmapP2Slist.FileName(fits_out_filename.string());
    logger.info()
            << "Added SNR patchestoSphere convergence data container and file";
}

// output 12
void DmOutput::createPatchtoSphereGalCountXml(
        std::unique_ptr<dpdTwoDMassConvergencePatchesToSphere> &product,
        const fs::path& fits_out_filename)
{
    auto output = createDataContainer(fits_out_filename);
    twoDMassGalCountPatchesToSphere galCountMap(output,
            "le3.wl.2dmass.output.patchesspheregalcount");
    product->Data().GalCount(
            twoDMassCollectConvergencePatchesToSphere::GalCount_type
            { galCountMap });
    auto& galcountmaplist = product->Data().GalCount().get().DataContainer();
    galcountmaplist.FileName(fits_out_filename.string());
    logger.info()
            << "Added patchestoSphere GalCount data container and file";
}

// output 13
void DmOutput::createPatchtoSphereProjCPosXml(
        std::unique_ptr<dpdTwoDMassConvergencePatchesToSphere> &product,
        const fs::path& fits_out_filename)
{
    auto output = createDataContainer(fits_out_filename);
    twoDMassProjPosPatchesToSphere posi(output,
            "le3.wl.2dmass.output.patchesspherepos");
    product->Data().ProjCenterPos(
            twoDMassCollectConvergencePatchesToSphere::ProjCenterPos_type
            { posi });
    auto& poislist = product->Data().ProjCenterPos().get().DataContainer();
    poislist.FileName(fits_out_filename.string());
    logger.info()
            << "Added patchestoSphere Center Positions of patches data container and file";
}

// output 14
/*
void DmOutput::createSingleClusterXml(
        std::unique_ptr<dpdTwoDMassConvergenceSingleCluster> &product,
        const fs::path& fits_out_filename)
{
    auto output = createDataContainer(fits_out_filename);
    twoDMassConvergenceSingleClusterFile clusterconvergence(
            output, "le3.wl.2dmass.output.clusterconvergence");
    product->Data().ClusterFile(
            twoDMassConvergenceSingleCluster::ClusterFile_type
            { clusterconvergence });
    auto& clusterconvergencelist =
            product->Data().ClusterFile().get().DataContainer();
    clusterconvergencelist.FileName(fits_out_filename.string());
    logger.info()
            << "Added create single cluster data container and file";
}
*/

// output 15
/*
void DmOutput::createClusterCatalogXml(
        const fs::path& out_dm_xml_filename,
        const fs::path& fits_out_filename)
{
    logger.info() << "Creating PF output XML Peak Catalog product file "
            << out_dm_xml_filename << "...";
    logger.warn() << "Create a data container pointing to the"
            << fits_out_filename << "......";
    auto output = createDataContainer(fits_out_filename);
    twoDMassConvergenceClusters data(output, fits_out_filename.string());
    const std::string produc = "DpdTwoDMassConvergenceClusters";
    auto header = getGenericHeader(produc);
    dpdTwoDMassConvergenceClusters product(*header, data);
    std::ofstream out;
    out.open(out_dm_xml_filename.string());
    DpdTwoDMassConvergenceClusters(out, product);
    logger.info() << "Added cluster catalogue data container and file";
    delete header;
}
*/

// output 16
void DmOutput::createSphereMCXml(
        std::unique_ptr<dpdTwoDMassConvergenceSphere> &product,
        const fs::path& fits_out_filename)
{
    auto output = createDataContainer(fits_out_filename);
    twoDMassMCSphere mcMap(output,
            "le3.wl.2dmass.output.sphereconvergence");
    twoDMassCollectConvergenceSphere::MonteCarloMap_sequence &mcmap(
            product->Data().MonteCarloMap());
    mcmap.push_back(twoDMassCollectConvergenceSphere::MonteCarloMap_type(mcMap));
    logger.info() << "Added Sphere MC convergence data container and file";
}

// output 17
void DmOutput::createPatchtoSphereMCXml(
        std::unique_ptr<dpdTwoDMassConvergencePatchesToSphere> &product,
        const fs::path& fits_out_filename)
{
    auto output = createDataContainer(fits_out_filename);
    twoDMassMCPatchesToSphere MCMap(output,
            "le3.wl.2dmass.output.patchessphereconvergence");
    twoDMassCollectConvergencePatchesToSphere::MonteCarloMap_sequence &mcmap(
            product->Data().MonteCarloMap());
    mcmap.push_back(
          twoDMassCollectConvergencePatchesToSphere::MonteCarloMap_type(MCMap));
    logger.info()
            << "Added patches to Sphere MC convergence data container and file";
}

} // namespace LE3_2D_MASS_WL_UTILITIES
