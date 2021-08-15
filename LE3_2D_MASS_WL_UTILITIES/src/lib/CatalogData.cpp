/**
 * @file src/lib/CatalogData.cpp
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

#include "LE3_2D_MASS_WL_UTILITIES/CatalogData.h"

using namespace Euclid;
using namespace FitsIO;

static Elements::Logging logger = Elements::Logging::getLogger("CatalogData");

namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {

CatalogData::CatalogData(): m_catType("LENSMC") {}
CatalogData::CatalogData(std::string catType): m_catType(catType) {}

//std::vector<std::vector<double> > CatalogData::getCatalogData(fs::path& workdir, fs::path& InCatalog) {
void CatalogData::getCatalogData(fs::path& workdir, fs::path& InCatalog, std::vector<std::vector<double> >& data) {
   fs::path datadir {workdir / "data"};
   // Case 1: When input catalog file is in XML format (fetch catalog filename from xml file)
   if (false == checkFileType(workdir /InCatalog, Euclid::WeakLensing::TwoDMass::signFITS)) {
    if (true == checkFileType(workdir /InCatalog, Euclid::WeakLensing::TwoDMass::signXML)) {
     logger.info(" Input catalog file is in XML format..");
     logger.info() << "Using file " << workdir / InCatalog << " as DM input product";
     if (true == fileHasField((workdir /InCatalog).native(), "DpdTwoDMassLensMCCatalog")) {
       logger.info() << "catalogue type: FITS: galaxies with LensMC shear values";
       // Read filename of the FITS catalog from the input XML file
       DmInput in_xml = DmInput::readLensMCCatalogXMLFile(workdir / InCatalog);
       // Get fits Catalog filename
       m_inputCatalog = (datadir /  in_xml.getFitsCatalogFilename()).string();
       m_catType = in_xml.getMethodType();
       readCatalog(m_inputCatalog, data);
       logger.info() << "Shear selection type: " << m_catType;
       logger.info("Done reading Input Shear Catalog");
     } else if (true == fileHasField((workdir /InCatalog).native(), "DpdTwoDMassMomentsMLCatalog")) {
       logger.info() << "catalogue type: FITS: galaxies with MomentsML shear values";
       DmInput in_xml = DmInput::readMomentsMLCatalogXMLFile (workdir / InCatalog);
       m_inputCatalog = (datadir /  in_xml.getFitsCatalogFilename()).string();
       m_catType = in_xml.getMethodType();
       readCatalog(m_inputCatalog, data);
       logger.info() << "Shear selection type: " << m_catType;
       logger.info("Done reading Input Shear Catalog");
     } else if (true == fileHasField((workdir /InCatalog).native(), "DpdTwoDMassKSBCatalog")) {
       logger.info() << "catalogue type: FITS: galaxies with KSB shear values";
       DmInput in_xml = DmInput::readKSBCatalogXMLFile (workdir / InCatalog);
       m_inputCatalog = (datadir /  in_xml.getFitsCatalogFilename()).string();
       m_catType = in_xml.getMethodType();
       readCatalog(m_inputCatalog, data);
       logger.info() << "Shear selection type: " << m_catType;
       logger.info("Done reading Input Shear Catalog");
     } else if (true == fileHasField((workdir /InCatalog).native(), "DpdTwoDMassRegaussCatalog")) {
       logger.info() << "catalogue type: FITS: galaxies with Regauss shear values";
       DmInput in_xml = DmInput::readRegaussCatalogXMLFile (workdir / InCatalog);
       m_inputCatalog = (datadir /  in_xml.getFitsCatalogFilename()).string();
       m_catType = in_xml.getMethodType();
       readCatalog(m_inputCatalog, data);
       logger.info() << "Shear selection type: " << m_catType;
       logger.info("Done reading Input Shear Catalog");
     } else if (true == fileHasField((workdir /InCatalog).native(), "DpdWLLE2Catalog")) {
       logger.info() << "catalogue type: FITS: galaxies with shear values";
       DmInput in_xml = DmInput::readLE2CatalogXMLFile (workdir / InCatalog);
       m_inputCatalog = (datadir /  in_xml.getFitsCatalogFilename()).string();
       m_catType = (in_xml.getMethodType());
       readCatalog(m_inputCatalog, data);
       logger.info() << "Shear selection type: " << m_catType;
       logger.info("Done reading Input Shear Catalog");
     } else if (true == fileHasField((workdir /InCatalog).native(), "DpdTwoDMassClusterCatalog")) {
       logger.info() << "catalogue type: FITS: galaxies clusters";
       m_catType = "CLUSTER";
       // Read filename of the FITS catalog from the input XML file
       DmInput in_xml = DmInput::readClusterCatalogXMLFile(workdir /InCatalog);
       // Get fits Catalog filename
       m_inputCatalog = (datadir /  in_xml.getFitsCatalogFilename()).string();
       //data = readCatalog(m_inputCatalog);
       readCatalog(m_inputCatalog, data);
       logger.info("Done reading Input Cluster Catalog");
     }
    }
   }
 //return data;
}

void CatalogData::readClusterCatalog(const std::string& filename, std::vector<std::vector<double> >& Cat_Data) {
    MefFile b(filename, MefFile::Permission::Read);
    //const auto id = b.access<BintableHdu>(2).readColumn<double>("ID"); //EL_FitsIO library 2.2.0
    //EL_FitsIO library > 3.0.0 (HDU index 0 based unlike previous version)
    const auto id = b.access<BintableHdu>(1).readColumn<double>("ID");
    const auto ra = b.access<BintableHdu>(1).readColumn<double>("RA");
    const auto dec = b.access<BintableHdu>(1).readColumn<double>("DEC");
    const auto radius = b.access<BintableHdu>(1).readColumn<double>("RADIUS");
    const auto z = b.access<BintableHdu>(1).readColumn<double>("Z");
    const auto richness = b.access<BintableHdu>(1).readColumn<double>("RICHNESS");
    Cat_Data.push_back(id.vector());
    Cat_Data.push_back(ra.vector());
    Cat_Data.push_back(dec.vector());
    Cat_Data.push_back(z.vector());
    Cat_Data.push_back(radius.vector());
    Cat_Data.push_back(richness.vector());
}

void CatalogData::readShearCatalog(const std::string& filename, std::vector<std::vector<double> >& Cat_Data) {
    MefFile b(filename, MefFile::Permission::Read);
    //MefFile b(filename, MefFile::Permission::Edit);
    const auto& ext = b.access<BintableHdu>(1);
//    std::vector<std::string> incolname = ext.readColumnNames();
    std::vector<std::string> colname = getcolumnNames (m_catType);
/*    bool status = checkFileIsEuclidize(filename, colname);
    if (status == false) {
      int g1 = getIndexCol(incolname, "gamma1_noise");
      ext.renameColumn(incolname[g1], getColName(colname, "G1"));
      int g2 = getIndexCol(incolname, "gamma2_noise");
      ext.renameColumn(incolname[g2], getColName(colname, "G2"));
      ext.renameColumn(getColName(incolname, "ra"), getColName(colname, "RA"));
      ext.renameColumn(getColName(incolname, "dec"), getColName(colname, "DEC"));
      ext.renameColumn(getColName(incolname, "redshift"), getColName(colname, "PHZ_MEDIAN"));
      ext.renameColumn(getColName(incolname, "mask", "weight"), getColName(colname, "WEIGHT"));
    }*/
    std::vector<std::string> upColname = ext.readColumnNames();
    const auto gamma1 = ext.readColumn<double>(getColName(colname, "G1"));
    const auto gamma2 = ext.readColumn<double>(getColName(colname, "G2"));
    const auto weight = ext.readColumn<double>(getColName(colname, "WEIGHT"));
    const auto phz_median = ext.readColumn<double>(getColName(colname, "PHZ_MEDIAN"));
    const auto ra = ext.readColumn<double>(getColName(colname, "RA"));
    const auto dec = ext.readColumn<double>(getColName(colname, "DEC"));
//    logger.info() << "First value of SHE_LENSMC_UPDATED_RA = " << ra.vector()[0];

    Cat_Data.push_back(ra.vector());
    Cat_Data.push_back(dec.vector());
    if (ext.hasColumn("KAPPA")) {
      const auto kappa = ext.readColumn<double>("KAPPA");
      Cat_Data.push_back(kappa.vector());
    } else {
      std::vector <double> kappa;
      for (size_t it = 0; it < ra.vector().size(); it++) {
        kappa.push_back(0.);
      }
      Cat_Data.push_back(kappa);
    }

    Cat_Data.push_back(gamma1.vector());
    Cat_Data.push_back(gamma2.vector());
    Cat_Data.push_back(phz_median.vector());
    Cat_Data.push_back(weight.vector());

    if (ext.hasColumn(getColName(colname, "CORRECTION"))) {
      const auto phz_correction = ext.readColumn<double>(getColName(colname, "CORRECTION"));
      Cat_Data.push_back(phz_correction.vector());
    } else {
      std::vector <double> phz_correction;
      for (size_t it = 0; it < ra.vector().size(); it++) {
        phz_correction.push_back(1.);
      }
      Cat_Data.push_back(phz_correction);
    }
}

//std::vector<std::vector<double> > CatalogData::readCatalog(const std::string& filename) {
void CatalogData::readCatalog(const std::string& filename, std::vector<std::vector<double> >& Cat_Data) {
  //std::vector<std::vector<double> > Cat_Data;

  if (true == checkFileType(fs::path(filename), Euclid::WeakLensing::TwoDMass::signFITS)){

 // read the catalogue table
  if (m_catType == "CLUSTER"){
    logger.info("Reading Input Cluster Catalog..");
    readClusterCatalog(filename, Cat_Data);
  } else {
    logger.info() << "Reading Input Shear Catalog of shear type " << m_catType << ". . .";
    readShearCatalog(filename, Cat_Data);
  }
 }
  //return Cat_Data;
}

std::string CatalogData::getCatalogType() {
  return m_catType;
}

std::string CatalogData::getCatalogFitsName() {
  return m_inputCatalog;
}

} /* namespace TwoDMass */
} /* namespace WeakLensing */
} /* namespace Euclid */
