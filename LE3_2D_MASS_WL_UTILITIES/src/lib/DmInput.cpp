/**
 * @file src/lib/DmInput.cpp
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

#include "LE3_2D_MASS_WL_UTILITIES/DmInput.h"


namespace fs = boost::filesystem;
// DM Input namespace and class
//using namespace dpd::le3::wl::twodmass::inp::shearcatalogs;
using namespace dpd::le3::wl::twodmass::inp::ksbcatalog;
using namespace dpd::le3::wl::twodmass::inp::clustercatalogs;

static Elements::Logging logger = Elements::Logging::getLogger("DmInput");

namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {

 // DmInput::DmInput(){}
 //std::vector< std::pair <double, double> > DmInput::m_vertices;

 DmInput DmInput::readCatalogXMLFile(const Euclid::WeakLensing::TwoDMass::catalogType catType,
                                     const fs::path& in_xml_catalog_file) {
  logger.info() << "Getting information from input product XML file " << in_xml_catalog_file << " ...";
  // Parse the XML file and create the binding object
  logger.debug() << "Parsing file " << in_xml_catalog_file << " ...";
  //Define Variable to store Catalog filename
  fs::path catalog_file {};
  std::string methodType = "LENSMC";
  try {
   if (catType == Euclid::WeakLensing::TwoDMass::catalogType::KSBCatalog) {
     auto in_xml = dpd::le3::wl::twodmass::inp::ksbcatalog::DpdTwoDMassKSBCatalog
                       (in_xml_catalog_file.string(), xml_schema::flags::dont_validate);
     logger.debug("Binding object created successfully");

    // Get the name of the FITS file containing the catalog
     catalog_file = in_xml->Data().ShearCatalog().DataContainer().FileName();
     logger.info()<< "filename: " << in_xml->Data().ShearCatalog().DataContainer().FileName();    
      methodType = "KSB";
    /* bas::imp::stc::polygonType::Vertex_sequence& vertices(in_xml->Data().SpatialCoverage().Polygon().Vertex());
       for (bas::imp::stc::polygonType::Vertex_iterator i (vertices.begin()); i!=vertices.end(); ++i) {
       bas::imp::stc::vertexType& ptr (*i);
       logger.info()<< "Vertex X1: " << ptr.Position().C1();
       logger.info()<< "Vertex Y1: " << ptr.Position().C2();
       m_vertices.push_back(std::make_pair(ptr.Position().C1(),ptr.Position().C2()));
     }*/
  //  m_catalogPath = in_xml->Data().CatalogDescription().PathToCatalogDefinition();
    // logger.info()<< "Path to catalog: " << in_xml->Data().CatalogDescription().PathToCatalogFile();
   }

   if (catType == Euclid::WeakLensing::TwoDMass::catalogType::clusterCatalog) {
     auto in_xml = dpd::le3::wl::twodmass::inp::clustercatalogs::DpdTwoDMassClusterCatalog
                              (in_xml_catalog_file.string(), xml_schema::flags::dont_validate);
     logger.debug("Binding object created successfully");
      methodType = "CLUSTER";
    // Get the name of the FITS file containing the catalog
     catalog_file = in_xml->Data().ClusterCatalog().DataContainer().FileName();
     logger.info()<< "filename: " << in_xml->Data().ClusterCatalog().DataContainer().FileName();    
   }

   if (catType == Euclid::WeakLensing::TwoDMass::catalogType::LensMCCatalog) {
     auto in_xml = dpd::le3::wl::twodmass::inp::lensmccatalog::DpdTwoDMassLensMCCatalog
                              (in_xml_catalog_file.string(), xml_schema::flags::dont_validate);
     logger.debug("Binding object created successfully");
      methodType = "LENSMC";
    // Get the name of the FITS file containing the catalog
     catalog_file = in_xml->Data().ShearCatalog().DataContainer().FileName();
     logger.info()<< "filename: " << in_xml->Data().ShearCatalog().DataContainer().FileName();    
   }

   if (catType == Euclid::WeakLensing::TwoDMass::catalogType::MomentsMLCatalog) {
     auto in_xml = dpd::le3::wl::twodmass::inp::momentsmlcatalog::DpdTwoDMassMomentsMLCatalog
                              (in_xml_catalog_file.string(), xml_schema::flags::dont_validate);
     logger.debug("Binding object created successfully");
      methodType = "MOMENTSML";
    // Get the name of the FITS file containing the catalog
     catalog_file = in_xml->Data().ShearCatalog().DataContainer().FileName();
     logger.info()<< "filename: " << in_xml->Data().ShearCatalog().DataContainer().FileName();    
   }

   if (catType == Euclid::WeakLensing::TwoDMass::catalogType::RegaussCatalog) {
     auto in_xml = dpd::le3::wl::twodmass::inp::regausscatalog::DpdTwoDMassRegaussCatalog
                              (in_xml_catalog_file.string(), xml_schema::flags::dont_validate);
     logger.debug("Binding object created successfully");
      methodType = "REGAUSS";
    // Get the name of the FITS file containing the catalog
     catalog_file = in_xml->Data().ShearCatalog().DataContainer().FileName();
     logger.info()<< "filename: " << in_xml->Data().ShearCatalog().DataContainer().FileName();    
   }
   if (catType == Euclid::WeakLensing::TwoDMass::catalogType::VisibilityMask) {
     auto in_xml = dpd::le3::wl::twodmass::inp::visibilitymask::DpdTwoDMassVisibilityMask
                              (in_xml_catalog_file.string(), xml_schema::flags::dont_validate);
     logger.debug("Binding object created successfully");
      methodType = "";
    // Get the name of the FITS file containing the catalog
     catalog_file = in_xml->Data().DataContainer().FileName();
     logger.info()<< "filename: " << in_xml->Data().DataContainer().FileName();    
   }
   if (catType == Euclid::WeakLensing::TwoDMass::catalogType::LE2Catalog) {
     auto in_xml = dpd::le3::wl::inp::inputle2catalog::DpdWLLE2Catalog
                              (in_xml_catalog_file.string(), xml_schema::flags::dont_validate);
     logger.debug("Binding object created successfully");
    // Get the name of the FITS file containing the catalog
     catalog_file = in_xml->Data().LE2Catalog().DataContainer().FileName();
     logger.info()<< "filename: " << in_xml->Data().LE2Catalog().DataContainer().FileName();   

     if (in_xml->Parameters().present()) {
        xsd::cxx::tree::optional<bas::ppr::genericKeyValueParameters>& gp(in_xml->Parameters());
        ::xsd::cxx::tree::sequence<::bas::ppr::genericKVParam> myParams=(gp->Parameter());
        //::bas::ppr::genericKVParam name = ::bas::ppr::genericKVParam(myParams.Key());
     for(auto i = myParams.begin(); i!= myParams.end(); ++i){
       ::bas::ppr::genericKVParam& ii(*i);
       std::string name = ii.Key();

       if (name.compare("MethodType") == 0 && ii.StringValue().present()){
         methodType = ii.StringValue().get();
       }
       logger.info()<< "MethodType: " << methodType;
     }
     }
   }
  } catch (const xml_schema::exception& e) {
   std::cerr << e << std::endl;
  }
  logger.debug() << "Found FITS catalog filename " << catalog_file << " in input XML";

  logger.info() << "Finished getting information from file " << in_xml_catalog_file;

  return DmInput(catalog_file, methodType);
 }

 DmInput DmInput::readKSBCatalogXMLFile (const boost::filesystem::path& in_xml_catalog_file) {
   return readCatalogXMLFile(Euclid::WeakLensing::TwoDMass::catalogType::KSBCatalog,
                             in_xml_catalog_file);
 }

 DmInput DmInput::readClusterCatalogXMLFile(const fs::path& in_xml_catalog) {
   return readCatalogXMLFile(Euclid::WeakLensing::TwoDMass::catalogType::clusterCatalog,
                             in_xml_catalog);
 }

 DmInput DmInput::readLensMCCatalogXMLFile (const boost::filesystem::path& in_xml_catalog_file) {
   return readCatalogXMLFile(Euclid::WeakLensing::TwoDMass::catalogType::LensMCCatalog,
                             in_xml_catalog_file);
 }

 DmInput DmInput::readMomentsMLCatalogXMLFile (const boost::filesystem::path& in_xml_catalog_file) {
   return readCatalogXMLFile(Euclid::WeakLensing::TwoDMass::catalogType::MomentsMLCatalog,
                             in_xml_catalog_file);
 }

 DmInput DmInput::readRegaussCatalogXMLFile (const boost::filesystem::path& in_xml_catalog_file) {
   return readCatalogXMLFile(Euclid::WeakLensing::TwoDMass::catalogType::RegaussCatalog,
                             in_xml_catalog_file);
 }

 DmInput DmInput::readVisibilityMaskXMLFile(const boost::filesystem::path& in_xml_catalog) {
   return readCatalogXMLFile(Euclid::WeakLensing::TwoDMass::catalogType::VisibilityMask,
                             in_xml_catalog);
 }

 DmInput DmInput::readLE2CatalogXMLFile(const boost::filesystem::path& in_xml_catalog) {
   return readCatalogXMLFile(Euclid::WeakLensing::TwoDMass::catalogType::LE2Catalog,
                             in_xml_catalog);
 }

 DmInput::DmInput(const fs::path& catalog_file, std::string& methodType)
      : m_catalog_file{catalog_file}, m_methodType{methodType} {
 }

 fs::path DmInput::getFitsCatalogFilename() const {
  return m_catalog_file;
 }

  std::string DmInput::getMethodType() const{
   return m_methodType;
  }
/*
 std::vector< std::pair <double, double> > DmInput::getVertices() const {
   return m_vertices;
 }*/

} /* namespace TwoDMass */
} /* namespace WeakLensing */
} /* namespace Euclid */
