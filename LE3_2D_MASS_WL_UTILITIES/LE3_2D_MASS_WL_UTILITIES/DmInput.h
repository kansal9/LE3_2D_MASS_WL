/**
 * @file LE3_2D_MASS_WL_UTILITIES/DmInput.h
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

#ifndef _DMINPUT_H
#define _DMINPUT_H

#include <string>
#include <boost/filesystem.hpp>
#include "ElementsKernel/Logging.h"
#include <utility>
#include "LE3_2D_MASS_WL_UTILITIES/Utils.h"

// Datamodel for INPUT products
//#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ShearCatalogs.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-RegaussCatalog.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-MomentsMLCatalog.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-LensMCCatalog.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-KSBCatalog.h"
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-ClusterCatalog.h"
// Input catalogue for all WL PFs
#include "ST_DataModelBindings/dpd/le3/wl/inp/euc-test-le3-wl-InputLE2Catalog.h"
#include "ST_DataModelBindings/dpd/le3/wl/inp/euc-test-le3-wl-KSBCatalog.h"
#include "ST_DataModelBindings/dpd/le3/wl/inp/euc-test-le3-wl-LensMCCatalog.h"
#include "ST_DataModelBindings/dpd/le3/wl/inp/euc-test-le3-wl-MomentsMLCatalog.h"
#include "ST_DataModelBindings/dpd/le3/wl/inp/euc-test-le3-wl-RegaussCatalog.h"
// Input Visibility Mask (TODO Need to be checked exact location wrt new Data Model)
#include "ST_DataModelBindings/dpd/le3/wl/twodmass/inp/euc-test-le3-wl-twodmass-VisibilityMask.h"
// Used to read generic Key/Value parameters that are not in Data Model (give some flexibility)
// Read the (optional) parameter "MethodType = LENSMC/KSB/MOMENTSML/REGAUSS" not in Data Model
#include "ST_DataModelBindings/dictionary/bas/ppr/euc-test-ppr.h"

namespace Euclid {
 namespace WeakLensing {
  namespace TwoDMass {

/**
 * @class DmInput
 * @brief This class reads the name of  input catalog file from input xml file
 *
 */
class DmInput {

public:

 // DmInput();

  /**
   * @brief Destructor
   */
  virtual ~DmInput() = default;

 /**
  * @brief     Read Catalog Filename and parameters from input KSB shear catalog XML file
  * @param     input XML filename, <filesystem::path> path and name of the file to parse
  * @return    <filesystem::path> Shear Catalog File name and path
 */
  static DmInput readKSBCatalogXMLFile (const boost::filesystem::path& in_xml_catalog_file);

 /**
  * @brief     Read Catalog Filename and parameters from input LensMC shear catalog XML file
  * @param     input XML filename, <filesystem::path> path and name of the file to parse
  * @return    <filesystem::path> Shear Catalog File name and path
 */
  static DmInput readLensMCCatalogXMLFile (const boost::filesystem::path& in_xml_catalog_file);

 /**
  * @brief     Read Catalog Filename and parameters from input MomentsML shear catalog XML file
  * @param     input XML filename, <filesystem::path> path and name of the file to parse
  * @return    <filesystem::path> Shear Catalog File name and path
 */
  static DmInput readMomentsMLCatalogXMLFile (const boost::filesystem::path& in_xml_catalog_file);

 /**
  * @brief     Read Catalog Filename and parameters from input Regauss shear catalog XML file
  * @param     input XML filename, <filesystem::path> path and name of the file to parse
  * @return    <filesystem::path> Shear Catalog File name and path
 */
  static DmInput readRegaussCatalogXMLFile (const boost::filesystem::path& in_xml_catalog_file);

 /**
  * @brief     Read Catalog Filename and parameters from input cluster catalog XML file
  * @param     input XML filename, <filesystem::path> path and name of the file to parse
  * @return    <filesystem::path> cluster Catalog File name and path
 */
  static DmInput readClusterCatalogXMLFile(const boost::filesystem::path& in_xml_catalog);

 /**
  * @brief     Read Catalog Filename and parameters from input LE2 catalog XML file
  * @param     input XML filename, <filesystem::path> path and name of the file to parse
  * @return    <filesystem::path> LE2 Catalog File name and path
 */
  static DmInput readLE2CatalogXMLFile(const boost::filesystem::path& in_xml_catalog);

 /**
  * @brief     Read Catalog Filename and parameters from input Proxy Shear catalog XML file
  * @param     input XML filename, <filesystem::path> path and name of the file to parse
  * @return    <filesystem::path> LE2 Catalog File name and path
 */
  //static DmInput readProxyShearCatalogXMLFile(const boost::filesystem::path& in_xml_catalog);

 /**
  * @brief     Read Filename from input Visibility Mask (XML file)
  * @param     input XML filename, <filesystem::path> path and name of the file to parse
  * @return    <filesystem::path> Visibility Mask File name and path
 */
  static DmInput readVisibilityMaskXMLFile(const boost::filesystem::path& in_xml_catalog);

 /**
  * @brief     gets Catalog Filename in Fits format
  * @return    <filesystem::path> Catalog File name
 */
  boost::filesystem::path getFitsCatalogFilename() const;

  //std::vector< std::pair <double, double> > getVertices() const;

 /**
  * @brief     returns type of method used to estimate shear
 */
  std::string getMethodType() const;

private:

 /**
  * @brief     Read Catalog Filename and parameters from input shear catalog XML file
  * @param     input XML filename, <filesystem::path> path and name of the file to parse
  * @return    <filesystem::path> Shear Catalog File name and path
 */
  static DmInput readCatalogXMLFile (const Euclid::WeakLensing::TwoDMass::catalogType catType,
                                     const boost::filesystem::path& in_xml_catalog_file);
  DmInput(const boost::filesystem::path& catalog_file, std::string& methodType);

  boost::filesystem::path m_catalog_file;

  std::string m_methodType;

  //static std::vector< std::pair <double, double> > m_vertices;

}; /* End of DmInput class */
} /* namespace TwoDMass */
} /* namespace WeakLensing */
} /* namespace Euclid */

#endif /* DMINPUT_H */
